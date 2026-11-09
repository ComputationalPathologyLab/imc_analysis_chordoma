#!/usr/bin/env python3
"""Export Steinbock single-cell and spatial graph outputs.

This script is intended to be called from the workflow notebook after intensity,
region property, and neighbor measurements have completed.

Steinbock commands:

    steinbock export csv intensities regionprops -o cells.csv
    steinbock export anndata --intensities intensities --data regionprops --neighbors neighbors -o cells.h5ad
    steinbock export graphs --format graphml --data intensities --data regionprops

Steinbock working outputs are written to `data/`. Numbered research copies are
written to `results/`.

After the export, the TMA key (data/roi_key_map.csv, see scripts/apply_tma_key.py)
is added to every cell and graph: roi_number, key_label (patient + core, e.g. E6),
patient_id, core, sample_type, mcd_label, acquired_position and repeat_of.
"""

from __future__ import annotations

import argparse
import csv
import shutil
import subprocess
from pathlib import Path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Export Steinbock analysis-ready data.")
    parser.add_argument(
        "--workflow-dir",
        type=Path,
        default=Path(__file__).resolve().parents[1],
        help="Workflow directory containing the data folder.",
    )
    parser.add_argument(
        "--steinbock-image",
        default="ghcr.io/bodenmillergroup/steinbock:0.16.1",
        help="Steinbock Docker image.",
    )
    parser.add_argument(
        "--numbered-cells-csv",
        type=Path,
        default=Path("results/08_cells.csv"),
        help="Numbered copy of exported cells.csv.",
    )
    parser.add_argument(
        "--numbered-cells-h5ad",
        type=Path,
        default=Path("results/09_cells.h5ad"),
        help="Numbered copy of exported cells.h5ad.",
    )
    parser.add_argument(
        "--numbered-graphs-dir",
        type=Path,
        default=Path("results/10_graphml"),
        help="Numbered directory containing copies of exported GraphML files.",
    )
    parser.add_argument(
        "--numbered-graph-inventory",
        type=Path,
        default=Path("results/10_graphml_inventory.csv"),
        help="Numbered inventory of exported GraphML files.",
    )
    parser.add_argument(
        "--log-file",
        type=Path,
        default=Path("logs/07_export_data.log"),
        help="Numbered export command log.",
    )
    return parser.parse_args()


def require_input(path: Path, description: str) -> None:
    if not path.exists():
        raise FileNotFoundError(f"Required {description} not found: {path}")


def docker_command(workflow_dir: Path, image: str, *steinbock_args: str) -> list[str]:
    data_dir = workflow_dir / "data"
    user_group = (
        f"{subprocess.check_output(['id', '-u'], text=True).strip()}:"
        f"{subprocess.check_output(['id', '-g'], text=True).strip()}"
    )
    return [
        "docker",
        "run",
        "--rm",
        "-v",
        f"{data_dir.resolve()}:/data",
        "-u",
        user_group,
        image,
        *steinbock_args,
    ]


def run_command(command: list[str], workflow_dir: Path) -> subprocess.CompletedProcess[str]:
    return subprocess.run(command, cwd=workflow_dir, text=True, capture_output=True)


def run_exports(workflow_dir: Path, image: str, log_file: Path) -> None:
    commands = [
        docker_command(workflow_dir, image, "export", "csv", "intensities", "regionprops", "-o", "cells.csv"),
        docker_command(
            workflow_dir,
            image,
            "export",
            "anndata",
            "--intensities",
            "intensities",
            "--data",
            "regionprops",
            "--neighbors",
            "neighbors",
            "-o",
            "cells.h5ad",
        ),
        docker_command(
            workflow_dir,
            image,
            "export",
            "graphs",
            "--format",
            "graphml",
            "--data",
            "intensities",
            "--data",
            "regionprops",
        ),
    ]

    log_file.parent.mkdir(parents=True, exist_ok=True)
    log_parts = []

    for command in commands:
        completed = run_command(command, workflow_dir)
        log_parts.append(
            "COMMAND\n"
            + " ".join(command)
            + "\n\nSTDOUT\n"
            + completed.stdout
            + "\n\nSTDERR\n"
            + completed.stderr
            + f"\n\nRETURN_CODE\n{completed.returncode}\n"
        )

        if completed.returncode != 0:
            log_file.write_text("\n\n---\n\n".join(log_parts), encoding="utf-8")
            raise RuntimeError(f"Steinbock export failed. See log: {log_file}")

    log_file.write_text("\n\n---\n\n".join(log_parts), encoding="utf-8")


KEY_COLUMNS = ["roi_number", "key_label", "patient_id", "core", "sample_type",
               "mcd_label", "acquired_position", "repeat_of"]


def add_tma_key(workflow_dir: Path) -> None:
    """Add the TMA key columns to cells.csv, cells.h5ad and the graphs in data/."""
    import re

    import anndata as ad
    import pandas as pd

    key_map = workflow_dir / "data" / "roi_key_map.csv"
    require_input(key_map, "TMA key map")
    info = pd.read_csv(key_map, dtype=str, keep_default_na=False).set_index("roi")

    cells_path = workflow_dir / "data" / "cells.csv"
    cells = pd.read_csv(cells_path, dtype=str, keep_default_na=False)
    roi = cells["Image"].str.removesuffix(".tiff")
    if not roi.isin(info.index).all():
        raise ValueError(f"images missing from {key_map}: {sorted(set(roi) - set(info.index))[:5]}")
    for col in KEY_COLUMNS:
        cells[col] = roi.map(info[col]).values
    cells.to_csv(cells_path, index=False)

    h5ad_path = workflow_dir / "data" / "cells.h5ad"
    adata = ad.read_h5ad(h5ad_path)
    roi = adata.obs["Image"].astype(str).str.removesuffix(".tiff")
    for col in KEY_COLUMNS:
        adata.obs[col] = roi.map(info[col]).values
    adata.write_h5ad(h5ad_path)

    for path in sorted((workflow_dir / "data" / "graphs").glob("*.graphml")):
        text = path.read_text(encoding="utf-8")
        if 'for="graph"' in text:
            continue
        row = info.loc[path.stem]
        attrs = {"roi": path.stem, **{col: row[col] for col in KEY_COLUMNS if col != "roi_number"}}
        keys = "".join(f'<key id="g_{k}" for="graph" attr.name="{k}" attr.type="string"/>\n' for k in attrs)
        data = "".join(f'<data key="g_{k}">{v}</data>' for k, v in attrs.items())
        text = text.replace("<graph ", keys + "<graph ", 1)
        text = re.sub(r"(<graph [^>]*>)", lambda match: match.group(1) + data, text, count=1)
        path.write_text(text, encoding="utf-8")


def copy_file(source: Path, destination: Path, description: str) -> None:
    require_input(source, description)
    destination.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(source, destination)


def copy_graphs_and_write_inventory(
    workflow_dir: Path,
    numbered_graphs_dir: Path,
    numbered_graph_inventory: Path,
) -> None:
    source_dir = workflow_dir / "data" / "graphs"
    require_input(source_dir, "Steinbock GraphML output directory")

    graph_destination_dir = workflow_dir / numbered_graphs_dir
    graph_destination_dir.mkdir(parents=True, exist_ok=True)

    rows = []
    for index, source in enumerate(sorted(source_dir.glob("*.graphml")), start=1):
        destination_name = f"10_{index:03d}_{source.name}"
        destination = graph_destination_dir / destination_name
        shutil.copy2(source, destination)
        rows.append(
            {
                "index": index,
                "source_file_name": source.name,
                "numbered_file_name": destination.name,
                "source_relative_path": str(source.relative_to(workflow_dir)),
                "numbered_relative_path": str(destination.relative_to(workflow_dir)),
                "size_bytes": destination.stat().st_size,
            }
        )

    inventory_destination = workflow_dir / numbered_graph_inventory
    inventory_destination.parent.mkdir(parents=True, exist_ok=True)
    with inventory_destination.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "index",
                "source_file_name",
                "numbered_file_name",
                "source_relative_path",
                "numbered_relative_path",
                "size_bytes",
            ],
        )
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    args = parse_args()
    workflow_dir = args.workflow_dir.resolve()

    require_input(workflow_dir / "data" / "intensities", "intensity table directory")
    require_input(workflow_dir / "data" / "regionprops", "regionprops table directory")
    require_input(workflow_dir / "data" / "neighbors", "neighbors table directory")

    run_exports(workflow_dir, args.steinbock_image, workflow_dir / args.log_file)
    add_tma_key(workflow_dir)
    copy_file(workflow_dir / "data" / "cells.csv", workflow_dir / args.numbered_cells_csv, "cells.csv export")
    copy_file(workflow_dir / "data" / "cells.h5ad", workflow_dir / args.numbered_cells_h5ad, "cells.h5ad export")
    copy_graphs_and_write_inventory(workflow_dir, args.numbered_graphs_dir, args.numbered_graph_inventory)

    print(f"Exported cells CSV: {workflow_dir / args.numbered_cells_csv}")
    print(f"Exported AnnData object: {workflow_dir / args.numbered_cells_h5ad}")
    print(f"Exported numbered graphs: {workflow_dir / args.numbered_graphs_dir}")
    print(f"Graph inventory: {workflow_dir / args.numbered_graph_inventory}")
    print(f"Command log: {workflow_dir / args.log_file}")


if __name__ == "__main__":
    main()
