#!/usr/bin/env python3
"""Run Steinbock spatial neighbor measurement.

This script is intended to be called from the workflow notebook after cell
segmentation has completed.

Steinbock command:

    steinbock measure neighbors --type expansion --dmax 4

Steinbock working outputs are written to `data/neighbors/`. Numbered research
outputs are written to `results/` and `logs/`.
"""

from __future__ import annotations

import argparse
import csv
import subprocess
from pathlib import Path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Measure spatial cell neighbors with Steinbock.")
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
        "--neighbor-type",
        default="expansion",
        help="Steinbock neighbor type.",
    )
    parser.add_argument(
        "--dmax",
        default="4",
        help="Maximum expansion distance in pixels.",
    )
    parser.add_argument(
        "--numbered-inventory",
        type=Path,
        default=Path("results/07_neighbors_table_inventory.csv"),
        help="Numbered inventory of generated neighbor CSV files.",
    )
    parser.add_argument(
        "--log-file",
        type=Path,
        default=Path("logs/06_measure_neighbors.log"),
        help="Numbered neighbors command log.",
    )
    return parser.parse_args()


def require_input(path: Path, description: str) -> None:
    if not path.exists():
        raise FileNotFoundError(f"Required {description} not found: {path}")


def run_steinbock_measure_neighbors(
    workflow_dir: Path,
    image: str,
    neighbor_type: str,
    dmax: str,
    log_file: Path,
) -> None:
    data_dir = workflow_dir / "data"
    user_group = (
        f"{subprocess.check_output(['id', '-u'], text=True).strip()}:"
        f"{subprocess.check_output(['id', '-g'], text=True).strip()}"
    )

    command = [
        "docker",
        "run",
        "--rm",
        "-v",
        f"{data_dir.resolve()}:/data",
        "-u",
        user_group,
        image,
        "measure",
        "neighbors",
        "--type",
        neighbor_type,
        "--dmax",
        dmax,
    ]

    log_file.parent.mkdir(parents=True, exist_ok=True)
    completed = subprocess.run(command, cwd=workflow_dir, text=True, capture_output=True)

    log_file.write_text(
        "COMMAND\n"
        + " ".join(command)
        + "\n\nSTDOUT\n"
        + completed.stdout
        + "\n\nSTDERR\n"
        + completed.stderr
        + f"\n\nRETURN_CODE\n{completed.returncode}\n",
        encoding="utf-8",
    )

    if completed.returncode != 0:
        raise RuntimeError(f"Steinbock neighbor measurement failed. See log: {log_file}")


def write_neighbors_inventory(workflow_dir: Path, numbered_inventory: Path) -> None:
    neighbors_dir = workflow_dir / "data" / "neighbors"
    require_input(neighbors_dir, "Steinbock neighbors output directory")

    rows = []
    for index, path in enumerate(sorted(neighbors_dir.glob("*.csv")), start=1):
        rows.append(
            {
                "index": index,
                "file_name": path.name,
                "relative_path": str(path.relative_to(workflow_dir)),
                "size_bytes": path.stat().st_size,
            }
        )

    destination = workflow_dir / numbered_inventory
    destination.parent.mkdir(parents=True, exist_ok=True)
    with destination.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=["index", "file_name", "relative_path", "size_bytes"])
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    args = parse_args()
    workflow_dir = args.workflow_dir.resolve()

    require_input(workflow_dir / "data" / "masks", "segmentation mask directory")

    run_steinbock_measure_neighbors(
        workflow_dir,
        args.steinbock_image,
        args.neighbor_type,
        args.dmax,
        workflow_dir / args.log_file,
    )
    write_neighbors_inventory(workflow_dir, args.numbered_inventory)

    print(f"Neighbor measurement complete: {workflow_dir / 'data' / 'neighbors'}")
    print(f"Numbered neighbor inventory: {workflow_dir / args.numbered_inventory}")
    print(f"Command log: {workflow_dir / args.log_file}")


if __name__ == "__main__":
    main()
