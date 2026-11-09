#!/usr/bin/env python3
"""Run Steinbock per-cell region property measurement.

This script is intended to be called from the workflow notebook after image
preprocessing and cell segmentation have completed.

Steinbock command:

    steinbock measure regionprops

Steinbock working outputs are written to `data/regionprops/`. Numbered research
outputs are written to `results/` and `logs/`.
"""

from __future__ import annotations

import argparse
import csv
import subprocess
from pathlib import Path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Measure per-cell region properties with Steinbock.")
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
        "--numbered-inventory",
        type=Path,
        default=Path("results/06_regionprops_table_inventory.csv"),
        help="Numbered inventory of generated regionprops CSV files.",
    )
    parser.add_argument(
        "--log-file",
        type=Path,
        default=Path("logs/05_measure_regionprops.log"),
        help="Numbered regionprops command log.",
    )
    return parser.parse_args()


def require_input(path: Path, description: str) -> None:
    if not path.exists():
        raise FileNotFoundError(f"Required {description} not found: {path}")


def run_steinbock_measure_regionprops(workflow_dir: Path, image: str, log_file: Path) -> None:
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
        "regionprops",
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
        raise RuntimeError(f"Steinbock regionprops measurement failed. See log: {log_file}")


def write_regionprops_inventory(workflow_dir: Path, numbered_inventory: Path) -> None:
    regionprops_dir = workflow_dir / "data" / "regionprops"
    require_input(regionprops_dir, "Steinbock regionprops output directory")

    rows = []
    for index, path in enumerate(sorted(regionprops_dir.glob("*.csv")), start=1):
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

    run_steinbock_measure_regionprops(workflow_dir, args.steinbock_image, workflow_dir / args.log_file)
    write_regionprops_inventory(workflow_dir, args.numbered_inventory)

    print(f"Region property measurement complete: {workflow_dir / 'data' / 'regionprops'}")
    print(f"Numbered regionprops inventory: {workflow_dir / args.numbered_inventory}")
    print(f"Command log: {workflow_dir / args.log_file}")


if __name__ == "__main__":
    main()
