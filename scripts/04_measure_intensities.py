#!/usr/bin/env python3
"""Run Steinbock per-cell marker intensity measurement.

This script is intended to be called from the workflow notebook after image
preprocessing and cell segmentation have completed.

Steinbock command:

    steinbock measure intensities

Steinbock working outputs are written to `data/intensities/`. Numbered research
outputs are written to `results/` and `logs/`.
"""

from __future__ import annotations

import argparse
import csv
import subprocess
from pathlib import Path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Measure per-cell marker intensities with Steinbock.")
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
        "--docker-platform",
        default="linux/amd64",
        help="Docker platform passed to `docker run --platform` for Steinbock compatibility on Apple Silicon.",
    )
    parser.add_argument(
        "--numbered-inventory",
        type=Path,
        default=Path("results/05_intensity_table_inventory.csv"),
        help="Numbered inventory of generated intensity CSV files.",
    )
    parser.add_argument(
        "--log-file",
        type=Path,
        default=Path("logs/04_measure_intensities.log"),
        help="Numbered intensity measurement command log.",
    )
    return parser.parse_args()


def require_input(path: Path, description: str) -> None:
    if not path.exists():
        raise FileNotFoundError(f"Required {description} not found: {path}")


def run_steinbock_measure_intensities(
    workflow_dir: Path,
    image: str,
    docker_platform: str,
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
        "--platform",
        docker_platform,
        "-v",
        f"{data_dir.resolve()}:/data",
        "-u",
        user_group,
        image,
        "measure",
        "intensities",
        "--mmap",
    ]

    log_file.parent.mkdir(parents=True, exist_ok=True)
    log_file.write_text(
        "WORKFLOW_DIR\n"
        + str(workflow_dir)
        + "\n\nDATA_DIR\n"
        + str(data_dir.resolve())
        + "\n\nCOMMAND\n"
        + " ".join(command)
        + "\n\nSTATUS\nrunning\n",
        encoding="utf-8",
    )
    completed = subprocess.run(command, cwd=workflow_dir, text=True, capture_output=True)

    log_file.write_text(
        "WORKFLOW_DIR\n"
        + str(workflow_dir)
        + "\n\nDATA_DIR\n"
        + str(data_dir.resolve())
        + "\n\nCOMMAND\n"
        + " ".join(command)
        + "\n\nSTDOUT\n"
        + completed.stdout
        + "\n\nSTDERR\n"
        + completed.stderr
        + f"\n\nRETURN_CODE\n{completed.returncode}\n",
        encoding="utf-8",
    )

    if completed.returncode != 0:
        raise RuntimeError(f"Steinbock intensity measurement failed. See log: {log_file}")


def write_intensity_inventory(workflow_dir: Path, numbered_inventory: Path) -> None:
    intensities_dir = workflow_dir / "data" / "intensities"
    require_input(intensities_dir, "Steinbock intensities output directory")

    rows = []
    for index, path in enumerate(sorted(intensities_dir.glob("*.csv")), start=1):
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

    require_input(workflow_dir / "data" / "img", "preprocessed image directory")
    require_input(workflow_dir / "data" / "masks", "segmentation mask directory")
    require_input(workflow_dir / "data" / "panel.csv", "Steinbock panel file")

    run_steinbock_measure_intensities(
        workflow_dir,
        args.steinbock_image,
        args.docker_platform,
        workflow_dir / args.log_file,
    )
    write_intensity_inventory(workflow_dir, args.numbered_inventory)

    print(f"Intensity measurement complete: {workflow_dir / 'data' / 'intensities'}")
    print(f"Numbered intensity inventory: {workflow_dir / args.numbered_inventory}")
    print(f"Command log: {workflow_dir / args.log_file}")


if __name__ == "__main__":
    main()
