#!/usr/bin/env python3
"""Run Cellpose segmentation for preprocessed IMC images.

This script is intended to be called from the workflow notebook.

Input images are read from:

    data/img/

Cellpose masks are written to:

    data/masks/

A numbered mask inventory and log are written to:

    results/
    logs/

The script uses panel.csv to identify nuclear and membrane markers.

Composite construction (review item 2.2): `--normalization percentile`
(default) is the predecessor behaviour: raw pixel-wise max, then per-image
1st-99.5th percentile stretch. `--normalization fixed` divides each channel by
a dataset-wide scale (`--scale-file`, column display_max) and clips to 0-1
before aggregation. Fixed scaling is kept for the segmentation benchmark
(review item 4.1) but is not the default: on 2026-10-09 it reduced ROI032_I13
from 615 to 23-26 cells and ROI107_A5 from 127 to 14-16, because channels with
a scale of a few counts (CD3 3, CD66b 4, SOX9 8) turn single-count noise into
25-33 % of full scale and bury the real DNA/S100 signal.

Reproducibility (review item 4.4): the log records the Cellpose model, package
versions and host; ROIs whose mask already exists are skipped unless
`--overwrite` is given, and a failing ROI is logged and does not stop the run.
"""

from __future__ import annotations

import argparse
import csv
import platform
import socket
import warnings
from datetime import datetime
from pathlib import Path

import numpy as np
import pandas as pd
import tifffile
from cellpose import models


warnings.filterwarnings(
    "ignore",
    message=r"Sparse invariant checks are implicitly disabled\..*",
    category=UserWarning,
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Segment IMC images with Cellpose."
    )

    parser.add_argument(
        "--workflow-dir",
        type=Path,
        default=Path(__file__).resolve().parents[1],
        help="Workflow directory containing the data folder.",
    )

    parser.add_argument(
        "--numbered-inventory",
        type=Path,
        default=Path("results/04_cellpose_mask_inventory.csv"),
        help="Numbered inventory of generated Cellpose mask files.",
    )

    parser.add_argument(
        "--log-file",
        type=Path,
        default=Path("logs/03_segment_cellpose.log"),
        help="Cellpose segmentation log.",
    )

    parser.add_argument(
        "--gpu",
        action="store_true",
        help="Use GPU if available.",
    )

    parser.add_argument(
        "--flow-threshold",
        type=float,
        default=0.4,
        help="Cellpose flow threshold.",
    )

    parser.add_argument(
        "--cellprob-threshold",
        type=float,
        default=0.0,
        help="Cellpose cell probability threshold.",
    )

    parser.add_argument(
        "--min-size",
        type=int,
        default=20,
        help="Minimum segmented object size in pixels.",
    )

    parser.add_argument(
        "--model",
        default="cpsam_v2",
        help="Cellpose pretrained model (cpsam_v2 is the Cellpose 4.2 default).",
    )

    parser.add_argument(
        "--normalization",
        choices=["fixed", "percentile"],
        default="percentile",
        help="percentile: predecessor per-image stretch (default); fixed: per-channel dataset-wide scale.",
    )

    parser.add_argument(
        "--scale-file",
        type=Path,
        default=Path("data/reference/QuPath_display_ranges.csv"),
        help="CSV with columns 'marker' and 'display_max' giving the fixed scale per channel.",
    )

    parser.add_argument(
        "--nuclear-aggregation",
        choices=["max", "mean"],
        default="max",
        help="How nuclear channels are combined with --normalization fixed (percentile always uses max).",
    )

    parser.add_argument(
        "--membrane-aggregation",
        choices=["max", "mean"],
        default="max",
        help="How membrane channels are combined with --normalization fixed (percentile always uses max).",
    )

    parser.add_argument(
        "--overwrite",
        action="store_true",
        help="Re-segment ROIs whose mask already exists.",
    )

    return parser.parse_args()


def require_input(path: Path, description: str) -> None:
    if not path.exists():
        raise FileNotFoundError(
            f"Required {description} not found: {path}"
        )


def normalize_channel(image: np.ndarray) -> np.ndarray:
    """Robustly normalize an image to approximately 0-1."""

    image = image.astype(np.float32)

    low = np.percentile(image, 1)
    high = np.percentile(image, 99.5)

    if high <= low:
        return np.zeros_like(image, dtype=np.float32)

    image = (image - low) / (high - low)
    image = np.clip(image, 0, 1)

    return image


def get_segmentation_indices(panel_df: pd.DataFrame):
    """Identify nuclear and membrane channels for Cellpose segmentation.

    The panel.csv 'cellpose' column uses:
        1 = nuclear channel
        2 = membrane/cell-boundary channel
    """

    if "cellpose" not in panel_df.columns:
        raise ValueError(
            "panel.csv must contain a 'cellpose' column with "
            "nuclear=1 and membrane/cell-boundary=2 assignments."
        )

    cellpose_values = pd.to_numeric(
        panel_df["cellpose"],
        errors="coerce",
    )

    nuclear_indices = panel_df.index[
        cellpose_values == 1
    ].tolist()

    membrane_indices = panel_df.index[
        cellpose_values == 2
    ].tolist()

    if not nuclear_indices:
        raise ValueError(
            "No nuclear channels found in panel.csv with cellpose=1."
        )

    if not membrane_indices:
        raise ValueError(
            "No membrane/cell-boundary channels found in panel.csv "
            "with cellpose=2."
        )

    return nuclear_indices, membrane_indices


def load_fixed_scales(scale_file: Path) -> dict[str, float]:
    require_input(scale_file, "fixed channel scale file")
    scales = pd.read_csv(scale_file)
    return dict(zip(scales["marker"], scales["display_max"].astype(float)))


def build_composite(
    stack: np.ndarray,
    indices: list[int],
    names: list[str],
    normalization: str,
    aggregation: str,
    scales: dict[str, float],
) -> np.ndarray:
    aggregate = np.max if aggregation == "max" else np.mean

    if normalization == "percentile":
        # Predecessor behaviour: raw pixel-wise max, then per-image stretch.
        return normalize_channel(np.max(stack[indices], axis=0))

    scaled = [
        np.clip(stack[i].astype(np.float32) / scales[name], 0, 1)
        for i, name in zip(indices, names)
    ]
    return aggregate(np.stack(scaled), axis=0).astype(np.float32)


def run_cellpose(
    workflow_dir: Path,
    log_file: Path,
    gpu: bool,
    flow_threshold: float,
    cellprob_threshold: float,
    min_size: int,
    model_name: str,
    normalization: str,
    scale_file: Path,
    nuclear_aggregation: str,
    membrane_aggregation: str,
    overwrite: bool,
) -> list[str]:
    import cellpose
    import torch

    data_dir = workflow_dir / "data"
    images_dir = data_dir / "img"
    masks_dir = data_dir / "masks"
    masks_dir.mkdir(parents=True, exist_ok=True)

    panel_df = pd.read_csv(data_dir / "panel.csv")
    nuclear_indices, membrane_indices = get_segmentation_indices(panel_df)
    nuclear_names = [panel_df.loc[i, "name"] for i in nuclear_indices]
    membrane_names = [panel_df.loc[i, "name"] for i in membrane_indices]

    scales = {}
    if normalization == "fixed":
        scales = load_fixed_scales(scale_file)
        missing = [n for n in nuclear_names + membrane_names if n not in scales]
        if missing:
            raise ValueError(f"No fixed scale in {scale_file} for: {missing}")

    image_paths = sorted(images_dir.glob("*.tiff"))
    if not image_paths:
        raise FileNotFoundError(f"No TIFF images found in: {images_dir}")

    model = models.CellposeModel(gpu=gpu, pretrained_model=model_name)

    log_file.parent.mkdir(parents=True, exist_ok=True)
    errors: list[str] = []

    def describe(names: list[str]) -> str:
        return ", ".join(
            f"{n} (scale {scales[n]:g})" if normalization == "fixed" else n
            for n in names
        )

    with log_file.open("a", encoding="utf-8") as log:
        def write(line: str = "") -> None:
            log.write(line + "\n")
            log.flush()

        write("CELLPOSE SEGMENTATION")
        write(f"run_started: {datetime.now().isoformat(timespec='seconds')}")
        write(f"host: {socket.gethostname()}")
        write(f"python: {platform.python_version()}")
        write(f"cellpose: {cellpose.version}")
        write(f"torch: {torch.__version__}")
        write(f"numpy: {np.__version__}")
        write(f"tifffile: {tifffile.__version__}")
        write(f"model_requested: {model_name}")
        write(f"model_loaded: {getattr(model, 'pretrained_model', 'unknown')}")
        write(f"gpu: {gpu}")
        write(f"flow_threshold: {flow_threshold}")
        write(f"cellprob_threshold: {cellprob_threshold}")
        write(f"min_size: {min_size}")
        write(f"diameter: none (model default)")
        write(f"normalization: {normalization}")
        if normalization == "fixed":
            write(f"scale_file: {scale_file}")
        write(f"nuclear_aggregation: {nuclear_aggregation if normalization == 'fixed' else 'max'}")
        write(f"membrane_aggregation: {membrane_aggregation if normalization == 'fixed' else 'max'}")
        write(f"nuclear_channels: {describe(nuclear_names)}")
        write(f"membrane_channels: {describe(membrane_names)}")
        write(f"images_found: {len(image_paths)}")
        write(f"overwrite_existing_masks: {overwrite}")
        write()

        for image_index, image_path in enumerate(image_paths, start=1):
            output_path = masks_dir / image_path.name
            prefix = f"[{image_index}/{len(image_paths)}] {image_path.name}"

            if output_path.exists() and not overwrite:
                print(f"{prefix}: mask exists, skipped")
                write(f"{image_path.name}\tskipped_existing")
                continue

            try:
                stack = tifffile.imread(image_path)
                if stack.ndim != 3:
                    raise ValueError(f"expected C x H x W, got {stack.shape}")
                if max(nuclear_indices + membrane_indices) >= stack.shape[0]:
                    raise ValueError(
                        f"panel.csv channel indices exceed the {stack.shape[0]} image channels"
                    )

                nuclear = build_composite(
                    stack, nuclear_indices, nuclear_names,
                    normalization, nuclear_aggregation, scales,
                )
                membrane = build_composite(
                    stack, membrane_indices, membrane_names,
                    normalization, membrane_aggregation, scales,
                )

                # Cellpose input: channel 0 = membrane, channel 1 = nuclear
                result = model.eval(
                    np.stack([membrane, nuclear], axis=-1),
                    channel_axis=-1,
                    flow_threshold=flow_threshold,
                    cellprob_threshold=cellprob_threshold,
                    min_size=min_size,
                )
                masks = result[0].astype(np.uint32)
                tifffile.imwrite(output_path, masks)

                cell_count = int(masks.max())
                print(f"{prefix}: {cell_count} cells")
                write(f"{image_path.name}\tcells={cell_count}\tshape={masks.shape}")

            except Exception as error:  # one bad ROI must not stop the batch
                print(f"{prefix}: ERROR {error!r}")
                write(f"{image_path.name}\tERROR\t{error!r}")
                errors.append(image_path.name)

        write(f"run_finished: {datetime.now().isoformat(timespec='seconds')}; errors={len(errors)}")
        write()

    return errors


def write_mask_inventory(
    workflow_dir: Path,
    numbered_inventory: Path,
) -> None:

    masks_dir = (
        workflow_dir
        / "data"
        / "masks"
    )

    require_input(
        masks_dir,
        "Cellpose masks output directory",
    )

    rows = []

    for index, path in enumerate(
        sorted(masks_dir.glob("*.tiff")),
        start=1,
    ):

        mask = tifffile.imread(
            path
        )

        rows.append(
            {
                "index": index,
                "file_name": path.name,
                "relative_path": str(
                    path.relative_to(
                        workflow_dir
                    )
                ),
                "size_bytes": path.stat().st_size,
                "cell_count": int(
                    mask.max()
                ),
                "height": int(
                    mask.shape[0]
                ),
                "width": int(
                    mask.shape[1]
                ),
            }
        )

    destination = (
        workflow_dir
        / numbered_inventory
    )

    destination.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    with destination.open(
        "w",
        newline="",
        encoding="utf-8",
    ) as handle:

        fieldnames = [
            "index",
            "file_name",
            "relative_path",
            "size_bytes",
            "cell_count",
            "height",
            "width",
        ]

        writer = csv.DictWriter(
            handle,
            fieldnames=fieldnames,
        )

        writer.writeheader()
        writer.writerows(
            rows
        )


def main() -> None:

    args = parse_args()

    workflow_dir = (
        args.workflow_dir.resolve()
    )

    require_input(
        workflow_dir / "data" / "img",
        "preprocessed image directory",
    )

    require_input(
        workflow_dir / "data" / "images.csv",
        "Steinbock images.csv file",
    )

    require_input(
        workflow_dir / "data" / "panel.csv",
        "Steinbock panel file",
    )

    log_file = (
        workflow_dir
        / args.log_file
    )

    scale_file = (
        args.scale_file
        if args.scale_file.is_absolute()
        else workflow_dir / args.scale_file
    )

    errors = run_cellpose(
        workflow_dir=workflow_dir,
        log_file=log_file,
        gpu=args.gpu,
        flow_threshold=args.flow_threshold,
        cellprob_threshold=args.cellprob_threshold,
        min_size=args.min_size,
        model_name=args.model,
        normalization=args.normalization,
        scale_file=scale_file,
        nuclear_aggregation=args.nuclear_aggregation,
        membrane_aggregation=args.membrane_aggregation,
        overwrite=args.overwrite,
    )

    write_mask_inventory(
        workflow_dir,
        args.numbered_inventory,
    )

    print(
        f"Cellpose segmentation complete: "
        f"{workflow_dir / 'data' / 'masks'}"
    )

    print(
        f"Numbered mask inventory: "
        f"{workflow_dir / args.numbered_inventory}"
    )

    print(
        f"Segmentation log: "
        f"{log_file}"
    )

    if errors:
        raise SystemExit(
            f"Segmentation failed for {len(errors)} ROI(s): {errors}. See {log_file}"
        )


if __name__ == "__main__":
    main()
