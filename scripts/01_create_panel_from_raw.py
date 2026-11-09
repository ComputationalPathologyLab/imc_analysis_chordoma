#!/usr/bin/env python3
"""Create a Steinbock-compatible panel.csv from raw IMC inputs.

For marker-channel OME-TIFF folder inputs, the script derives channel and marker
names directly from files named like `142Nd_CD138.ome.tiff` inside each ROI
folder.

The panel is the union of marker-channel files detected across all ROI folders,
so the workflow can handle one ROI, many ROIs, CH-only inputs, CHS-only inputs,
and mixed CH/CHS subsets with partially overlapping panels.

For legacy MCD/text workflows, if no OME-TIFF ROI folders are present, it scans
ROI `.txt` files, extracts channel and marker names from the first header,
validates that all ROI files share the same channel order, and writes a dynamic
Steinbock panel.

Default Cellpose assignments for the chordoma (CH) panel are:

- Cellpose channel 1 (nuclear):
    DNA1, DNA2, Brachiury

- Cellpose channel 2 (membrane / cell-boundary):
    S100, CD68, aSMA, CD31, CD3, CD66b

They are adapted from the chondrosarcoma (CHS) workflow (review item 2.2,
decision V.G. 06/10/2026). The CH panel has no SOX9; Brachyury, the nuclear
tumour marker of chordoma, takes its place in the nuclear composite. The
marker is spelled "Brachiury" in the exported filenames and is kept as such.
The membrane set is unchanged from CHS. Both choices are provisional until
confirmed for CH.

`--roi` restricts the panel and availability table to selected ROI folders,
so a single ROI can be processed as a test run.

The `analysis` column (review item 2.3) is 0 for technical channels (iodine,
xenon, barium, lead, argon dimer) and 1 for antibody/DNA channels. `keep`
stays 1 for every channel so that column names stay aligned with the stacked
image channels; 138Ba and 208Pb remain available for artifact masking.

These assignments can be changed using command-line arguments when applying
the workflow to another IMC dataset.
"""

from __future__ import annotations

import argparse
import csv
import re
from pathlib import Path


METADATA_COLUMNS = {
    "Start_push",
    "End_push",
    "Pushes_duration",
    "X",
    "Y",
    "Z",
}

CHANNEL_PATTERN = re.compile(
    r"^(?P<marker>.+)\((?P<channel>.+?)(?:Di|Ir)\)$"
)

OME_TIFF_PATTERN = re.compile(
    r"^(?P<channel>[^_]+)_(?P<marker>.+)\.ome\.tiff$"
)


def parse_args() -> argparse.Namespace:

    parser = argparse.ArgumentParser(
        description=(
            "Create data/panel.csv from raw OME-TIFF ROI folders "
            "or IMC ROI text-file headers."
        )
    )

    parser.add_argument(
        "--workflow-dir",
        type=Path,
        default=Path(__file__).resolve().parents[1],
        help=(
            "Workflow directory containing data/, results/, "
            "scripts/, and notebooks/."
        ),
    )

    parser.add_argument(
        "--raw-dir",
        type=Path,
        default=Path("data/raw"),
        help=(
            "Directory containing ROI folders with per-channel "
            ".ome.tiff files, or legacy raw files."
        ),
    )

    parser.add_argument(
        "--output",
        type=Path,
        default=Path("data/panel.csv"),
        help="Output Steinbock panel CSV path.",
    )

    parser.add_argument(
        "--numbered-copy",
        type=Path,
        default=Path("results/01_panel.csv"),
        help=(
            "Chronologically numbered research copy. "
            "Use 'none' to skip creating a numbered copy."
        ),
    )

    parser.add_argument(
        "--availability-output",
        type=Path,
        default=Path("data/roi_marker_availability.csv"),
        help=(
            "ROI-by-marker availability table for heterogeneous "
            "panel inputs."
        ),
    )

    parser.add_argument(
        "--numbered-availability-copy",
        type=Path,
        default=Path(
            "results/01_roi_marker_availability.csv"
        ),
        help=(
            "Numbered copy of the ROI-by-marker availability table. "
            "Use 'none' to skip."
        ),
    )

    # ---------------------------------------------------------
    # Cellpose nuclear channels
    # ---------------------------------------------------------

    parser.add_argument(
        "--nuclear-marker",
        nargs="+",
        default=[
            "DNA1",
            "DNA2",
            "Brachiury",
        ],
        help=(
            "Marker(s) or channel label(s) assigned as "
            "Cellpose nuclear channel (cellpose=1)."
        ),
    )

    # ---------------------------------------------------------
    # Cellpose membrane / cell-boundary channels
    # ---------------------------------------------------------

    parser.add_argument(
        "--membrane-marker",
        nargs="+",
        default=[
            "S100",
            "CD68",
            "aSMA",
            "CD31",
            "CD3",
            "CD66b",
        ],
        help=(
            "Marker(s) or channel label(s) assigned as "
            "Cellpose membrane/cell-boundary channel "
            "(cellpose=2)."
        ),
    )

    parser.add_argument(
        "--technical-channel",
        nargs="+",
        default=[
            "127I",
            "131Xe",
            "134Xe",
            "138Ba",
            "208Pb",
            "80ArAr",
        ],
        help=(
            "Channel(s) or marker label(s) marked analysis=0 "
            "(technical, not biological markers)."
        ),
    )

    parser.add_argument(
        "--roi",
        nargs="+",
        default=None,
        help=(
            "ROI folder name(s) to include, e.g. ROI001_D13. "
            "Default: all ROI folders in --raw-dir."
        ),
    )

    return parser.parse_args()


def workflow_path(
    workflow_dir: Path,
    path: Path,
) -> Path:

    return (
        path
        if path.is_absolute()
        else workflow_dir / path
    )


def read_header(path: Path) -> list[str]:

    with path.open(
        "r",
        encoding="utf-8",
        errors="replace",
    ) as handle:

        return (
            handle
            .readline()
            .rstrip("\n")
            .split("\t")
        )


def extract_marker_channels(
    header: list[str],
) -> list[tuple[str, str]]:

    marker_channels: list[
        tuple[str, str]
    ] = []

    for column in header:

        if column in METADATA_COLUMNS:
            continue

        match = CHANNEL_PATTERN.match(
            column
        )

        if match is None:
            continue

        marker = (
            match
            .group("marker")
            .strip()
        )

        channel = (
            match
            .group("channel")
            .strip()
        )

        marker_channels.append(
            (
                channel,
                marker,
            )
        )

    if not marker_channels:

        raise ValueError(
            "No IMC marker/channel columns were "
            "detected in the ROI header."
        )

    return marker_channels


def find_roi_text_files(
    raw_dir: Path,
) -> list[Path]:

    roi_files = sorted(
        raw_dir.glob("*.txt")
    )

    if not roi_files:

        raise FileNotFoundError(
            f"No ROI .txt files were found in {raw_dir}"
        )

    return roi_files


def find_ome_tiff_roi_dirs(
    raw_dir: Path,
    selected: list[str] | None = None,
) -> list[Path]:

    roi_dirs = sorted(
        path
        for path in raw_dir.iterdir()
        if (
            path.is_dir()
            and any(
                path.glob(
                    "*.ome.tiff"
                )
            )
        )
    )

    if selected:
        missing = sorted(set(selected) - {path.name for path in roi_dirs})
        if missing:
            raise FileNotFoundError(
                f"Selected ROI folder(s) not found in {raw_dir}: {missing}"
            )
        roi_dirs = [path for path in roi_dirs if path.name in selected]

    return roi_dirs


def extract_marker_channels_from_roi_dir(
    roi_dir: Path,
) -> list[tuple[str, str]]:

    marker_channels: list[
        tuple[str, str]
    ] = []

    ignored_files: list[str] = []

    for path in sorted(
        roi_dir.glob("*.ome.tiff")
    ):

        match = OME_TIFF_PATTERN.match(
            path.name
        )

        if match is None:

            ignored_files.append(
                path.name
            )

            continue

        marker_channels.append(
            (
                match.group(
                    "channel"
                ),
                match.group(
                    "marker"
                ),
            )
        )

    if ignored_files:

        ignored = ", ".join(
            ignored_files
        )

        raise ValueError(
            f"OME-TIFF files in {roi_dir} do not match "
            f"'<channel>_<marker>.ome.tiff': {ignored}"
        )

    if not marker_channels:

        raise ValueError(
            f"No marker-channel OME-TIFF files "
            f"were detected in {roi_dir}"
        )

    return marker_channels


def collect_ome_tiff_roi_marker_channels(
    roi_dirs: list[Path],
) -> tuple[
    list[tuple[str, str]],
    dict[
        str,
        set[tuple[str, str]],
    ],
]:

    union: set[
        tuple[str, str]
    ] = set()

    roi_marker_channels: dict[
        str,
        set[tuple[str, str]],
    ] = {}

    for roi_dir in roi_dirs:

        marker_channels = set(
            extract_marker_channels_from_roi_dir(
                roi_dir
            )
        )

        roi_marker_channels[
            roi_dir.name
        ] = marker_channels

        union.update(
            marker_channels
        )

    return (
        sorted(union),
        roi_marker_channels,
    )


def validate_consistent_headers(
    roi_files: list[Path],
) -> list[tuple[str, str]]:

    reference = (
        extract_marker_channels(
            read_header(
                roi_files[0]
            )
        )
    )

    mismatches: list[str] = []

    for roi_file in roi_files[1:]:

        current = (
            extract_marker_channels(
                read_header(
                    roi_file
                )
            )
        )

        if current != reference:

            mismatches.append(
                str(
                    roi_file
                )
            )

    if mismatches:

        mismatch_list = "\n".join(
            mismatches
        )

        raise ValueError(
            "Not all ROI files share the same "
            "marker/channel order. "
            f"Mismatching files:\n{mismatch_list}"
        )

    return reference


# -------------------------------------------------------------
# Cellpose channel assignment
# -------------------------------------------------------------

def cellpose_value(
    marker: str,
    channel: str,
    nuclear_markers: list[str],
    membrane_markers: list[str],
) -> str:

    if (
        marker in nuclear_markers
        or channel in nuclear_markers
    ):
        return "1"

    if (
        marker in membrane_markers
        or channel in membrane_markers
    ):
        return "2"

    return ""


def write_panel(
    marker_channels: list[
        tuple[str, str]
    ],
    output: Path,
    nuclear_markers: list[str],
    membrane_markers: list[str],
    technical_channels: list[str],
) -> None:

    output.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    with output.open(
        "w",
        newline="",
        encoding="utf-8",
    ) as handle:

        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "channel",
                "name",
                "keep",
                "ilastik",
                "deepcell",
                "cellpose",
                "analysis",
            ],
        )

        writer.writeheader()

        for index, (
            channel,
            marker,
        ) in enumerate(
            marker_channels,
            start=1,
        ):

            writer.writerow(
                {
                    "channel": channel,
                    "name": marker,
                    "keep": "1",

                    # Channel order for ilastik
                    "ilastik": str(index),

                    # We are no longer using Mesmer
                    "deepcell": "",

                    # Cellpose channel definition
                    "cellpose": cellpose_value(
                        marker,
                        channel,
                        nuclear_markers,
                        membrane_markers,
                    ),

                    "analysis": (
                        "0"
                        if (
                            channel in technical_channels
                            or marker in technical_channels
                        )
                        else "1"
                    ),
                }
            )


def write_availability(
    roi_marker_channels: dict[
        str,
        set[tuple[str, str]],
    ],
    marker_channels: list[
        tuple[str, str]
    ],
    output: Path,
) -> None:

    output.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    fieldnames = [
        "roi",
        "channel",
        "marker",
        "channel_marker",
        "available",
    ]

    with output.open(
        "w",
        newline="",
        encoding="utf-8",
    ) as handle:

        writer = csv.DictWriter(
            handle,
            fieldnames=fieldnames,
        )

        writer.writeheader()

        for roi in sorted(
            roi_marker_channels
        ):

            roi_set = (
                roi_marker_channels[
                    roi
                ]
            )

            for (
                channel,
                marker,
            ) in marker_channels:

                writer.writerow(
                    {
                        "roi": roi,
                        "channel": channel,
                        "marker": marker,
                        "channel_marker": (
                            f"{channel}_{marker}"
                        ),
                        "available": (
                            "1"
                            if (
                                channel,
                                marker,
                            ) in roi_set
                            else "0"
                        ),
                    }
                )


def copy_text_file(
    source: Path,
    destination: Path | None,
) -> None:

    if destination is None:
        return

    destination.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    destination.write_text(
        source.read_text(
            encoding="utf-8"
        ),
        encoding="utf-8",
    )


def main() -> None:

    args = parse_args()

    workflow_dir = (
        args.workflow_dir.resolve()
    )

    raw_dir = workflow_path(
        workflow_dir,
        args.raw_dir,
    ).resolve()

    output = workflow_path(
        workflow_dir,
        args.output,
    ).resolve()

    numbered_copy = (
        None
        if str(
            args.numbered_copy
        ).casefold() == "none"
        else workflow_path(
            workflow_dir,
            args.numbered_copy,
        ).resolve()
    )

    availability_output = (
        workflow_path(
            workflow_dir,
            args.availability_output,
        ).resolve()
    )

    numbered_availability_copy = (
        None
        if str(
            args.numbered_availability_copy
        ).casefold() == "none"
        else workflow_path(
            workflow_dir,
            args.numbered_availability_copy,
        ).resolve()
    )

    # ---------------------------------------------------------
    # Detect input format
    # ---------------------------------------------------------

    roi_dirs = (
        find_ome_tiff_roi_dirs(
            raw_dir,
            args.roi,
        )
    )

    if roi_dirs:

        (
            marker_channels,
            roi_marker_channels,
        ) = (
            collect_ome_tiff_roi_marker_channels(
                roi_dirs
            )
        )

        write_panel(
            marker_channels,
            output,
            args.nuclear_marker,
            args.membrane_marker,
            args.technical_channel,
        )

        write_availability(
            roi_marker_channels,
            marker_channels,
            availability_output,
        )

        copy_text_file(
            availability_output,
            numbered_availability_copy,
        )

        input_description = (
            f"OME-TIFF ROI folders checked: "
            f"{len(roi_dirs)}"
        )

    else:

        roi_files = (
            find_roi_text_files(
                raw_dir
            )
        )

        marker_channels = (
            validate_consistent_headers(
                roi_files
            )
        )

        write_panel(
            marker_channels,
            output,
            args.nuclear_marker,
            args.membrane_marker,
            args.technical_channel,
        )

        write_availability(
            {},
            marker_channels,
            availability_output,
        )

        copy_text_file(
            availability_output,
            numbered_availability_copy,
        )

        input_description = (
            f"ROI text files checked: "
            f"{len(roi_files)}"
        )

    # ---------------------------------------------------------
    # Save research copies
    # ---------------------------------------------------------

    copy_text_file(
        output,
        numbered_copy,
    )

    # ---------------------------------------------------------
    # Summary
    # ---------------------------------------------------------

    print(
        input_description
    )

    print(
        f"Panel rows written: "
        f"{len(marker_channels)}"
    )

    print(
        f"Steinbock panel path: "
        f"{output}"
    )

    if numbered_copy is not None:

        print(
            f"Numbered research copy: "
            f"{numbered_copy}"
        )

    print(
        f"ROI marker availability table: "
        f"{availability_output}"
    )

    if (
        numbered_availability_copy
        is not None
    ):

        print(
            f"Numbered availability copy: "
            f"{numbered_availability_copy}"
        )

    print(
        "Cellpose nuclear markers/channels "
        f"(cellpose=1): {args.nuclear_marker}"
    )

    print(
        "Cellpose membrane markers/channels "
        f"(cellpose=2): {args.membrane_marker}"
    )

    print(
        "Technical channels "
        f"(analysis=0): {args.technical_channel}"
    )


if __name__ == "__main__":
    main()