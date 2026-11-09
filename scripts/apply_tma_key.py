#!/usr/bin/env python3
"""Rename all CH outputs from slide positions to the TMA key (patient letter + core number).

Old names carry the position printed on the slide, e.g. ROI013_F11 (row F,
column 11). The TMA key (data/reference/Chordoma - TMA Ref.pdf) labels the same
core by patient and core number: columns are patients O..A from left to right
(column 1 = O, column 15 = A) and rows A..K are cores 1..11. ROI013_F11 becomes
ROI013_E6. The ROI number (= CH.mcd acquisition ID) is kept.

The key label is taken from the stage position where the ROI was acquired
(data/mcd_acquisition_positions.csv, from scripts/mcd_positions.py), not from
the typed label, because two labels were typed wrong (ROI097 "F9" was acquired
at D9, ROI141 "J8" at J9). ROI117 re-acquired the core of ROI070 (same
position) and is marked as a repeat.

Only names and added label columns change; pixel data and measurements are not
recomputed. Default is a dry run that writes data/roi_key_map.csv and prints
every rename. Run with --apply to rename. Raw data and CH.mcd are not touched.

    /Users/rashid/Installations/venve/bin/python scripts/apply_tma_key.py
    /Users/rashid/Installations/venve/bin/python scripts/apply_tma_key.py --apply
"""

from __future__ import annotations

import argparse
import json
import re
import shutil
from datetime import datetime
from pathlib import Path

import pandas as pd

WORKFLOW = Path(__file__).resolve().parents[1]
ROWS = "ABCDEFGHIJK"
PATIENT_BY_COLUMN = "ONMLKJIHGFEDCBA"  # slide column 1..15 -> patient letter
SAMPLE_TYPE = {1: "High grade", 2: "High grade", 3: "High grade", 4: "Low grade",
               5: "High infiltration", 6: "High infiltration", 7: "Tumour periphery",
               8: "Tumour periphery", 9: "Tumour border", 10: "Tumour border", 11: "Non-neoplastic"}
EXPECTED_POSITION_FIXES = {97: "D9", 141: "J9"}
EXPECTED_REPEATS = {117: 70}
GRAPH_INVENTORY = dict(pd.read_csv(WORKFLOW / "results/10_graphml_inventory.csv", dtype=str)
                       [["source_file_name", "numbered_file_name"]].values)
KEY_COLUMNS = ["roi_number", "key_label", "patient_id", "core", "sample_type",
               "mcd_label", "acquired_position", "repeat_of"]


def key_from_position(position: str) -> tuple[str, int, str]:
    row, column = position[0], int(position[1:])
    core = ROWS.index(row) + 1
    patient = PATIENT_BY_COLUMN[column - 1]
    return patient, core, SAMPLE_TYPE[core]


def build_map() -> pd.DataFrame:
    images = pd.read_csv(WORKFLOW / "data/images.csv")
    positions = pd.read_csv(WORKFLOW / "data/mcd_acquisition_positions.csv")
    key = pd.read_csv(WORKFLOW / "data/reference/CH_TMA_key_from_TMA_Ref.csv")

    m = pd.DataFrame({"old_name": images["image"].str.removesuffix(".tiff")})
    if m["old_name"].str.contains(r"^ROI\d{3}_[A-K]\d{1,2}$").sum() != len(m):
        raise SystemExit("data/images.csv does not hold old-style names only; already renamed?")
    m["roi_number"] = m["old_name"].str[3:6].astype(int)
    m["mcd_label"] = m["old_name"].str.split("_").str[1]
    m = m.merge(positions[["acquisition_id", "mcd_label", "acquired_position", "distance_to_grid_point_um"]],
                left_on="roi_number", right_on="acquisition_id", how="left", suffixes=("", "_mcd"))
    assert m["acquisition_id"].notna().all(), "ROI without .mcd acquisition"
    assert (m["mcd_label"] == m["mcd_label_mcd"]).all(), "file label differs from .mcd label"

    fixes = {int(r.roi_number): r.acquired_position for r in m.itertuples() if r.mcd_label != r.acquired_position}
    assert fixes == EXPECTED_POSITION_FIXES, f"unexpected position corrections: {fixes}"

    keyed = m["acquired_position"].map(key_from_position)
    m["patient_id"] = [k[0] for k in keyed]
    m["core"] = [k[1] for k in keyed]
    m["sample_type"] = [k[2] for k in keyed]
    m["key_label"] = m["patient_id"] + m["core"].astype(str)
    m["new_name"] = [f"ROI{n:03d}_{lab}" for n, lab in zip(m["roi_number"], m["key_label"])]

    # independent check against the key table transcribed from the scan
    ref = key.set_index("printed_position")
    for r in m.itertuples():
        row = ref.loc[r.acquired_position]
        assert (row["label"], row["patient_id"], int(row["core"]), row["sample_type"]) == \
            (r.key_label, r.patient_id, r.core, r.sample_type), f"key mismatch for {r.old_name}"

    first = m.sort_values("roi_number").drop_duplicates("acquired_position").set_index("acquired_position")["roi_number"]
    m["repeat_of"] = [
        "" if first[r.acquired_position] == r.roi_number else f"ROI{first[r.acquired_position]:03d}"
        for r in m.itertuples()
    ]
    repeats = {int(r.roi_number): int(r.repeat_of[3:]) for r in m.itertuples() if r.repeat_of}
    assert repeats == EXPECTED_REPEATS, f"unexpected repeats: {repeats}"
    assert m["new_name"].is_unique and m["roi_number"].is_unique
    assert m.loc[m["repeat_of"] == "", "key_label"].is_unique, "two non-repeat ROIs share a key label"

    m["roi"] = m["new_name"]  # column read by scripts/02 --roi-name-map
    cols = ["roi", "old_name", "new_name"] + KEY_COLUMNS + ["distance_to_grid_point_um"]
    return m[cols].sort_values("roi_number").reset_index(drop=True)


def file_renames(m: pd.DataFrame) -> list[tuple[Path, Path]]:
    pairs = []
    for r in m.itertuples():
        if r.old_name == r.new_name:
            continue
        for folder, ext in [("data/img", ".tiff"), ("data/masks", ".tiff"), ("data/intensities", ".csv"),
                            ("data/regionprops", ".csv"), ("data/neighbors", ".csv"), ("data/graphs", ".graphml"),
                            ("qupath", ".ome.tiff"), ("qupath", "_cells.geojson")]:
            src = WORKFLOW / folder / f"{r.old_name}{ext}"
            if folder == "qupath" and not src.exists():
                continue
            pairs.append((src, WORKFLOW / folder / f"{r.new_name}{ext}"))
        numbered = GRAPH_INVENTORY[f"{r.old_name}.graphml"]  # e.g. 10_012_ROI013_F11.graphml (running index)
        pairs.append((WORKFLOW / "results/10_graphml" / numbered,
                      WORKFLOW / "results/10_graphml" / numbered.replace(r.old_name, r.new_name)))
    preview = WORKFLOW / "results/01b_channel_preview_ROI001_D13.png"
    if preview.exists():
        pairs.append((preview, WORKFLOW / "results/01b_channel_preview_ROI001_C4.png"))
    return pairs


def read_str(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, dtype=str, keep_default_na=False)


def add_key_columns(df: pd.DataFrame, roi_old: pd.Series, m: pd.DataFrame) -> pd.DataFrame:
    info = m.set_index("old_name")
    for col in KEY_COLUMNS:
        df[col] = roi_old.map(info[col]).astype(str).values
    return df


def rename_in_text(text: str, names: dict[str, str]) -> str:
    pattern = re.compile(r"ROI\d{3}_[A-K]\d{1,2}(?![0-9])")
    return pattern.sub(lambda mt: names.get(mt.group(0), mt.group(0)), text)


def apply(m: pd.DataFrame, pairs: list[tuple[Path, Path]]) -> None:
    stamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    backup = WORKFLOW / "archive" / f"pre_tma_key_{stamp}"
    backup.mkdir(parents=True)
    names = dict(zip(m["old_name"], m["new_name"]))

    # 1. back up every table and project file that is rewritten (small; images and masks are only renamed)
    tables = ["data/images.csv", "results/02_images.csv", "data/roi_marker_availability.csv",
              "results/01_roi_marker_availability.csv", "results/03_preprocessed_tiff_inventory.csv",
              "results/04_cellpose_mask_inventory.csv", "results/05_intensity_table_inventory.csv",
              "results/06_regionprops_table_inventory.csv", "results/07_neighbors_table_inventory.csv",
              "results/10_graphml_inventory.csv", "results/11_review_checks_per_roi.csv",
              "results/12_mcd_vs_tiff_per_channel.csv", "data/cells.csv", "results/08_cells.csv",
              "data/cells.h5ad", "results/09_cells.h5ad"]
    for t in tables:
        dst = backup / t
        dst.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(WORKFLOW / t, dst)
    shutil.copytree(WORKFLOW / "qupath_project", backup / "qupath_project")
    m.to_csv(backup / "roi_key_map.csv", index=False)

    # 2. rename files (content unchanged)
    for src, dst in pairs:
        if dst.exists():
            raise SystemExit(f"target exists, stopping: {dst}")
    log = []
    for src, dst in pairs:
        src.rename(dst)
        log.append(f"{src.relative_to(WORKFLOW)} -> {dst.relative_to(WORKFLOW)}")

    # 3. graphml: add the key as graph attributes (nodes and edges untouched)
    for folder in ["data/graphs", "results/10_graphml"]:
        for r in m.itertuples():
            name = f"{r.new_name}.graphml" if folder == "data/graphs" else \
                GRAPH_INVENTORY[f"{r.old_name}.graphml"].replace(r.old_name, r.new_name)
            path = WORKFLOW / folder / name
            text = path.read_text(encoding="utf-8")
            attrs = {"roi": r.new_name, "key_label": r.key_label, "patient_id": r.patient_id,
                     "core": str(r.core), "sample_type": r.sample_type, "mcd_label": r.mcd_label,
                     "acquired_position": r.acquired_position, "repeat_of": r.repeat_of}
            keys = "".join(f'<key id="g_{k}" for="graph" attr.name="{k}" attr.type="string"/>\n' for k in attrs)
            data = "".join(f'<data key="g_{k}">{v}</data>' for k, v in attrs.items())
            assert text.count("<graph ") == 1 and 'for="graph"' not in text
            text = text.replace("<graph ", keys + "<graph ", 1)
            text = re.sub(r"(<graph [^>]*>)", lambda mt: mt.group(1) + data, text, count=1)
            path.write_text(text, encoding="utf-8")

    # 4. tables
    for t in ["data/images.csv", "results/02_images.csv"]:
        df = read_str(WORKFLOW / t)
        df["image"] = df["image"].str.removesuffix(".tiff").map(names) + ".tiff"
        df.to_csv(WORKFLOW / t, index=False)
    for t in ["data/roi_marker_availability.csv", "results/01_roi_marker_availability.csv",
              "results/12_mcd_vs_tiff_per_channel.csv"]:
        df = read_str(WORKFLOW / t)
        df["roi"] = df["roi"].map(names)
        df.to_csv(WORKFLOW / t, index=False)
    for t in ["results/03_preprocessed_tiff_inventory.csv", "results/04_cellpose_mask_inventory.csv",
              "results/05_intensity_table_inventory.csv", "results/06_regionprops_table_inventory.csv",
              "results/07_neighbors_table_inventory.csv", "results/10_graphml_inventory.csv"]:
        df = read_str(WORKFLOW / t)
        for col in df.columns:
            if col.endswith("file_name") or col.endswith("relative_path"):
                df[col] = df[col].map(lambda s: rename_in_text(s, names))
        if t == "results/10_graphml_inventory.csv":  # graph attributes were added in step 3
            df["size_bytes"] = [str((WORKFLOW / p).stat().st_size) for p in df["numbered_relative_path"]]
        df.to_csv(WORKFLOW / t, index=False)
    df = read_str(WORKFLOW / "results/11_review_checks_per_roi.csv")
    old = df["roi"].copy()
    df["roi"] = old.map(names)
    add_key_columns(df, old, m).to_csv(WORKFLOW / "results/11_review_checks_per_roi.csv", index=False)
    for t in ["data/cells.csv", "results/08_cells.csv"]:
        df = read_str(WORKFLOW / t)
        old = df["Image"].str.removesuffix(".tiff")
        df["Image"] = old.map(names) + ".tiff"
        add_key_columns(df, old, m).to_csv(WORKFLOW / t, index=False)

    # 5. AnnData
    import anndata as ad
    for t in ["data/cells.h5ad", "results/09_cells.h5ad"]:
        a = ad.read_h5ad(WORKFLOW / t)
        old = a.obs["Image"].astype(str).str.removesuffix(".tiff")
        a.obs_names = [rename_in_text(o, names) for o in a.obs_names]
        a.obs["Image"] = (old.map(names) + ".tiff").values
        a.obs["image"] = a.obs["Image"].values
        info = m.set_index("old_name")
        for col in KEY_COLUMNS:
            a.obs[col] = old.map(info[col]).astype(str).values
        a.write_h5ad(WORKFLOW / t)

    # 6. QuPath project: point entries to the renamed files
    proj = WORKFLOW / "qupath_project"
    for path in [proj / "project.qpproj", *proj.glob("data/*/server.json"), *proj.glob("data/*/summary.json")]:
        text = path.read_text(encoding="utf-8")
        new = rename_in_text(text, names)
        if new != text:
            path.write_text(new, encoding="utf-8")

    m.to_csv(WORKFLOW / "data/roi_key_map.csv", index=False)
    (WORKFLOW / "logs").mkdir(exist_ok=True)
    (WORKFLOW / "logs" / f"apply_tma_key_{stamp}.log").write_text(
        f"backup: {backup.relative_to(WORKFLOW)}\n" + "\n".join(log) + "\n", encoding="utf-8")
    print(f"renamed {len(pairs)} files; tables updated; backup {backup.relative_to(WORKFLOW)}")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--apply", action="store_true", help="rename files and rewrite tables")
    args = parser.parse_args()

    m = build_map()
    pairs = file_renames(m)
    missing = [str(s.relative_to(WORKFLOW)) for s, _ in pairs if not s.exists()]
    if missing:
        raise SystemExit(f"{len(missing)} source files missing, e.g. {missing[:5]}")
    out = WORKFLOW / ("data/roi_key_map.csv" if args.apply else "logs/roi_key_map_DRYRUN.csv")
    out.parent.mkdir(exist_ok=True)
    m.to_csv(out, index=False)
    print(m[["old_name", "new_name", "patient_id", "core", "sample_type", "repeat_of"]].to_string(index=False))
    print(f"\n{len(m)} ROIs, {int((m.old_name == m.new_name).sum())} names unchanged, {len(pairs)} files to rename")
    print(f"map: {out.relative_to(WORKFLOW)}")
    if args.apply:
        apply(m, pairs)
    else:
        print("dry run: nothing renamed (use --apply)")


if __name__ == "__main__":
    main()
