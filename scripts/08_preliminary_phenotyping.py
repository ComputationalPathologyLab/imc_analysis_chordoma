#!/usr/bin/env python3
"""Preliminary phenotyping of the chordoma (CH) TMA by marker gating.

One threshold per marker is used for the whole array.

Steps
1. DNA status per core: a cell is low-DNA when arcsinh(DNA1) + arcsinh(DNA2) < 1;
   a core is DNA-negative when more than 50% of its cells are low-DNA.
   DNA-negative cores and the repeat acquisition ROI117 are left out.
   Low-DNA cells in the remaining cores are labelled "Low DNA" and not phenotyped.
2. A cell is positive for a marker when its mean count is above the threshold
   chosen from the arcsinh histograms (THRESHOLD).
3. Cells positive for 4 or more lineage markers are labelled "Mixed signal".
4. Gating rules are applied from top to bottom; the first match wins.
   CD56 is not used: most cells, including tumour cells, are CD56 positive.

    /Users/rashid/Installations/venve/bin/python scripts/08_preliminary_phenotyping.py
Output: results/phenotyping_preliminary/
"""

from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import tifffile

W = Path(__file__).resolve().parents[1]
OUT = W / "results" / "phenotyping_preliminary"
OUT.mkdir(parents=True, exist_ok=True)
plt.rcParams.update({"font.family": "Arial", "font.size": 10,
                     "axes.spines.top": False, "axes.spines.right": False})

THRESHOLD = {"Brachiury": 1, "S100": 3}          # every other marker: 2 counts
thr = lambda m: THRESHOLD.get(m, 2)
LINEAGE = ["Brachiury", "S100", "CD68", "CD163", "CD3", "CD4", "CD8", "CD31", "aSMA", "CD66b", "CD11c"]
ORDER = ["Tumour", "Macrophage", "Neutrophil", "T cell", "Dendritic cell", "Endothelium",
         "aSMA+ stroma", "Other cells", "Mixed signal", "Low DNA"]
COLOURS = {"Tumour": "#d62728", "Macrophage": "#2ca02c", "Neutrophil": "#ff7f0e", "T cell": "#1f77b4",
           "Dendritic cell": "#9467bd", "Endothelium": "#e377c2", "aSMA+ stroma": "#8c564b",
           "Other cells": "#c7c7c7", "Mixed signal": "#7f7f7f", "Low DNA": "#000000"}
CELL_MAP_ROIS = ["ROI003_K9", "ROI142_K1"]
SAMPLE_ORDER = ["High grade", "Low grade", "High infiltration", "Tumour periphery", "Tumour border", "Non-neoplastic"]

# 1. cells and DNA status
cells = pd.read_csv(W / "data/cells.csv", low_memory=False)
cells = cells[cells["repeat_of"].isna()].copy()
cells["Image"] = cells["Image"].str.replace(".tiff", "", regex=False)
low_dna = np.arcsinh(cells["DNA1"]) + np.arcsinh(cells["DNA2"]) < 1
frac = low_dna.groupby(cells["Image"]).mean()
status = pd.cut(frac, [-1, 0.2, 0.5, 1.1], labels=["valid", "intermediate", "dna_negative"]).astype(str)
cells["dna_status"] = cells["Image"].map(status)
cells["low_dna"] = low_dna
df = cells[cells["dna_status"] != "dna_negative"].copy()
print(f"{len(df)} cells in {df['Image'].nunique()} cores "
      f"({(status == 'dna_negative').sum()} DNA-negative cores left out)")

# 2. positive / negative
pos = lambda m: df[m].to_numpy() > thr(m)
neg = lambda m: ~pos(m)
myeloid = pos("CD68") | pos("CD163")

# 3-4. gating rules, first match wins
rules = {
    "Low DNA": df["low_dna"].to_numpy(),
    "Mixed signal": (df[LINEAGE] > 2).sum(axis=1).to_numpy() >= 4,
    "T cell": (df["CD3"].to_numpy() > 1) & (pos("CD4") | pos("CD8")) & neg("CD68") & neg("Brachiury"),
    "Neutrophil": pos("CD66b") & neg("CD68") & neg("Brachiury"),
    "Macrophage": myeloid & neg("Brachiury"),
    "Dendritic cell": pos("CD11c") & neg("Brachiury") & neg("S100"),
    "Endothelium": pos("CD31") & neg("aSMA") & neg("Brachiury") & neg("S100"),
    "Tumour": pos("Brachiury") | (pos("S100") & neg("aSMA")),
    "aSMA+ stroma": pos("aSMA") & neg("Brachiury"),
}
df["phenotype"] = np.select(list(rules.values()), list(rules.keys()), "Other cells")

keep_cols = ["Image", "Object", "key_label", "patient_id", "core", "sample_type", "dna_status",
             "centroid-0", "centroid-1", "area", "phenotype"]
df[keep_cols].to_csv(OUT / "preliminary_phenotypes.csv", index=False)
print(df["phenotype"].value_counts().to_string())

# per-core table (% of phenotyped cells, Low DNA and Mixed signal excluded from the denominator)
ph = df[~df["phenotype"].isin(["Low DNA", "Mixed signal"])]
core = pd.crosstab(ph["Image"], ph["phenotype"], normalize="index").mul(100).round(1)
core = core.reindex(columns=[o for o in ORDER if o in core.columns])
info = df.groupby("Image")[["key_label", "patient_id", "core", "sample_type", "dna_status"]].first()
core = info.join(df.groupby("Image").size().rename("cells")).join(core)
core.to_csv(OUT / "phenotypes_per_core.csv")

thr_table = pd.DataFrame({"threshold_counts": [thr(m) for m in LINEAGE + ["GranzymeB"]],
                          "percent_positive": [round((df[m] > thr(m)).mean() * 100, 1)
                                               for m in LINEAGE + ["GranzymeB"]]},
                         index=LINEAGE + ["GranzymeB"])
thr_table.to_csv(OUT / "thresholds.csv")

# figure: histograms with thresholds
fig, axs = plt.subplots(3, 4, figsize=(12, 7.5))
for ax, m in zip(axs.ravel(), LINEAGE + ["GranzymeB"]):
    ax.hist(np.arcsinh(df[m]), bins=150, color="grey")
    ax.set_yscale("log")
    ax.axvline(np.arcsinh(thr(m)), color="red", lw=1, ls="--")
    ax.set_title(f"{m.replace('Brachiury', 'Brachyury')} (> {thr(m)})")
    ax.set_xlabel("arcsinh(mean count)")
fig.tight_layout()
fig.savefig(OUT / "histograms_thresholds.png", dpi=200, bbox_inches="tight")
plt.close(fig)

# figure: composition by sample type
comp = pd.crosstab(ph["sample_type"], ph["phenotype"], normalize="index").mul(100)
comp = comp.reindex(index=[s for s in SAMPLE_ORDER if s in comp.index],
                    columns=[o for o in ORDER if o in comp.columns])
n = ph.groupby("sample_type").size()
fig, ax = plt.subplots(figsize=(8, 4.5))
left = np.zeros(len(comp))
for p in comp.columns:
    ax.barh(range(len(comp)), comp[p], left=left, color=COLOURS[p], label=p)
    left += comp[p].to_numpy()
ax.set_yticks(range(len(comp)))
ax.set_yticklabels([f"{s} (n={n[s]:,})" for s in comp.index])
ax.invert_yaxis()
ax.set_xlabel("% of phenotyped cells")
ax.set_xlim(0, 100)
ax.legend(bbox_to_anchor=(1.01, 1), loc="upper left", frameon=False)
fig.tight_layout()
fig.savefig(OUT / "composition_by_sample_type.png", dpi=200, bbox_inches="tight")
plt.close(fig)

# figure: mean marker expression per phenotype (scaled 0-1 per marker)
hm_markers = LINEAGE + ["GranzymeB", "CD56", "Ki67"]
hm = np.arcsinh(df.groupby("phenotype")[hm_markers].mean())
hm = hm.reindex([o for o in ORDER if o in hm.index])
hm = (hm - hm.min()) / (hm.max() - hm.min())
fig, ax = plt.subplots(figsize=(8, 4.5))
im = ax.imshow(hm.to_numpy(), cmap="Reds", aspect="auto", vmin=0, vmax=1)
ax.set_xticks(range(len(hm_markers)))
ax.set_xticklabels([m.replace("Brachiury", "Brachyury") for m in hm_markers], rotation=45, ha="right")
ax.set_yticks(range(len(hm)))
ax.set_yticklabels(hm.index)
fig.colorbar(im, ax=ax, label="mean arcsinh, scaled per marker")
fig.tight_layout()
fig.savefig(OUT / "marker_heatmap.png", dpi=200, bbox_inches="tight")
plt.close(fig)

# figure: cell maps (phenotype painted on the segmentation mask)
for roi in CELL_MAP_ROIS:
    mask = tifffile.imread(W / "data/masks" / f"{roi}.tiff")
    sub = df[df["Image"] == roi].set_index("Object")["phenotype"]
    lut = np.ones((int(mask.max()) + 1, 3))
    for obj, p in sub.items():
        lut[obj] = matplotlib.colors.to_rgb(COLOURS[p])
    rgb = lut[mask]
    fig, ax = plt.subplots(figsize=(7, 6))
    ax.imshow(rgb)
    ax.set_axis_off()
    present = [p for p in ORDER if p in set(sub)]
    ax.legend(handles=[matplotlib.patches.Patch(color=COLOURS[p], label=p) for p in present],
              bbox_to_anchor=(1.01, 1), loc="upper left", frameon=False)
    ax.set_title(roi.split("_")[1])
    fig.tight_layout()
    fig.savefig(OUT / f"cell_map_{roi}.png", dpi=200, bbox_inches="tight")
    plt.close(fig)

print(OUT)
