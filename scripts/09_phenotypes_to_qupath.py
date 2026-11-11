#!/usr/bin/env python3
"""Export cores with their phenotypes for viewing in QuPath.

For each core two files are written to results/phenotyping_preliminary/qupath/:
    <roi>.ome.tiff             selected marker channels, with marker names, 1 px = 1 um
    <roi>_phenotypes.geojson   one polygon per cell, classified by its phenotype
                               (from results/phenotyping_preliminary/preliminary_phenotypes.csv)

    /Users/rashid/Installations/venve/bin/python scripts/09_phenotypes_to_qupath.py
    /Users/rashid/Installations/venve/bin/python scripts/09_phenotypes_to_qupath.py --roi ROI003_K9 ROI142_K1

In QuPath: open the .ome.tiff, then drag the matching _phenotypes.geojson onto the image.
"""

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
import tifffile
from skimage.measure import find_contours, regionprops

W = Path(__file__).resolve().parents[1]
OUT = W / "results" / "phenotyping_preliminary" / "qupath"
CHANNELS = ["DNA1", "DNA2", "Brachiury", "S100", "CD68", "CD163", "CD66b", "CD3", "CD4", "CD8",
            "aSMA", "CD31", "CD11c"]
COLOURS = {"Tumour": "#d62728", "Macrophage": "#2ca02c", "Neutrophil": "#ff7f0e", "T cell": "#1f77b4",
           "Dendritic cell": "#9467bd", "Endothelium": "#e377c2", "aSMA+ stroma": "#8c564b",
           "Other cells": "#c7c7c7", "Mixed signal": "#7f7f7f", "Low DNA": "#000000"}
rgb = lambda h: [int(h[i:i + 2], 16) for i in (1, 3, 5)]


def cell_polygons(mask):
    """Outline of every cell in the mask as a closed ring of (x, y) points."""
    for region in regionprops(mask):
        r0, c0, _, _ = region.bbox
        contours = find_contours(np.pad(region.image, 1).astype(float), 0.5)
        if not contours:
            continue
        ring = [[float(c + c0 - 0.5), float(r + r0 - 0.5)] for r, c in max(contours, key=len)]
        if ring[0] != ring[-1]:
            ring.append(ring[0])
        if len(ring) >= 4:
            yield region.label, ring


parser = argparse.ArgumentParser()
parser.add_argument("--roi", nargs="*", default=["ROI003_K9", "ROI142_K1"])
args = parser.parse_args()

panel = pd.read_csv(W / "data/panel.csv").query("keep == 1").reset_index(drop=True)
index = [panel.index[panel["name"] == c][0] for c in CHANNELS]
names = [c.replace("Brachiury", "Brachyury") for c in CHANNELS]
pheno = pd.read_csv(W / "results/phenotyping_preliminary/preliminary_phenotypes.csv")
OUT.mkdir(parents=True, exist_ok=True)

for roi in args.roi:
    stack = tifffile.imread(W / "data/img" / f"{roi}.tiff")[index].astype(np.float32)
    tifffile.imwrite(OUT / f"{roi}.ome.tiff", stack, ome=True, compression="zlib",
                     metadata={"axes": "CYX", "Channel": {"Name": names},
                               "PhysicalSizeX": 1.0, "PhysicalSizeXUnit": "µm",
                               "PhysicalSizeY": 1.0, "PhysicalSizeYUnit": "µm"})

    mask = tifffile.imread(W / "data/masks" / f"{roi}.tiff")
    ph = pheno[pheno["Image"] == roi].set_index("Object")["phenotype"]
    features = []
    for label, ring in cell_polygons(mask):
        if label not in ph.index:
            continue
        p = ph[label]
        features.append({"type": "Feature",
                         "geometry": {"type": "Polygon", "coordinates": [ring]},
                         "properties": {"objectType": "detection", "name": f"cell {label}",
                                        "classification": {"name": p, "color": rgb(COLOURS[p])}}})
    (OUT / f"{roi}_phenotypes.geojson").write_text(json.dumps({"type": "FeatureCollection", "features": features}))
    print(f"{roi}: {len(features)} cells, {ph.value_counts().to_dict()}")
print(OUT)
