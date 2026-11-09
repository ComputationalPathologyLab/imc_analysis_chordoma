#!/usr/bin/env python3
"""Stage position of every acquisition in CH.mcd, and the TMA grid position it sits on.

The ROI label typed at acquisition (e.g. "F9") can be wrong. The stage
coordinates stored in the .mcd show where each ROI was really acquired. A
tilted grid (x, y = a + b*column + c*row) is fitted to all acquisitions, with
outliers (> 400 um) left out of the fit, and every acquisition is assigned to
the nearest grid point. The .mcd is opened read-only.

    /Users/rashid/Installations/venve/bin/python scripts/mcd_positions.py
Output: data/mcd_acquisition_positions.csv
"""

from pathlib import Path

import numpy as np
import pandas as pd
from readimc import MCDFile

WORKFLOW = Path(__file__).resolve().parents[1]
MCD = Path("/Volumes/ricerca/Unita_Imaging/Shared/Hyperion/Projects/Cananzi/20260116_TMA/IMC_data/CH/CH.mcd")
ROWS = "ABCDEFGHIJK"

rows = []
with MCDFile(MCD) as f:
    for slide in f.slides:
        for a in slide.acquisitions:
            m = a.metadata
            rows.append({"acquisition_id": a.id, "mcd_label": a.description,
                         "start_x": float(m["ROIStartXPosUm"]) / 1000, "start_y": float(m["ROIStartYPosUm"]) / 1000,
                         "end_x": float(m["ROIEndXPosUm"]), "end_y": float(m["ROIEndYPosUm"])})
d = pd.DataFrame(rows)
d = d[d["mcd_label"].str.fullmatch(r"[A-K]\d{1,2}")].copy()  # drops the "placeholder" acquisition
d["centre_x_um"] = (d["start_x"] + d["end_x"]) / 2
d["centre_y_um"] = (d["start_y"] + d["end_y"]) / 2
r = d["mcd_label"].str[0].map(ROWS.index).to_numpy()
c = d["mcd_label"].str[1:].astype(int).to_numpy()

keep = np.ones(len(d), bool)
X = np.c_[np.ones(len(d)), c, r]
for _ in range(3):
    bx = np.linalg.lstsq(X[keep], d["centre_x_um"].to_numpy()[keep], rcond=None)[0]
    by = np.linalg.lstsq(X[keep], d["centre_y_um"].to_numpy()[keep], rcond=None)[0]
    keep = np.hypot(X @ bx - d["centre_x_um"], X @ by - d["centre_y_um"]).to_numpy() < 400

grid = [(ri, ci) for ri in range(11) for ci in range(1, 16)]
gx = np.array([bx @ [1, ci, ri] for ri, ci in grid])
gy = np.array([by @ [1, ci, ri] for ri, ci in grid])
pos, dist = [], []
for x, y in zip(d["centre_x_um"], d["centre_y_um"]):
    dd = np.hypot(gx - x, gy - y)
    i = int(np.argmin(dd))
    pos.append(f"{ROWS[grid[i][0]]}{grid[i][1]}")
    dist.append(round(float(dd[i]), 1))
d["acquired_position"] = pos
d["distance_to_grid_point_um"] = dist
out = WORKFLOW / "data/mcd_acquisition_positions.csv"
d[["acquisition_id", "mcd_label", "acquired_position", "distance_to_grid_point_um",
   "centre_x_um", "centre_y_um"]].round(1).to_csv(out, index=False)
print(f"{len(d)} acquisitions; label differs from stage position: "
      f"{d.loc[d['mcd_label'] != d['acquired_position'], ['acquisition_id', 'mcd_label', 'acquired_position']].values.tolist()}")
print(out)
