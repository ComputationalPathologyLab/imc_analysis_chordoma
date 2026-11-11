# IMC Analysis of the Chordoma TMA

![Workflow](https://img.shields.io/badge/workflow-Steinbock%20IMC-2f6f73)
![Status](https://img.shields.io/badge/status-segmentation%20complete-3b7a57)
![Data](https://img.shields.io/badge/data-chordoma%20TMA%2C%2015%20patients-5b6f95)
![Segmentation](https://img.shields.io/badge/segmentation-Cellpose%20cpsam__v2-7a5c99)
![Input](https://img.shields.io/badge/pixels-read%20from%20.mcd-8a6f3d)
![Phenotyping](https://img.shields.io/badge/phenotyping-preliminary%20gating-9a5c5c)

## Overview

This repository holds the imaging mass cytometry (IMC) analysis of the chordoma
(CH) tissue microarray (TMA) of the Humanitas sarcoma cohort. It takes the raw
Hyperion acquisitions to single-cell tables, applies the TMA key so that every
core is named by patient and core, and adds a preliminary rule-based
phenotyping.

The workflow is the one used for the chondrosarcoma TMA
([imc-analysis-cellpose-refined](https://github.com/ComputationalPathologyLab/imc-analysis-cellpose-refined)),
adapted to chordoma: Brachyury replaces SOX9 as the tumour nuclear marker, and
pixels are read directly from the `.mcd` file.

## Workflow

```text
CH.mcd (Hyperion acquisitions)
-> panel generation
-> hot-pixel filter, channel stacks (pixels read from .mcd)
-> TMA key: stage positions, ROI names by patient and core
-> Cellpose segmentation (nuclear + membrane composites)
-> intensity, morphology and neighbour measurement
-> CSV / AnnData / GraphML export with TMA key columns
-> DNA status per core
-> preliminary rule-based phenotyping
-> QuPath and napari views of the phenotypes
```

## Principles

- The `.mcd` file on the imaging facility share is opened read-only.
- Every step is a numbered script in `scripts/`, called from a notebook cell.
- Each script writes a log to `logs/` and numbered outputs to `results/`.
- Notebooks are kept unexecuted; they document the steps and their parameters.
- Raw data, archived runs and exploratory material are kept in `data/raw/`,
  which is not part of the repository.

## Notebooks

| Notebook | Role |
|---|---|
| `notebooks/01_Complete IMC Data Analysis Workflow.ipynb` | Panel, preprocessing, segmentation, measurement and export (scripts 01–07) |
| `notebooks/02_Preliminary phenotyping.ipynb` | DNA status per core, marker thresholds, gating rules, tables and figures, napari view |

## Dataset

The CH TMA holds 11 cores from each of 15 patients. Each core was punched from
an area chosen on H&E.

| Core | Sampled area |
|---|---|
| 1–3 | Tumour (high grade) |
| 4 | Tumour (low-grade position; chordoma is not graded) |
| 5–6 | High inflammatory infiltrate |
| 7–8 | Tumour periphery, tumour side |
| 9–10 | Tumour border, non-neoplastic side |
| 11 | Non-neoplastic tissue far from the tumour |

**TMA key.** Patients A–O are the slide columns read right to left (column 15
is patient A); cores 1–11 are the slide rows A–K read top to bottom. ROIs are
named `ROI<acquisition number>_<patient><core>`, for example `ROI013_E6`
(acquisition 13, patient E, core 6). The key is in
`data/reference/CH_TMA_key_from_TMA_Ref.csv`; the key per ROI is in
`data/roi_key_map.csv`.

**Acquisition.**

| Item | Value |
|---|---:|
| Acquisitions in `CH.mcd` | 153 |
| Repeat acquisition (ROI117 of core N7, first acquired as ROI070) | 1 |
| Cores acquired | 152 |
| Positions not acquired (A11, B11, C11, D5, D11, E11, F3, F5, F11, G11, H11, I11, M11) | 13 |
| Channels per ROI (35 antibodies, 2 DNA, 6 technical) | 43 |
| Pixel size | 1 µm |

**Label corrections.** The stage coordinates stored in the `.mcd`
(`scripts/mcd_positions.py`) were fitted to the TMA grid. Two ROIs had been
labelled with the wrong position at acquisition and were renamed by their
stage position: ROI097 (labelled F9, acquired at D9, key G4) and ROI141
(labelled J8, acquired at J9, key G10). All other ROIs sit at their labelled
position.

## Directory Structure

| Path | Purpose |
|---|---|
| `data/panel.csv` | Steinbock panel: channels, names, segmentation roles |
| `data/roi_marker_availability.csv` | Channels present in each ROI |
| `data/images.csv` | Image metadata (size, acquisition, stage coordinates) |
| `data/roi_key_map.csv` | Old name, new name, patient, core and sample type for every ROI |
| `data/mcd_acquisition_positions.csv` | Stage position and nearest grid point of every acquisition |
| `data/reference/` | TMA reference sheet and the key as a table |
| `data/img/` | Multi-channel TIFF stacks, one per ROI (files over 100 MB excluded) |
| `data/masks/` | Cellpose masks |
| `data/intensities/`, `data/regionprops/`, `data/neighbors/` | Per-cell measurements |
| `data/cells.h5ad` | Single-cell AnnData object with TMA key columns |
| `data/graphs/` | Spatial neighbour graphs (GraphML) |
| `notebooks/` | Documented workflow |
| `scripts/` | Numbered scripts |
| `logs/` | Run logs |
| `results/` | Numbered tables and figures |
| `docs/presentations/` | Slides and figures of the preliminary results |

`data/cells.csv`, `results/08_cells.csv` and the largest TIFF stacks exceed the
GitHub file limit and are regenerated by steps 2 and 5.

## Segmentation Panel

| Role | Channels | `cellpose` value |
|---|---|---:|
| Nuclear | DNA1 (191Ir), DNA2 (193Ir), Brachyury | 1 |
| Membrane / cytoplasm | S100, CD68, aSMA, CD31, CD3, CD66b | 2 |
| Technical, not analysed | 127I, 131Xe, 134Xe, 138Ba, 208Pb, 80ArAr | – |

Brachyury, a nuclear marker of chordoma cells, is added to the nuclear image so
that tumour cells are segmented in cores with a weak DNA signal.

## Workflow Modules

### 1. Panel Generation

Script: `scripts/01_create_panel_from_raw.py`

Outputs: `data/panel.csv`, `data/roi_marker_availability.csv`,
`results/01_panel.csv`, `results/01_roi_marker_availability.csv`

- Build the panel from the channel names of the raw export.
- Assign the nuclear and membrane segmentation markers.
- Flag technical channels with `analysis = 0`.

### 2. Preprocessing

Script: `scripts/02_preprocess_imc_images.py`

Outputs: `data/img/*.tiff`, `data/images.csv`, `results/02_images.csv`,
`results/03_preprocessed_tiff_inventory.csv`

- Read each acquisition from `CH.mcd` (float32, no rounding).
- Apply the steinbock hot-pixel filter (`hpf = 50`): a pixel more than 50 counts
  above all eight neighbours is set to the brightest neighbour.
- Stack the channels in panel order and name each ROI by the TMA key.

### 3. TMA Key

Scripts: `scripts/mcd_positions.py`, `scripts/apply_tma_key.py`

Outputs: `data/mcd_acquisition_positions.csv`, `data/roi_key_map.csv`

- Fit a tilted grid to the stage coordinates of all acquisitions and assign
  each acquisition to its nearest grid point.
- Check the expected corrections (ROI097, ROI141) and the repeat (ROI117).
- Rename every output file and add the key columns to the tables.

### 4. Cellpose Segmentation

Script: `scripts/03_segment_cellpose.py`

Outputs: `data/masks/*.tiff`, `results/04_cellpose_mask_inventory.csv`

- Build a nuclear and a membrane composite, each the pixel-wise maximum of its
  markers scaled between the 1st and 99.5th percentiles.
- Run Cellpose (`cpsam_v2`, flow threshold 0.4, cell probability threshold 0)
  on GPU, without retraining.

### 5. Single-Cell Measurement and Export

Scripts: `scripts/04_measure_intensities.py`, `scripts/05_measure_regionprops.py`,
`scripts/06_measure_neighbors.py`, `scripts/07_export_data.py`

Outputs: `data/intensities/`, `data/regionprops/`, `data/neighbors/`,
`data/cells.h5ad`, `data/graphs/`, `results/05`–`results/09`

- Mean intensity of every channel over the pixels of each cell.
- Area, centroid, axis lengths and eccentricity.
- Neighbours: cells whose outlines, expanded by 4 px, touch.
- Export as CSV, AnnData and GraphML, with patient, core and sample type.

### 6. DNA Status and Preliminary Phenotyping

Script: `scripts/08_preliminary_phenotyping.py` (same steps as notebook 02)

Outputs: `results/phenotyping_preliminary/`

- A cell is low-DNA when arcsinh(DNA1) + arcsinh(DNA2) < 1. A core is valid
  below 20% low-DNA cells, intermediate at 20–50% and DNA-negative above 50%.
  DNA-negative cores and the repeat ROI117 are left out.
- One positivity threshold per marker, read from the arcsinh histograms:
  Brachyury 1 count, S100 3 counts, all others 2 counts. CD56 is not used,
  because most cells, including tumour cells, are CD56 positive.
- Cells positive for 4 or more lineage markers are labelled *Mixed signal*.
- Gating rules, first match wins:

| Cell type | Rule |
|---|---|
| T cell | CD3 > 1 and (CD4 or CD8), CD68 negative, Brachyury negative |
| Neutrophil | CD66b, CD68 negative, Brachyury negative |
| Macrophage | CD68 or CD163, Brachyury negative |
| Dendritic cell | CD11c, Brachyury and S100 negative |
| Endothelium | CD31, aSMA, Brachyury and S100 negative |
| Tumour | Brachyury, or S100 with aSMA negative |
| aSMA+ stroma | aSMA, Brachyury negative |
| Other cells | no rule matched |

### 7. Viewing the Phenotypes

Script: `scripts/09_phenotypes_to_qupath.py`

Outputs: `results/phenotyping_preliminary/qupath/`

- One OME-TIFF with 13 named channels and one GeoJSON of classified cell
  outlines per core, for QuPath.
- Notebook 02, step 7, opens a core in napari with the phenotypes over the
  marker images.

## Numbered Results

| Output | Description |
|---|---|
| `results/01_panel.csv` | Panel copy |
| `results/01_roi_marker_availability.csv` | Channels per ROI |
| `results/01b_channel_preview_ROI001_C4.png` | Channel preview of one core |
| `results/02_images.csv` | Image metadata copy |
| `results/03_preprocessed_tiff_inventory.csv` | TIFF stack inventory |
| `results/04_cellpose_mask_inventory.csv` | Mask inventory |
| `results/05_intensity_table_inventory.csv` | Intensity table inventory |
| `results/06_regionprops_table_inventory.csv` | Region property inventory |
| `results/07_neighbors_table_inventory.csv` | Neighbour table inventory |
| `results/09_cells.h5ad` | Exported AnnData object |
| `results/phenotyping_preliminary/preliminary_phenotypes.csv` | One row per cell with its phenotype |
| `results/phenotyping_preliminary/phenotypes_per_core.csv` | Cell-type percentages per core |
| `results/phenotyping_preliminary/thresholds.csv` | Thresholds and percentage of positive cells |
| `results/phenotyping_preliminary/*.png` | Histograms, composition by sample type, marker heatmap, cell maps |

## Current Status

| Stage | Status |
|---|---|
| Panel generation | Complete |
| Preprocessing from `.mcd` | Complete |
| TMA key and label corrections | Complete |
| Cellpose segmentation | Complete |
| Single-cell measurement and export | Complete |
| DNA status per core | Complete |
| Preliminary phenotyping | Complete (first pass) |
| H&E review and artifact annotation by the pathologist | In progress |
| Phenotyping review on the images | Planned |
| Spatial analysis | Planned |

## Current Summary

| Metric | Value |
|---|---:|
| Cores with cells (ROI117 excluded) | 151 |
| Segmented cells (ROI117 excluded) | 169,708 |
| Median cells per core | 1,046 |
| Cores with valid / intermediate / negative DNA signal | 122 / 5 / 24 |
| Cores used for phenotyping | 127 |
| Cells used for phenotyping | 166,082 |
| Tumour cells, tumour cores 1–4 → border → non-neoplastic | 63% → 13% → 3% |

18 of the 24 DNA-negative cores come from two patients (all cores of patient E,
8 of patient B). In core E6, for example, Brachyury-positive tumour cells are
visible although the DNA signal is almost absent.

## Interpretation Boundaries

- No artifact annotation has been applied yet; the phenotypes include cells
  touched by metal deposits or antibody aggregates.
- One threshold per marker is used for the whole array.
- The phenotypes have not yet been reviewed on the images by a pathologist.
- Core 4 is a tumour core; chordoma is not graded, so it is pooled with
  cores 1–3 in comparisons between areas.
- Cell-type percentages are exploratory and are not yet compared with the IHC
  counts of the same cores.

## Suggested Next Step

Apply the pathologist's H&E review and artifact annotation, as done for
chondrosarcoma, then repeat the phenotyping on clean cells and review the cell
types on the images before any spatial analysis.

## References

- Windhager J, Zanotelli VRT, Schulz D, et al. An end-to-end workflow for
  multiplexed image processing and analysis. *Nat Protoc*. 2023;18(11):3565–3613.
  https://doi.org/10.1038/s41596-023-00881-0
- Stringer C, Wang T, Michaelos M, Pachitariu M. Cellpose: a generalist algorithm
  for cellular segmentation. *Nat Methods*. 2021;18(1):100–106.
  https://doi.org/10.1038/s41592-020-01018-x
