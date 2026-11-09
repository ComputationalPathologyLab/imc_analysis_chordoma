# IMC analysis of the chordoma TMA (Cellpose)

Imaging mass cytometry (IMC) analysis of the chordoma (CH) tissue microarray, processed with the same workflow as the chondrosarcoma (CHS) project `261007_imcAnalysis_cellpose_Refined`, adapted to the CH panel.

## Data

- `data/raw/Tiffs/`: 153 ROI folders (ROI001–ROI154, ROI011 absent), copied from `Cananzi/20260116_TMA/IMC_data/CH/Tiffs/` on the imaging facility share. Each folder holds 43 single-channel `.ome.tiff` files (1000 × 1000 px, uint16). Read-only.
- `data/reference/CH_ROI_positions.png`: ROI positions on the TMA.

## Workflow

`notebooks/01_Complete IMC Data Analysis Workflow.ipynb` calls the numbered scripts in `scripts/`, writes logs to `logs/` and numbered outputs to `results/`. `ROI_SELECTION` in Step 0 sets the ROIs processed (`None` = all).

Changes relative to the CHS workflow:

- `--raw-dir` and `--roi` options in `01_create_panel_from_raw.py` and `02_preprocess_imc_images.py`.
- Nuclear composite DNA1, DNA2, Brachyury (no SOX9 in the CH panel); membrane composite S100, CD68, aSMA, CD31, CD3, CD66b. Provisional.
- TMA key applied (2026-11-09): images and all derived files are named `ROI<number>_<patient><core>`, e.g. raw folder `ROI013_F11` (slide position F11) -> `ROI013_E6` (patient E, core 6, high infiltration). Key source: `data/reference/Chordoma - TMA Ref.pdf`; per-ROI map with patient, core, sample type, typed .mcd label and acquired position: `data/roi_key_map.csv` (`scripts/apply_tma_key.py`, stage positions from `scripts/mcd_positions.py`). Patients are slide columns 1-15 read as O..A; cores are slide rows A-K read as 1-11 (1-3 high grade, 4 low grade, 5-6 high infiltration, 7-8 tumour periphery, 9-10 tumour border, 11 non-neoplastic).
- `percentile` composite scaling only.
- `acquisition_id` taken from the ROI number (review item 1.6).

Positions typed twice in the .mcd were resolved from stage coordinates: ROI097 (typed F9) was acquired at D9 = G4, ROI141 (typed J8) at J9 = G10. ROI117 is a second acquisition of the core of ROI070 (N7) and is marked `repeat_of = ROI070`. Not acquired: A11, B11, C11, D5, D11, E11, F3, F5, F11, G11, H11, I11, M11.

Steps 5–8 (Steinbock) require Docker Desktop to be running.
