# Cell Overlap & ROI Annotation Tools (MATLAB)

This repository contains two App Designer applications for microscopy image inspection and analysis. Both are designed for **reproducible, config‑driven** workflows in cellular imaging pipelines.

* **Cell_overlap_V1_1** — two‑channel cell segmentation, manual mask refinement, ROI‑level overlap quantification, **batch processing**, and export to PNG/XLSX.
* **Image_Segament_Annotation_V1** — multi‑channel viewer with LUT control, manual ROI drawing/annotation (freehand/polygon), area computation (µm²), and export of merged overlays, per‑ROI masked TIFFs, and an ROI summary.

---

## What’s New (Batch Processing)

**Cell_overlap_V1_1** now supports a guided batch workflow:

* **Batch Import** multiple `.tif/.tiff` files; optional **Randomize** to blind file order.
* **Save & Next** exports per‑image outputs (3‑panel PNG + per‑image XLSX) and advances.
* **Skip** bypasses current image; **Quit** stops the session without losing collected results.
* **Save all** writes an aggregate summary Excel file: `Summary_of_all_tif_overlaping.xlsx` with columns:
  `ImagePrefix, Ch1, Ch2, Ch1_CellCount, Ch2_CellCount, ROI_Overlap_Count, ROI_Overlap_Percent, Overlap_Threshold_Pct, Dice, IoU`.
* **Config persistence:** calibration, segmentation, morphology, and overlap parameters **persist across images** during the batch and only change if you adjust controls or **Load** a new JSON.

---

## System Requirements

* **MATLAB** R2021b or later (App Designer)
* **Image Processing Toolbox**
* Excel I/O via MATLAB’s built‑in `writetable` for `.xlsx`
* OS: Windows, macOS, or Linux
* Optional utilities (if present on the MATLAB path, otherwise built‑in fallbacks are used by *Cell_overlap_V1_1*):

  * `op_threshold`, `op_sizeFilter`, `op_morph`

---

## Installation

1. **Clone or download** this repository.
2. In MATLAB, **add the project folder to the path**:

   ```matlab
   addpath(genpath('path/to/repo'))
   ```
3. Launch either app:

   ```matlab
   app = Cell_overlap_V1_1;            % segmentation + overlap + batch
   app = Image_Segament_Annotation_V1; % viewer + ROI annotation
   ```

---

## Data Requirements

* Input format: `.tif` / `.tiff`

  * Single image with multiple channels (SamplesPerPixel > 1) **or** multi‑page TIFF (one page per channel)
* Pixel intensities are normalized to double `[0,1]` for display; exports preserve original bit depth where applicable

---

## Quick Start

* Use **Image_Segament_Annotation_V1** for visual ROI definition and export of ROI‑masked TIFFs.
* Use **Cell_overlap_V1_1** for segmentation‑based cell detection, overlap quantification, config‑based reproducible analyses, and **batch processing**.

---

## App 1 — Cell_overlap_V1_1

### Overview

Two‑channel cell segmentation and overlap quantification with interactive mask editing and **batch mode**.

### Key Features

* Multi‑channel TIFF import; per‑channel processing
* µm/pixel **calibration** (drives area threshold units)
* Per‑channel segmentation: **auto/local or manual threshold**, size filter, morphology (erode/dilate), optional **watershed** splitting
* Manual editing tools: **Brush**, **Eraser**, **Delete‑lasso** (single‑interaction, click‑to‑start / click‑to‑finish)
* Overlap rule: **“Ch2 ROI inside Ch1 ROI ≥ Overlap%”**
* Pixel‑level metrics: **Dice** and **IoU**
* **Batch mode** (see *What’s New* above): Batch Import, Randomize, Save & Next, Skip, Quit, Save all
* **Config persistence across images** in a batch; **Save/Load** config JSON for reproducibility

### Typical Workflow (Single Image)

1. **Import** TIFF
2. (Optional) Enter **µm/pixel** calibration
3. Select **Ch1** and **Ch2** for comparison
4. Adjust segmentation (threshold/area/morphology/watershed)
5. Refine with manual tools
6. Set **Overlap%**
7. **Export** PNG + XLSX
8. (Optional) **Save** config JSON for later reuse / batch consistency

### Typical Workflow (Batch)

1. Click **Batch Import** and select multiple TIFFs (optionally **Randomize**)
2. For each image:

   * Adjust settings (if desired) and/or use manual tools
   * Click **Save & Next** to export current outputs and move on, or **Skip** to bypass
3. At any time, **Quit** to stop; already exported items remain
4. When finished, click **Save all** to write the aggregate summary `Summary_of_all_tif_overlaping.xlsx`

> **Note:** During batch, the current settings (including any slider/checkbox tweaks) are carried forward to the next image, unless you explicitly **Load** a different config.

### Output Files

* `*_ROIs_post_process.png` — 1×3 panel: Ch1, Ch2, Overlay & Results
* `*_Summary_CellCounting.xlsx` — per‑image counts: Ch1/Ch2 cell counts, ROI overlap counts/%, Overlap_Threshold_Pct, Dice, IoU
* `Summary_of_all_tif_overlaping.xlsx` — aggregate summary across the batch (one row per saved image)

---

## App 2 — Image_Segament_Annotation_V1

### Overview

Multi‑channel viewer and annotation tool for drawing/naming ROIs and exporting ROI‑masked data.

### Key Features

* Per‑channel LUT and intensity adjustment; merged composite view (up to 4 visible channels)
* ROI tools: **freehand**, **polygon**, **edit**, **undo**, **delete**
* ROI table: ID, name, **area (µm²)**, visibility
* Export merged overlay PNG, per‑ROI masked TIFFs, and ROI summary XLSX

### Typical Workflow

1. **Import** TIFF
2. (Optional) Set **µm/pixel** calibration
3. Adjust channel LUTs and intensities; choose the active view
4. Draw ROI (closed automatically), enter name, **Add**
5. Edit/Undo/Delete or manage via the table
6. **Export** outputs

### Output Files

* `*_Merged_ROIs.png` — composite with ROI overlays
* `masked/` — per‑ROI multi‑page TIFFs (channels preserved; outside ROI set to black)
* `*_ROIs_Summary.xlsx` — ROI list, display parameters, image metadata

---

## Troubleshooting

* **Export disabled:** Import an image first
* **Incorrect ROI area:** Check **µm/pixel** calibration
* **Edits not applying:** Ensure the correct axes/view is active
* **Channel missing (Annotation app):** Verify channel visibility controls
* **Save all empty:** In batch, use **Save & Next** at least once to collect rows before **Save all**

---

## Reproducibility

* Use **Save/Load config** in *Cell_overlap_V1_1* to capture **all** analysis settings
* During batch, the loaded/adjusted settings **persist** across images to ensure consistent processing
* Record MATLAB and toolbox versions alongside PNG/XLSX outputs for traceability

---

## License

This project is licensed under the **MIT License**.

---

## Citation

If used in publications, cite **MATLAB** and the **Image Processing Toolbox**. Acknowledge the use of these custom tools in Methods sections as appropriate.

---

## Contact

**Ping Dong**
Email: [superdongping@gmail.com](mailto:superdongping@gmail.com), [pingdong@unc.edu](mailto:pingdong@unc.edu)
