# Cell Overlap & ROI Annotation Tools (MATLAB)

This repository contains two App Designer applications for microscopy image inspection and analysis:

1. **Cell\_overlap\_V1** — two-channel cell segmentation, manual mask refinement, and ROI-level overlap quantification with export to PNG/XLSX.
2. **Image\_Segament\_Annotation\_V1** — multi-channel viewer with LUT control, manual ROI drawing/annotation (freehand/polygon), area computation (µm²), and export of merged overlays, masked per-ROI TIFFs, and an ROI summary.

Both apps are designed for reproducible workflows in cellular imaging pipelines.

---

## System Requirements

* **MATLAB** R2021b or later (App Designer).
* **Image Processing Toolbox**.
* **Excel I/O** via MATLAB’s built-in `writetable` for `.xlsx`.
* OS: Windows, macOS, or Linux.

> Optional: if utility functions `op_threshold`, `op_sizeFilter`, `op_morph` exist on the MATLAB path, **Cell\_overlap\_V1** will use them; otherwise, built-in fallbacks are used.

---

## Installation

1. Clone or download this repository.
2. Open MATLAB and **add the project folder to the path** (e.g., `addpath(genpath('path/to/repo'))`).
3. Launch either app by running in the Command Window:

   ```matlab
   app = Cell_overlap_V1;
   app = Image_Segament_Annotation_V1;
   ```

---

## Data Requirements

* **Input format**: `.tif` / `.tiff`

  * Single image with 3D planes (SamplesPerPixel > 1) or multi-page TIFF (one page per channel).
* Pixel intensities are normalized to **double \[0,1]** for display; export preserves **8/16-bit** depth where applicable.

---

## Quick Start

* Use **Image\_Segament\_Annotation\_V1** for visual ROI definition and export of ROI-masked TIFFs.
* Use **Cell\_overlap\_V1** for segmentation-based cell detection, overlap quantification, and config-based reproducible analyses.

---

## App 1 — Cell\_overlap\_V1

### Overview

Two-channel cell segmentation and overlap quantification with interactive mask editing.

### Key Features

* Multi-channel TIFF import.
* µm/pixel calibration for area thresholds.
* Per-channel segmentation: thresholding, area, morphology, watershed.
* Manual editing: brush, eraser, delete (lasso).
* Overlap rule: Ch2 ROI inside Ch1 ROI ≥ threshold %.
* Reports Dice, IoU, and overlap counts.
* Save/Load config JSON for reproducibility.

### Workflow

1. Import TIFF.
2. (Optional) Enter µm/pixel.
3. Select **Ch1** and **Ch2** for comparison.
4. Adjust segmentation settings.
5. Refine with manual tools.
6. Set **Overlap%**.
7. Export PNG + XLSX.
8. Save/Load config for batch consistency.

### Output Files

* `*_ROIs_post_process.png` — 1×3 panel: Ch1, Ch2, Overlay & Results.
* `*_Summary_CellCounting.xlsx` — cell counts, overlap counts, overlap %, Dice, IoU.

---

## App 2 — Image\_Segament\_Annotation\_V1

### Overview

Multi-channel viewer and annotation tool for drawing/naming ROIs and exporting ROI-masked data.

### Key Features

* Per-channel LUT and intensity adjustment.
* Merged composite view with up to 4 visible channels.
* ROI tools: freehand, polygon, edit, undo, delete.
* ROI table: ID, name, area (µm²), visibility.
* Export merged overlay PNG, per-ROI masked TIFFs, and ROI summary XLSX.

### Workflow

1. Import TIFF.
2. (Optional) Set µm/pixel.
3. Adjust channel LUTs and intensities.
4. Select view (click to activate).
5. Draw ROI (closed automatically), enter name, Add.
6. Edit/Undo/Delete or manage via table.
7. Export outputs.

### Output Files

* `*_Merged_ROIs.png` — composite with ROI overlays.
* `masked/` — per-ROI multi-page TIFFs, channels preserved, outside ROI set to black.
* `*_ROIs_Summary.xlsx` — ROI list, display parameters, image metadata.

---

## Troubleshooting

* **Export disabled**: Import an image first.
* **Incorrect ROI area**: Check µm/pixel calibration.
* **Edits not applying**: Ensure you clicked the correct active axes.
* **Channel missing**: Verify Ch1–Ch4 checkboxes in Annotation app.

---

## Reproducibility

* Use **Save/Load config** in **Cell\_overlap\_V1** to capture all analysis settings.
* Record MATLAB and toolbox versions alongside PNG/XLSX outputs for traceability.

---

## License

This project is licensed under the [MIT License](./LICENSE).

---

## Citation

If used in publications, cite MATLAB and the Image Processing Toolbox. Acknowledge the use of these custom tools in Methods sections where appropriate.

---

## Contact

For questions or feedback, contact **Ping Dong** at:

* [superdongping@gmail.com](mailto:superdongping@gmail.com)
* [pingdong@unc.edu](mailto:pingdong@unc.edu)
