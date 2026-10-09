# CGExplorer NEWS

## CGExplorer 0.3.2 (2026-10-09)

### New Features
* **Multi-Page Growth Curve PDF Export**: Added `export_batch_growth_curves_pdf()` to export comprehensive, publication-ready growth curve PDFs for all strains across all plates in a batch:
  - Generates one page per strain with facet panels for each plate and individual trajectories for replicate wells.
  - Controls (`WT`, `Blank`, etc.) are rendered first, followed by mutants alphabetically.
  - Uses dynamic per-page y-axis scaling with small margins (~8%) fitted to the strain's own data range, preventing blank upper space for slow-growing or low-yield strains.
* **Batch Overview Export UI**: Added a dedicated export card in `mod_batch_overview` allowing interactive downloading of Raw OD or Corrected OD multi-page growth curve reports with a live Shiny progress indicator.

### Improvements & Bug Fixes
* **Data Fidelity (No Normalization)**:
  - Raw OD mode outputs unmodified absorbance readings directly from instrument assays.
  - Corrected OD mode performs direct blank subtraction ($OD_{raw} - OD_{blank\_median}$) without clipping, ceiling, or artificial floors (e.g. `pmax(..., 0.01)`). Negative blank-corrected values are preserved accurately.
  - Retains all control wells (`Blank`, `Empty`) across both modes.

---

## CGExplorer 0.3.1 (2026-10-08)

### Quality Control & Plate Validation
* **Automated Failure Detection**: `compute_metrics()` now flags plate QC failure if all blank wells are contaminated or if all wild-type (WT) control wells fail to grow.
* **Manual QC Override**: Users can inspect and toggle a plate's status (`PASSED` / `FAILED`) directly within the Plate Card interface.
* **Batch Quarantine**: Added `validate_batch_plates()` to forbid including `FAILED` plates in downstream batches.
* **Layout Consistency Enforcement**: `validate_batch_plates()` verifies that all plates within a single batch share an identical 384-well layout.

### Scoring & UI Enhancements
* **Intra-Plate Scoring Gating**: Dynamic UI validation prevents selecting the intra-plate (`NoDrug_Wells`) strategy when layouts lack `_NoDrug` control wells.
* **Plate Card Recalculation**: Plate edits automatically trigger metric and QC recomputation.
* **Navigation**: Restored and activated **Batch Overview** and **Scoring** dashboard menu items.

---

## CGExplorer 0.3.0 (2026-08-09)

### Initial Architecture & Core Features
* **S4 Data Architecture**: Core classes `Plate`, `Batch`, and persistent container `PlateRegistry`.
* **Tecan Reader & Layout Parsing**: Built tidy parsers for Tecan OD raw matrix outputs and 96/384-well plate layouts.
* **Staining Jump Detection**: Automatic staining detection (`detect_staining_jump()`) for crystal violet biofilm assays.
* **Interactive Shiny Dashboard**: Golem-structured modules for plate browsing, batch building, quality control visualization, and chemo-genomic interaction scoring.
