# CGExplorer NEWS

## CGExplorer 0.4.0 (2026-11-22)

### New Features
* **Metric-Level Combining Scoring (`Metric_Combined`)**: Added `score_batch_metric_combined()`, `huber_pheno_stats()`, and `filter_low_growth_mutants()`, and extended `score_batch_dataset()` / `score_batch()` to combine plate-level phenotypes across replicates before deriving relative fitness and interaction scores:
  - **Average (Mean)**: Mean plate-level phenotype across replicates.
  - **Min (Worst-Case)**: Minimum plate-level phenotype across replicates.
  - **Max (Best-Case)**: Maximum plate-level phenotype across replicates.
  - **Robust Variance-Aware**: Huber IRLS robust location with Empirical Bayes moderated phenotypic variance (`d0 = 3`), delta-method SE propagation to `se_log2_cg`, screen-wide MAD scale, James-Stein shrinkage, and total-variance Z-scores.
* **Redesigned Scoring Module (`mod_scoring`)**:
  - Parallel mode selector for **One-vs-One** (choose control basis and replicate pair to compare) vs. **Merge (Metric-Level)** (choose `Average`, `Min`, `Max`, or `Variance-Aware`, with optional low-growth QC toggle off by default).
  - Explicit **Run Scoring** button that generates, plots (**Ranked Hits Plot**), tabulates (**Scoring Table**), and exports (**Download All Results (CSV)**) all unfiltered results.
  - Added **4-Strategy Comparison** tab showing Sensitive vs. Resistant hit counts and cross-strategy overlap across all four merge strategies.

---

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
