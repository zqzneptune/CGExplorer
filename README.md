# CGExplorer: Chemo-Genomic Screening Data Preprocessor & Visualizer

**CGExplorer** is an R package and Shiny application designed for lab researchers to quickly preprocess, perform quality control, analyze, and visualize high-throughput Chemo-Genomic Screening datasets (such as Tecan 384-well plate reader outputs for bacterial mutant libraries under chemical/drug treatments).

The package provides S4 data structures (`Plate`, `Batch`, `PlateRegistry`), automated data preprocessing pipelines, chemo-genomic interaction scoring, interactive `plotly` visualizations, and a full-featured Shiny GUI.

---

## Key Features

- **S4 Data Architecture**:
  - `Plate`: Manages raw Tecan reader data (growth, staining), 384-well layout metadata, computed kinetic metrics, and QC flags.
  - `Batch`: Groups plates across experimental replicates and stores multi-replicate interaction scoring results.
  - `PlateRegistry`: Manages persistent local state across user sessions and package upgrades.
- **Automated Data Assembly**: Pairs raw Tecan plate reader output files (`.txt`) with 384-well layout templates (`.xlsx`/`.tsv`).
- **Growth & Biofilm Metrics**: Calculates kinetic parameters (AUC, Max OD, Max Rate, Duration, Endpoint OD) and stained biofilm metrics (Mean OD, Relative Biofilm index).
- **Staining Detection**: Automated jump detection algorithm (`detect_staining_jump`) to identify staining timepoints.
- **Quality Control Pipeline**: Identifies contaminated blank wells, evaluates wild-type (WT) control growth, and flags anomalous wells.
- **Chemo-Genomic Interaction Scoring**: Scores chemical-gene interactions (linear & $\log_2$ interaction scores, $z$-scores, FDR adjustment) with single-replicate or batch-merged strategies (`mean`, `min`, `max`).
- **Interactive Visualizations**: Powered by `plotly` for interactive growth curves, rank-ordered hit plots, plate layout heatmaps, and mutant-level replicate overlays.
- **Data Export**: Direct UI export of hit tables, QC logs, and well-level metrics to CSV or Excel (`.xlsx`) files.

---

## Installation

```R
# Install pak if not already installed
if (!requireNamespace("pak", quietly = TRUE)) {
  install.packages("pak")
}

# Install CGExplorer package
pak::pak("zqzneptune/CGExplorer")
```

---

## Launching the Interactive Shiny App

```R
library(CGExplorer)

# Launch the Shiny application
run_cgexplorer_app()
```

---

## License

Distributed under the MIT License. Copyright (c) 2026 QINGZHOU ZHANG.
