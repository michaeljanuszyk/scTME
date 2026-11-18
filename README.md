# scTME

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

`scTME` is an R package that annotates the stromal and immune cells in a single-cell RNA-seq dataset using a pan-cancer atlas of the tumor microenvironment (TME). The reference atlas is built into the package, and annotations are available at two resolutions: 7 broad cell populations or 45 transcriptionally distinct subpopulations.

This repository also contains the code used to build and analyze the atlas in the accompanying article.

* **Article:** Januszyk, M., Lu, J., et al. Conservation of the tumour microenvironment across organs and species. *Nature* (2026). <https://doi.org/10.1038/s41586-026-11270-4>
* **Interactive atlas:** <https://tmeatlas.stanford.edu>
* **Archived version of record:** <https://doi.org/10.5281/zenodo.23027610>

This code was developed by the **Longaker Laboratory at Stanford University** ([Longaker Lab Website](https://www.longakerlab.com)). Please contact **Michael Januszyk** (januszyk@stanford.edu) with any questions or concerns.

---

## 📦 Installation

Install the Bioconductor dependencies, then install `scTME` from GitHub:

```r
if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c("SingleR", "SingleCellExperiment"))

# install.packages("devtools")
devtools::install_github("michaeljanuszyk/scTME")
```

---

## 🚀 Quick Start: Annotate a New Dataset

`scTME()` accepts a Seurat or SingleCellExperiment object and returns the predicted cell type for each cell (or for each cluster, if clusters are supplied).

```r
library(Seurat)
library(scTME)

# Example query: a 1,000-cell subset of the Nieto et al. (2021) tumor immune atlas, included with the package
query_path <- system.file("extdata", "nieto_2021_downsampled1k.rds", package = "scTME")
query_dataset <- readRDS(query_path)

# Broad populations (level = 1), one prediction per cluster
res <- scTME(obj = query_dataset, clusters = query_dataset$cluster, level = 1)
query_dataset$prediction <- unname(res)

# Subpopulations (level = 2, the default), one prediction per cell
res <- scTME(obj = query_dataset, level = 2)
```

A full walkthrough, including a comparison against the dataset's published labels, is in the package vignette (`vignettes/scTME_vignette_v2.Rmd`) and online at <https://tmeatlas.stanford.edu/Vignette>.

---

## 🧬 Reproducing the Atlas

The code used to generate and analyze the atlas is in the `code/` folder:

* `code/code_to_generate_atlas/`: dataset processing, quality control, integration, and annotation for the human and mouse atlases
* `code/code_to_analyze_atlas/`: downstream analyses (cross-species integration, CellChat, NicheNet, inferCNV, SCENIC, scMetabolism, and subpopulation enrichment)
* `resources/datasets_human.csv` and `resources/datasets_mouse.csv`: the datasets used in the article

To run this code in full, first download the publicly available GEO or ArrayExpress datasets and convert each sample from each study into a separate Seurat object. Save each object as the accession code + `.object.rds` (e.g., `GSM7114948.object.rds`), in folders named by study, inside a parent folder called `readyForSeurat`.

For example:

```text
readyForSeurat/
└── GSE190597/
    └── GSM7114948.object.rds
```

Instructions for running this code on a very small number of samples can be found in `code/code_to_generate_atlas/exampleSetup.R`, including expected output. Running the code in `tmeMetaAnalysis_human.R` on these samples takes approximately four hours.

### Hardware Requirements & Performance

We recommend using a large server to run the full pipeline. It took us roughly 6 weeks to run our code in its entirety using all datasets from `resources/datasets_human.csv` on Dual AMD EPYC 7713 processors (2.0 GHz, 64-cores/128-threads) with 2TB of NVMe SSD memory. We believe on average the runtime scales in an $O(n \log n)$ fashion, meaning that using only 50% of the samples would require considerably less than 50% of the runtime.

### Software & Environment Versions

The atlas pipeline relies on the following R and Python package versions. Annotating a new dataset with `scTME()` requires only the packages installed above.

#### R Packages
* `R` v4.2.2
* `Seurat` v5.0.1
* `Seurat` v4.9.9.9050 (for nichenetr and reticulate only)
* `SeuratWrappers` v0.4.0
* `BPCells` v0.1.0
* `SingleR` v2.0.0
* `scATOMIC` v2.0.2
* `zellkonverter` v1.8.0
* `reticulate` v1.35.0
* `nichenetr` v2.0.4
* `biomaRt` v2.54.1
* `CellChat` v1.5.0
* `infercnv` v1.14.2
* `scMetabolism` v0.2.1
* `SCENIC` v1.3.1
* `TCGAbiolinks` v3.18

#### Python Environments
* `Python` 3.9.19
* `scanpy` v1.10.1
* `scrublet` v0.2.3
* `scib` v1.1.5

---

## 📖 Citation

If you use `scTME` or the TME atlas in your work, please cite:

> Januszyk, M., Lu, J., et al. Conservation of the tumour microenvironment across organs and species. *Nature* (2026). https://doi.org/10.1038/s41586-026-11270-4
