# scTME

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

This code was developed by the **Longaker Laboratory at Stanford University** ([Longaker Lab Website](https://www.longakerlab.com)).

Please contact **Michael Januszyk** (januszyk@stanford.edu) directly with any questions or concerns.

---

## 📋 Software & Environment Versions

The pipeline relies on the following R and Python package versions:

### R Packages
* `R` v4.2.2
* `Seurat` v5.0.1
* `Seurat` v4.9.9.9050 (for nichenetr and reticulate only)
* `SeuratWrampers` v0.4.0
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

### Python Environments
* `Python` 3.9.19
* `scanpy` v1.10.1
* `scrublet` v0.2.3
* `scib` v1.1.5

---

## 💻 Hardware Requirements & Performance

We recommend using a large size server to run this code. It took us roughly 6 weeks to run our code in its entirety using all datasets from resources/datasets_human.csv on Dual AMD EPYC 7713 processors (2.0 GHz, 64-cores/128-threads) with 2TB of NVMe SSD memory. We believe on average the runtime scales in an $O(n \log n)$ fashion, meaning that using only 50% of the samples would require considerably less than 50% of the runtime.

Instructions for running this code on a very small number of samples can by found in code/code_to_generate-atlas/exampleSetup.R, including expected output. Running the code in tmeMetaAnalysis_human.R on these samples takes approximately four hours.

---

## 📦 Installation

You can install the development version of `scTME` directly from GitHub using `devtools`:

```r
# install.packages("devtools")
devtools::install_github("michaeljanuszyk/scTME")

---

## 🚀 Getting Started & Pipeline Execution

To run this code in full, first download publicly available GEO or ArrayExpress datasets and convert each sample from each study into a separate Seurat object. A list of the datasets used in the accompanying manuscript can be found in `datasets_human.csv`. Save each object as the accession code + `.object.rds` (e.g., `GSM7114948.object.rds`) and save them into folders named by study into a parent folder called `readyForSeurat`. 

For example:
```text
readyForSeurat/
└── GSE190597/
    └── GSM7114948.object.rds
    


