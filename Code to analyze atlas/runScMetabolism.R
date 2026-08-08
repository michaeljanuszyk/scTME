# Load required libraries
library(scMetabolism)
library(ggplot2)
library(rsvd)
library(VISION)

# Run scMetabolism on Seurat object with cell type names in the cellType_fine_abbr slot
# Optional argument subset will focus only on those cell types listed

runScMetabolism(
  obj, # A Seurat object containing scRNA-seq data
  subset = NULL, # An optional vector of cell types (from cellType_fine_abbr) to use
  name = "allCells", # A name to post-pend to saved file names
  input.pathway = c("Glycolysis / Gluconeogenesis", "Oxidative phosphorylation", "Citrate cycle (TCA cycle)")
  # A list of native metabolic pathways provided by scMetabolism
)

