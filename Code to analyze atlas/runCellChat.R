# Load required libraries
library(Seurat)
library(CellChat)
library(patchwork)
library(ggplot2)

# Runs CellChat on a Seurat object with multiple customization parameters
# 
# Before running, ensure the Seurat object contains only the specific cell types you wish to
# evaluate. Recommended downsampling rate is 1,000 cells per cell type

runCellChat(
  obj, # A Seurat object containing scRNA-seq data
  group.by = "ident", # The name of the metadata column containing cell type labels
  org = "human", # Organism (only supports human and mouse)
  db.use = "Secreted Signaling", # Databases to use ("Secreted Signaling", "ECM-Receptor", "Cell-Cell Contact", or "all")
  out_dir = "CellChat_Output", # Directory to save output files
  min.cells = 10, # Minimum cell number required in a cluster for communication analysis
  pathways_to_show = NULL # Optional. Specific pathways for which to generate detailed plots (e.g., c("CXCL", "CCL"))
) 

