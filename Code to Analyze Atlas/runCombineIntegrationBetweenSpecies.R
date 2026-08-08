# Load required libraries
library(Seurat)
library(nichenetr)

# Runs the combineIntegrationBetweenSpecies, which takes a human and a mouse integrated (from meta-analysis) 
# Seurat object and combines them into a single Seurat object human-space using human orthologs for mouse genes. 
# This is most effective when using the target_celltype parameter to combine objects one cell type at a time. 

runCombineIntegrationBetweenSpecies (
  obj_human, # An integrated Seurat object containing human scRNA-seq data
  obj_mouse, # An integrated Seurat object containing mouse scRNA-seq data
  target_celltype = NULL, # An optional cell type to subset on (e.g., "Fibroblast")
  group.by = "cellType" # A metadata column containing cell annotations
)
