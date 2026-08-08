# Load required libraries
library(nichenetr)
library(Seurat)
library(tidyverse)

# Runs NicheNet on a Seurat object with specific sender and receiver cell types
#
# Before running, calculate the DEGs using Seurat's FindMarkers function. These will need to be
# calculated separately for each unique set of receiver cell types

runNicheNet (
  seurat_obj, # A Seurat object containig scRNA-seq data
  senders, # A character vector of cell types acting as senders
  receivers, # A character vector of cell types asking as receivers
  degs, # Differentially expressed genes among the receiver cell types
  cell_type_col, # Metadata column of Seurat object containing cell type annotations
  ligand_target_matrix_path = "nichenet_models/ligand_target_matrix.rds", # Path to the downloaded NicheNet file
  lr_network_path = "nichenet_models/lr_network.rds", # Path to the downloaded NicheNetFile
  blacklist = NULL, # Obtional vector of genes to exclude
  out_dir = "NicheNet_Output" # Directory path in which to save the resulting heatmap
) 


