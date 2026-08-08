# Load required libraries
library(infercnv)
library(Seurat)

# Runs inferCNV on a Seurat object with multiple parameters
# 
# Recommended downsampling Seurat object prior to running

runInferCnv (
  obj, # A Seurat object containing scRNA-seq data
  group.by = "cellType", # The metadata column containing cell type annotations
  ref_group_names = c("Normal"), # A set of annotations from cellType to use as a normal reference
  gene_order_file = "hg38_gencode_v27.txt", # Path to a gene positions file
  out_dir = "infercnv_output", # Directory path to save outputs and plots
  cutoff = 0.1, # Numeric cutoff for minimum gene expression
  num_threads = 4 # Number of cores to use for parallel processing
)


