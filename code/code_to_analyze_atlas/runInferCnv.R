# ==============================================================================
# Script: runInferCnv.R
# Description: Pipeline for inferring copy number variations from scRNA-seq data.
# ==============================================================================

# ------------------------------------------------------------------------------
# Required Libraries
# ------------------------------------------------------------------------------
library(infercnv)
library(Seurat)

# ==============================================================================
# InferCNV Pipeline Function
# ==============================================================================

#' Run InferCNV Analysis
#'
#' @description Runs inferCNV on a Seurat object with multiple parameters. 
#' Recommended downsampling Seurat object prior to running.
#'
#' @param obj A Seurat object containing scRNA-seq data.
#' @param group.by Character. The metadata column containing cell type annotations. Defaults to "cellType".
#' @param ref_group_names Character vector. A set of annotations from cellType to use as a normal reference. Defaults to c("Normal").
#' @param gene_order_file Character. Path to a gene positions file. Defaults to "hg38_gencode_v27.txt".
#' @param out_dir Character. Directory path to save outputs and plots. Defaults to "infercnv_output".
#' @param cutoff Numeric. Numeric cutoff for minimum gene expression. Defaults to 0.1.
#' @param num_threads Numeric. Number of cores to use for parallel processing. Defaults to 4.
#' 
#' @return An InferCNV object.
#' @export
runInferCnv <- function(
  obj, 
  group.by = "cellType", 
  ref_group_names = c("Normal"), 
  gene_order_file = "hg38_gencode_v27.txt", 
  out_dir = "infercnv_output", 
  cutoff = 0.1,
  num_threads = 4
) {
  
  # Ensure output directory exists
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  
  # Extract raw counts as a sparse matrix
  counts_matrix <- GetAssayData(obj, assay = "RNA", slot = "counts")
  
  # Prepare the annotations dataframe
  if (group.by == "ident") {
    annotations <- data.frame(cell_type = Idents(obj), row.names = colnames(obj))
  } else {
    annotations <- data.frame(cell_type = obj[[group.by]][, 1], row.names = colnames(obj))
  }
  
  # Validate that the requested reference groups actually exist in the data
  missing_refs <- setdiff(ref_group_names, unique(annotations$cell_type))
  if (length(missing_refs) > 0) {
    warning(sprintf("The following ref_group_names were not found in the object annotations: %s", paste(missing_refs, collapse = ", ")))
  }
  
  # Initialize the InferCNV object
  infercnv_obj <- CreateInfercnvObject(
    raw_counts_matrix = counts_matrix,
    annotations_file = annotations,
    delim = "\t",
    gene_order_file = gene_order_file,
    ref_group_names = ref_group_names
  )
  
  # Run the core InferCNV pipeline
  infercnv_obj <- infercnv::run(
    infercnv_obj,
    cutoff = cutoff, 
    out_dir = out_dir,
    cluster_by_groups = TRUE,
    denoise = TRUE,
    HMM = TRUE,
    num_threads = num_threads
  )
  
  return(infercnv_obj)
}
