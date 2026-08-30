# ==============================================================================
# Script: runCombineIntegrationBetweenSpecies.R
# Description: Merges human and mouse single-cell datasets via ortholog conversion.
# ==============================================================================

# ------------------------------------------------------------------------------
# Required Libraries
# ------------------------------------------------------------------------------
library(Seurat)
library(nichenetr)

# ==============================================================================
# Cross-Species Integration Function
# ==============================================================================

#' Combine Integration Between Species
#'
#' @description Takes a human and a mouse integrated (from meta-analysis) Seurat object
#' and combines them into a single Seurat object in human-space using human orthologs 
#' for mouse genes. This is most effective when using the target_celltype 
#' parameter to combine objects one cell type at a time.
#'
#' @param obj_human An integrated Seurat object containing human scRNA-seq data.
#' @param obj_mouse An integrated Seurat object containing mouse scRNA-seq data.
#' @param target_celltype Character. An optional cell type to subset on (e.g., "Fibroblast"). Defaults to NULL.
#' @param group.by Character. A metadata column containing cell annotations. Defaults to "cellType".
#' 
#' @return A merged Seurat object with mouse genes converted to human orthologs.
#' @export
runCombineIntegrationBetweenSpecies <- function(
  obj_human, 
  obj_mouse, 
  target_celltype = NULL, 
  group.by = "cellType"
) {
  
  # Subset objects if a target cell type is specified
  if (!is.null(target_celltype)) {
    message(sprintf("Subsetting data to target cell type: %s", target_celltype))
    
    # Extract metadata columns for subsetting
    human_labels <- obj_human[[group.by]][, 1]
    mouse_labels <- obj_mouse[[group.by]][, 1]
    
    obj_human <- subset(obj_human, cells = colnames(obj_human)[which(human_labels == target_celltype)])
    obj_mouse <- subset(obj_mouse, cells = colnames(obj_mouse)[which(mouse_labels == target_celltype)])
  }
  
  if (ncol(obj_human) == 0 || ncol(obj_mouse) == 0) {
    stop("One or both objects have 0 cells after subsetting.")
  }
  
  # Convert mouse gene symbols to human orthologs
  message("Converting mouse genes to human orthologs...")
  
  # Get raw counts from the mouse object
  mouse_counts <- GetAssayData(obj_mouse, assay = "RNA", slot = "counts")
  mouse_genes <- rownames(mouse_counts)
  
  # Convert symbols using nichenetr
  human_orthologs <- convert_mouse_to_human_symbols(mouse_genes)
  
  # Filter out genes that did not map (NA values)
  valid_mapping <- !is.na(human_orthologs)
  mouse_counts <- mouse_counts[valid_mapping, ]
  human_orthologs <- human_orthologs[valid_mapping]
  
  # Apply the new human names to the mouse counts matrix
  rownames(mouse_counts) <- human_orthologs
  
  # Remove duplicated gene names that can arise from ortholog mapping (many-to-one)
  unique_genes <- !duplicated(rownames(mouse_counts))
  mouse_counts <- mouse_counts[unique_genes, ]
  
  # Recreate the mouse object with human gene names, preserving original metadata
  obj_mouse_converted <- CreateSeuratObject(counts = mouse_counts, meta.data = obj_mouse@meta.data)
  
  # Find common genes between the two datasets
  message("Identifying overlapping genes...")
  common_genes <- intersect(rownames(obj_human), rownames(obj_mouse_converted))
  message(sprintf("Found %d common genes.", length(common_genes)))
  
  # Subset both objects to only the common genes
  human_counts_final <- GetAssayData(obj_human, assay = "RNA", slot = "counts")[common_genes, ]
  mouse_counts_final <- GetAssayData(obj_mouse_converted, assay = "RNA", slot = "counts")[common_genes, ]
  
  obj_human_final <- CreateSeuratObject(counts = human_counts_final, meta.data = obj_human@meta.data)
  obj_mouse_final <- CreateSeuratObject(counts = mouse_counts_final, meta.data = obj_mouse_converted@meta.data)
  
  # Add species metadata for downstream identification
  obj_human_final$species <- "human"
  obj_mouse_final$species <- "mouse"
  
  if ("study" %in% colnames(obj_human_final@meta.data)) {
    obj_human_final$speciesStudy <- paste0("human_", obj_human_final$study)
  }
  if ("study" %in% colnames(obj_mouse_final@meta.data)) {
    obj_mouse_final$speciesStudy <- paste0("mouse_", obj_mouse_final$study)
  }
  
  # Merge objects
  message("Merging human and mouse objects...")
  merged_obj <- merge(obj_human_final, y = obj_mouse_final)
  
  if (packageVersion("Seurat") >= "5.0.0") {
    merged_obj <- JoinLayers(merged_obj)
  }
  
  message("Integration complete.")
  return(merged_obj)
}
