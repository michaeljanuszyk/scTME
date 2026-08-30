# ==============================================================================
# Script: runScMetabolism.R
# Description: Pipeline for single-cell metabolism pathway analysis.
# ==============================================================================

# ------------------------------------------------------------------------------
# Required Libraries
# ------------------------------------------------------------------------------
library(scMetabolism)
library(ggplot2)
library(rsvd)
library(VISION)
library(Seurat)

# ==============================================================================
# scMetabolism Pipeline Function
# ==============================================================================

#' Run scMetabolism Pathway Analysis
#'
#' @description Run scMetabolism on Seurat object with cell type names in the 
#' cellType_fine_abbr slot. Optional argument subset will focus only on those 
#' cell types listed.
#'
#' @param obj A Seurat object containing scRNA-seq data.
#' @param subset An optional vector of cell types (from cellType_fine_abbr) to use.
#' @param name A name to post-pend to saved file names.
#' @param input.pathway A list of native metabolic pathways provided by scMetabolism.
#' @param out_dir Directory path to save outputs.
#' @export
runScMetabolism <- function(
  obj, 
  subset = NULL, 
  name = "allCells", 
  input.pathway = c("Glycolysis / Gluconeogenesis", "Oxidative phosphorylation", "Citrate cycle (TCA cycle)"),
  out_dir = "figures"
) {
  
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  
  # Check if subset is provided
  if (!is.null(subset)) {
    obj <- subset(obj, cellType_fine_abbr %in% subset)
  }
  
  # Convert to v3 Seurat object
  obj_v3 <- obj
  obj_v3[['RNA']] <- as(obj[['RNA']], 'Assay')
  
  # Run scMetabolism
  countexp.Seurat <- sc.metabolism.Seurat(
    obj = obj_v3, 
    method = "AUCell", 
    imputation = FALSE, 
    ncores = 2, 
    metabolism.type = "KEGG"
  )
  
  countexp.Seurat_VISION <- sc.metabolism.Seurat(
    obj = obj_v3, 
    method = "VISION", 
    imputation = FALSE, 
    ncores = 2, 
    metabolism.type = "KEGG"
  )
  
  metabolism.matrix <- countexp.Seurat@assays$METABOLISM$score
  
  # Generate plot
  p <- DotPlot.metabolism(
    obj = countexp.Seurat, 
    pathway = input.pathway, 
    phenotype = "cellType_fine", 
    norm = "y"
  )
  
  # Save plot
  ggsave(file.path(out_dir, paste0("metabolic_fine_", name, ".jpg")), plot = p, width = 12, height = 6)
}
