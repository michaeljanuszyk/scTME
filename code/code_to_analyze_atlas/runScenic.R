# ==============================================================================
# Script: runScenic.R
# Description: Pipeline for SCENIC single-cell regulatory network analysis.
# ==============================================================================

# ------------------------------------------------------------------------------
# Required Libraries
# ------------------------------------------------------------------------------
library(Seurat)
library(SCENIC)
library(AUCell)
library(RcisTarget)
library(ggplot2)

# ==============================================================================
# SCENIC Pipeline Functions
# ==============================================================================

#' Run SCENIC Analysis Core
#'
#' @description Run SCENIC analysis on scRNA-seq Seurat object.
#' Note: Before running you must download the appropriate cisTarget databases 
#' (e.g., hg38 refseq-r80) and place them in your relative db_dir.
#'
#' @param obj A Seurat object containing scRNA-seq data.
#' @param org Character. Organism code (e.g., "hgnc").
#' @param db_dir Character. Relative path to cisTarget databases.
#' @param num_cores Numeric. Number of cores for parallel processing.
#' 
#' @return A Seurat object with SCENIC AUC values added as an assay.
#' @export
runScenicAnalysis <- function(obj, org = "hgnc", db_dir = "data/cisTarget", num_cores = 4) {
  
  exprMat <- as.matrix(GetAssayData(obj, assay = "RNA", slot = "counts"))
  scenicOptions <- initializeScenic(org = org, dbDir = db_dir, nCores = num_cores)
  
  genesKept <- geneFiltering(
    exprMat, 
    scenicOptions = scenicOptions, 
    minCountsPerGene = 3 * 0.01 * ncol(exprMat), 
    minSamples = ncol(exprMat) * 0.01
  )
  exprMat_filtered <- exprMat[genesKept, ]
  
  runCorrelation(exprMat_filtered, scenicOptions)
  exprMat_filtered_log <- log2(exprMat_filtered + 1)
  runGenie3(exprMat_filtered_log, scenicOptions)
  
  scenicOptions@settings$verbose <- TRUE
  scenicOptions <- runSCENIC_1_coexNetwork2modules(scenicOptions)
  scenicOptions <- runSCENIC_2_createRegulons(scenicOptions)
  scenicOptions <- runSCENIC_3_scoreCells(scenicOptions, exprMat_filtered_log)
  
  auc_matrix <- readRDS(getIntName(scenicOptions, "aucell_regulonAUC"))
  
  # Add the AUC matrix to the Seurat object as a new assay
  obj[['AUC']] <- CreateAssayObject(data = getAUC(auc_matrix))
  
  return(obj)
}


#' Generate SCENIC Feature Plots
#'
#' @description Make SCENIC plots placing cell-level regulon data into original object 
#' and produce the resulting FeaturePlots.
#'
#' @param objSce Seurat object output from runScenicAnalysis.
#' @param objAll Original Seurat object.
#' @param name Character. Prefix for saved plot files.
#' @param out_dir Character. Directory to save outputs.
#' @export
makeScenicPlots <- function(objSce, objAll, name, out_dir = "SCENIC") {
  
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  
  obj <- subset(objAll, cells = colnames(objAll)[colnames(objAll) %in% colnames(objSce)])
  
  objSce[['regulon']] <- objSce[['AUC']]
  obj[['regulon']] <- objSce[['regulon']]
  DefaultAssay(obj) <- 'regulon'
  
  for (i in 1:nrow(obj[['regulon']])) {
    regulon <- substr(rownames(obj[['regulon']])[i], 1, nchar(rownames(obj[['regulon']])[i]))
    print(regulon)
    
    p <- FeaturePlot(obj, regulon, min.cutoff = 'q05', max.cutoff = 'q95', alpha = 0.5) & 
         scale_color_viridis_c(option = "viridis", direction = 1)
    
    file_path <- file.path(out_dir, paste0(name, '_', substr(regulon, 1, nchar(regulon) - 3), '.jpg'))
    ggsave(file_path, plot = p, width = 8, height = 8, units = 'in', dpi = 300)
  }
}
