# ==============================================================================
# Script: runCellChat.R
# Description: Modular pipeline for executing CellChat intercellular communication.
# ==============================================================================

# ------------------------------------------------------------------------------
# Required Libraries
# ------------------------------------------------------------------------------
library(Seurat)
library(CellChat)
library(patchwork)
library(ggplot2)

# ==============================================================================
# CellChat Pipeline Function
# ==============================================================================

#' Run CellChat Intercellular Communication Pipeline
#'
#' @description Runs CellChat on a Seurat object with multiple customization parameters.
#' Before running, ensure the Seurat object contains only the specific cell types you wish to
#' evaluate. Recommended downsampling rate is 1,000 cells per cell type.
#'
#' @param obj A Seurat object containing scRNA-seq data.
#' @param group.by Character. The name of the metadata column containing cell type labels. Defaults to "ident".
#' @param org Character. Organism (only supports "human" and "mouse"). Defaults to "human".
#' @param db.use Character. Databases to use ("Secreted Signaling", "ECM-Receptor", "Cell-Cell Contact", or "all"). Defaults to "Secreted Signaling".
#' @param out_dir Character. Directory to save output files. Defaults to "CellChat_Output".
#' @param min.cells Numeric. Minimum cell number required in a cluster for communication analysis. Defaults to 10.
#' @param pathways_to_show Character vector. Optional. Specific pathways for which to generate detailed plots (e.g., c("CXCL", "CCL")).
#' 
#' @return A processed CellChat object.
#' @export
runCellChat <- function(
  obj, 
  group.by = "ident", 
  org = "human", 
  db.use = "Secreted Signaling", 
  out_dir = "CellChat_Output", 
  min.cells = 10, 
  pathways_to_show = NULL
) {
  
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  
  # Set cell identities
  obj$labels <- if (group.by == "ident") Idents(obj) else obj[[group.by]][, 1]
  obj$labels <- droplevels(as.factor(obj$labels))
  
  cellchat <- createCellChat(object = obj, group.by = "labels")
  cellchat <- setIdent(cellchat, ident.use = "labels")
  groupSize <- as.numeric(table(cellchat@idents))
  
  # Configure database
  if (org == "human") {
    CellChatDB <- CellChatDB.human
    ppi_data <- PPI.human
  } else if (org == "mouse") {
    CellChatDB <- CellChatDB.mouse
    ppi_data <- PPI.mouse
  } else stop("Unsupported organism. Please specify 'human' or 'mouse'.")
  
  CellChatDB.use <- if (db.use != "all") subsetDB(CellChatDB, search = db.use) else CellChatDB
  cellchat@DB <- CellChatDB.use
  
  # Preprocessing and Pipeline execution
  cellchat <- subsetData(cellchat)
  cellchat <- identifyOverExpressedGenes(cellchat)
  cellchat <- identifyOverExpressedInteractions(cellchat)
  cellchat <- projectData(cellchat, ppi_data)
  cellchat <- computeCommunProb(cellchat, raw.use = TRUE)
  cellchat <- filterCommunication(cellchat, min.cells = min.cells)
  cellchat <- computeCommunProbPathway(cellchat)
  cellchat <- aggregateNet(cellchat)
  
  # Save object and broad network plots
  saveRDS(cellchat, file = file.path(out_dir, "cellchat_analyzed.rds"), compress = FALSE)
  
  jpeg(file.path(out_dir, 'netVisual_circle_all_counts.jpg'), width = 4600, height = 2875, res = 300)
  netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = TRUE, label.edge = FALSE, title.name = "Number of interactions")
  dev.off()
  
  jpeg(file.path(out_dir, 'netVisual_circle_all_weights.jpg'), width = 4600, height = 2875, res = 300)
  netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = TRUE, label.edge = FALSE, title.name = "Interaction weights/strength")
  dev.off()
  
  # Optional Pathway-Specific Plots
  if (!is.null(pathways_to_show)) {
    cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")
    available_pathways <- cellchat@netP$pathways
    
    for (pathway in pathways_to_show) {
      if (pathway %in% available_pathways) {
        jpeg(file.path(out_dir, paste0('Circle_', pathway, '.jpg')), width = 4600, height = 2875, res = 300)
        netVisual_aggregate(cellchat, signaling = pathway, layout = "circle")
        dev.off()
        
        jpeg(file.path(out_dir, paste0('Chord_', pathway, '.jpg')), width = 4600, height = 2875, res = 300)
        netVisual_aggregate(cellchat, signaling = pathway, layout = "chord")
        dev.off()
        
        gg <- netAnalysis_contribution(cellchat, signaling = pathway)
        ggsave(filename = file.path(out_dir, paste0("Contribution_", pathway, "_LR.pdf")), plot = gg, width = 5, height = 3, units = 'in', dpi = 300)
        
        jpeg(file.path(out_dir, paste0('RoleNetwork_', pathway, '.jpg')), width = 4600, height = 2875, res = 300)
        netAnalysis_signalingRole_network(cellchat, signaling = pathway, width = 12, height = 4, font.size = 10)
        dev.off()
      } else {
        warning(paste("Pathway", pathway, "was not found as significant. Skipping plots."))
      }
    }
  }
  return(cellchat)
}
