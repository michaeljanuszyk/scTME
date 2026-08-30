# ==============================================================================
# Script: runNicheNet.R
# Description: Pipeline for predicting ligand-target links between cells.
# ==============================================================================

# ------------------------------------------------------------------------------
# Required Libraries
# ------------------------------------------------------------------------------
library(nichenetr)
library(Seurat)
library(tidyverse)
library(ggplot2)

# ==============================================================================
# NicheNet Pipeline Function
# ==============================================================================

#' Run NicheNet Cell-Cell Interaction Modeling
#'
#' @description Runs NicheNet on a Seurat object with specific sender and receiver cell types.
#' Before running, calculate the DEGs using Seurat's FindMarkers function. These will need to be
#' calculated separately for each unique set of receiver cell types.
#'
#' @param seurat_obj A Seurat object containing scRNA-seq data.
#' @param senders A character vector of cell types acting as senders.
#' @param receivers A character vector of cell types acting as receivers.
#' @param degs Differentially expressed genes among the receiver cell types.
#' @param cell_type_col Metadata column of Seurat object containing cell type annotations.
#' @param ligand_target_matrix_path Path to the downloaded NicheNet matrix file.
#' @param lr_network_path Path to the downloaded NicheNet network file.
#' @param blacklist Optional vector of genes to exclude.
#' @param out_dir Directory path in which to save the resulting heatmap.
#' 
#' @return A list containing top ligands, ligand activities, and the generated heatmap.
#' @export
runNicheNet <- function(
  seurat_obj, 
  senders, 
  receivers, 
  degs, 
  cell_type_col, 
  ligand_target_matrix_path = "nichenet_models/ligand_target_matrix.rds",
  lr_network_path = "nichenet_models/lr_network.rds",
  blacklist = NULL,
  out_dir = "NicheNet_Output"
) {
  
  Idents(seurat_obj) <- cell_type_col
  
  if (!file.exists(ligand_target_matrix_path) || !file.exists(lr_network_path)) {
    stop("NicheNet model files not found.")
  }
  
  message("Loading NicheNet prior models...")
  ligand_target_matrix <- readRDS(ligand_target_matrix_path)
  lr_network <- readRDS(lr_network_path)
  
  message(sprintf("Calculating expressed genes for receiver(s): %s", paste(receivers, collapse = ", ")))
  expressed_genes_receiver <- get_expressed_genes(receivers, seurat_obj, pct = 0.10)
  background_expressed_genes <- expressed_genes_receiver %>% .[. %in% rownames(ligand_target_matrix)]
  
  message(sprintf("Calculating expressed genes for sender(s): %s", paste(senders, collapse = ", ")))
  expressed_genes_sender <- get_expressed_genes(senders, seurat_obj, pct = 0.10)
  
  message("Filtering potential ligands and receptors...")
  ligands <- lr_network %>% pull(from) %>% unique()
  receptors <- lr_network %>% pull(to) %>% unique()
  
  expressed_ligands <- intersect(ligands, expressed_genes_sender)
  expressed_receptors <- intersect(receptors, expressed_genes_receiver)
  
  potential_ligands <- lr_network %>% 
    filter(from %in% expressed_ligands & to %in% expressed_receptors) %>% 
    pull(from) %>% unique()
  
  if (!is.null(blacklist)) {
    message("Removing blacklisted ambient RNA genes...")
    potential_ligands <- setdiff(potential_ligands, blacklist)
  }
  
  message("Predicting ligand activities based on provided DEGs...")
  ligand_activities <- predict_ligand_activities(
    geneset = degs, 
    background_expressed_genes = background_expressed_genes, 
    ligand_target_matrix = ligand_target_matrix, 
    potential_ligands = potential_ligands
  )
  
  ligand_activities <- ligand_activities %>% arrange(-pearson) %>% mutate(rank = rank(desc(pearson)))
  top_ligands <- ligand_activities %>% top_n(20, pearson) %>% pull(test_ligand) %>% unique()
  
  message("Inferring specific ligand-target links...")
  active_ligand_target_links_df <- top_ligands %>% 
    lapply(get_weighted_ligand_target_links, geneset = degs, ligand_target_matrix = ligand_target_matrix, n = 200) %>% 
    bind_rows() %>% drop_na()
  
  active_ligand_target_links <- prepare_ligand_target_visualization(
    ligand_target_df = active_ligand_target_links_df, 
    ligand_target_matrix = ligand_target_matrix, 
    cutoff = 0.33
  )
  
  message("Generating heatmap visualization...")
  order_ligands <- intersect(top_ligands, colnames(active_ligand_target_links)) %>% rev()
  order_targets <- active_ligand_target_links_df$target %>% unique() %>% intersect(rownames(active_ligand_target_links))
  vis_ligand_target <- active_ligand_target_links[order_targets, order_ligands] %>% t()
  
  p_ligand_target_network <- vis_ligand_target %>% 
    make_heatmap_ggplot(
      y_name = "Prioritized Sender Ligands", 
      x_name = "Predicted Target Genes", 
      color = "purple", legend_position = "top", 
      x_axis_position = "top", legend_title = "Regulatory potential"
    ) +
    theme(axis.text.x = element_text(angle = 45, hjust = 0, size = 6))
  
  message("Saving heatmap to output directory...")
  safe_senders <- paste(gsub("[^A-Za-z0-9]", "_", senders), collapse = "_")
  safe_receivers <- paste(gsub("[^A-Za-z0-9]", "_", receivers), collapse = "_")
  file_name <- paste0(safe_senders, "__", safe_receivers, "_.jpg")
  
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  file_path <- file.path(out_dir, file_name)
  
  ggsave(filename = file_path, plot = p_ligand_target_network, width = 14, height = 8, dpi = 300)
  message(sprintf("Saved plot successfully to: %s", file_path))
  
  return(list(top_ligands = top_ligands, ligand_activities = ligand_activities, heatmap = p_ligand_target_network))
}
