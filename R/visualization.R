#' Generate a Confusion Matrix Heatmap for scTME Predictions
#'
#' This function evaluates cell type predictions against ground truth 
#' annotations by generating a row-normalized confusion matrix heatmap. 
#' Color intensity represents the proportion of predictions, while the 
#' printed numbers display the raw cell counts.
#'
#' @param dataset A data frame containing the prediction and annotation metadata.
#' @param pred_col A character string specifying the column name for predictions. Default is "prediction".
#' @param annot_col A character string specifying the column name for true annotations. Default is "lv1_annot".
#' @param title A character string for the main title of the heatmap.
#'
#' @return A pheatmap object.
#' @importFrom pheatmap pheatmap
#' @importFrom grDevices colorRampPalette
#' @export
scTME_generateHeatmap <- function(dataset, 
                                  pred_col = "prediction", 
                                  annot_col = "lv1_annot",
                                  title = "scTME Predictions vs. True Annotations") {
  
  # 1. Create the base matrix dynamically from the specified columns
  conf_table <- table(dataset[[pred_col]], dataset[[annot_col]])
  conf_matrix <- as.matrix(conf_table)
  
  # 2. Row Normalization (convert to proportions from 0.0 to 1.0)
  # A safeguard is included to prevent division by zero for empty rows
  row_sums <- rowSums(conf_matrix)
  row_sums[row_sums == 0] <- 1 
  norm_matrix <- conf_matrix / row_sums
  
  # 3. Create a custom palette (Pure white for 0, then noticeable pink to dark red)
  my_colors <- c("white", colorRampPalette(c("lightpink", "firebrick3"))(49))
  
  # 4. Generate the heatmap
  heatmap_plot <- pheatmap::pheatmap(
    mat = norm_matrix, 
    display_numbers = conf_matrix, 
    cluster_rows = FALSE,
    cluster_cols = FALSE,
    color = my_colors,
    main = title
  )
  
  return(heatmap_plot)
}


#' Plot scTME Predictions on a UMAP
#'
#' @param seurat_obj A Seurat object containing scTME predictions in the metadata.
#' @param reduction A character string specifying the dimensional reduction to use. Default is "umap".
#' @param pred_col A character string specifying the metadata column containing predictions.
#' @export
plotPredictionUmap <- function(seurat_obj,
                               reduction = "umap",
                               pred_col = "prediction") {
 
  # 1. Verify the reduction exists using Seurat's official accessor
  if (!reduction %in% SeuratObject::Reductions(seurat_obj)) {
    message(paste("Notice: The reduction", reduction, "was not found. Skipping plot."))
    return(invisible(NULL))
  }
 
  # 2. Generate the plot
  plot <- Seurat::DimPlot(seurat_obj,
                          reduction = reduction,
                          group.by = pred_col,
                          label = TRUE) +
          ggplot2::ggtitle("scTME Automated Predictions")

  return(plot)
}
