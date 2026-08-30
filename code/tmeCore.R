# ==============================================================================
# Script: tmeCore.R
# Description: Core single-cell analytical and plotting functions for scTME.
# ==============================================================================

# ------------------------------------------------------------------------------
# Required Libraries
# ------------------------------------------------------------------------------
library(Seurat)
library(ggplot2)
library(viridis)
library(scCustomize)
library(patchwork)
library(scales)
library(stringr)

# ==============================================================================
# Marker & Plotting Functions
# ==============================================================================

#' Get Top Markers
#'
#' @param obj A Seurat object.
#' @param only.pos Logical. Only return positive markers.
#' @param assay Character. Assay to use.
#' @param out_dir Character. Directory to save outputs.
#' @export
getTopMarkers <- function(obj, only.pos = TRUE, assay = "RNA", out_dir = "results") {
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  
  if (packageVersion("Seurat") >= "5.0.0") obj <- JoinLayers(obj)
  obj[[assay]]$data <- as(object = obj[[assay]]$data, Class = "dgCMatrix")
  Idents(obj) <- factor(x = Idents(obj), levels = sort(levels(obj), decreasing = FALSE))
  
  pbmc.markers <- FindAllMarkers(obj, only.pos = only.pos, min.pct = -1.10, logfc.threshold = 0.25, max.cells.per.ident = 5000)
  write.csv(pbmc.markers, file.path(out_dir, "pbmc_markers.csv"), quote = FALSE)
  
  return(obj)
}


#' Generate Dimensionality Reduction Plot
#'
#' @param obj A Seurat object.
#' @param raster Logical. Whether to rasterize the plot.
#' @param out_dir Character. Output directory.
#' @export
makeDimPlot <- function(obj, raster = TRUE, out_dir = "results") {
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  
  DimPlot(obj, group.by = 'cellType', raster = raster, label = FALSE)
  ggsave(file.path(out_dir, "featurePlot.jpg"), height = 8, width = 9, limitsize = FALSE)
}


#' Generate Customized DotPlot
#'
#' @param obj A Seurat object.
#' @param features Vector of genes.
#' @param dotColor Color scheme name.
#' @param w Numeric. Width of plot.
#' @param h Numeric. Height of plot.
#' @param out_dir Character. Output directory.
#' @export
makeDotPlot <- function(obj, features, dotColor = 'black', w = 6, h = 6, out_dir = "results") {
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  
  Idents(obj) <- factor(x = Idents(obj), levels = sort(levels(obj), decreasing = TRUE))
  
  lowHex <- switch(dotColor,
    'black'  = "#DAE2F0", 'blue'   = "#EEF4FB", 'green'  = "#E2F6EC",
    'orange' = "#FAF0EA", 'pink'   = "#FEEFF9", 'purple' = "#F6EFFC",
    'yellow' = "#F9F7D7"
  )
  highHex <- switch(dotColor,
    'black'  = "#496079", 'blue'   = "#9FBCE1", 'green'  = "#A5DEC9",
    'orange' = "#F5D5C0", 'pink'   = "#E5BCDA", 'purple' = "#BCAADD",
    'yellow' = "#F7E97D"
  )

  p <- DotPlot(obj, features = features, scale.min = 0, scale.max = 100) +
    NoLegend() + scale_colour_gradient(low = lowHex, high = highHex) +
    theme(
      axis.title.x = element_blank(), axis.title.y = element_blank(),
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
      axis.text.y = element_text(hjust = 0), axis.line = element_blank(),
      panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.75)
    )
  
  ggsave(file.path(out_dir, "featurePlot.jpg"), height = h, width = w, limitsize = FALSE)
  return(p)
}

# ==============================================================================
# Analytical & Pipeline Workflows
# ==============================================================================

#' Standard Seurat Analytical Pipeline
#'
#' @param obj A Seurat object.
#' @param name Character. Prefix for saving intermediate steps.
#' @param skipStudies Vector of studies to exclude.
#' @param out_dir Character. Output directory for intermediate objects.
#' @export
analyzeObj <- function(obj, name, skipStudies = c(), out_dir = "sketch_v5_all") {
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  if (length(Layers(obj)) > 3) obj <- JoinLayers(obj)
  
  temp <- table(obj$study)
  badStudies <- names(temp)[temp <= 9]
  `%notin%` <- Negate(`%in%`)
  cellsToKeep <- colnames(obj)[obj$study %notin% c(badStudies, skipStudies)]
  
  obj <- subset(obj, cells = cellsToKeep)
  obj@reductions <- list()
  DefaultAssay(obj) <- "RNA"
  
  if (length(Layers(obj)) > 1) obj[["RNA"]]$data <- NULL
  obj[["sketch"]] <- NULL
  obj$leverage.score <- NULL
  obj$umap1 <- NULL
  obj$umap2 <- NULL
  
  counts.mat <- as(object = obj[['RNA']]$counts, Class = "dgCMatrix")
  obj <- CreateSeuratObject(counts = counts.mat, meta.data = obj@meta.data)
  
  sketch_dir <- file.path(out_dir, paste0(name, "_sketch1000"))
  write_matrix_dir(mat = obj[["RNA"]]$counts, dir = sketch_dir, overwrite = TRUE)
  counts.mat <- open_matrix_dir(dir = sketch_dir)
  counts.mat <- as(counts.mat, Class = "dgCMatrix")
  obj[["RNA"]]$counts <- counts.mat
  saveRDS(obj, file.path(out_dir, paste0(name, "_sketch1000.rds")), compress = FALSE)

  options(future.globals.maxSize = 10 * 1024^3)
  obj <- NormalizeData(obj)
  obj[["RNA"]] <- split(obj[["RNA"]], f = obj$study)
  obj <- FindVariableFeatures(obj, verbose = TRUE)
  obj <- ScaleData(obj, verbose = TRUE)
  obj <- RunPCA(obj, npcs = 30, verbose = TRUE)
  obj <- IntegrateLayers(obj, method = HarmonyIntegration, orig.reduction = "pca", new.reduction = "harmony", kmeans_init_nstart = 20, kmeans_init_iter_max = 5000, verbose = TRUE)
  obj <- RunUMAP(obj, reduction = "harmony", dims = 1:30, reduction.name = "umap.harmony", return.model = TRUE, min.dist = 0.001)
  obj <- FindNeighbors(obj, reduction = "harmony", dims = 1:30)
  
  for (i in c(0.06, 0.07, 0.09, 0.1, 0.12, 0.14, 0.16, 0.18, 0.2, 0.3, 0.4, 0.5)) {
    obj <- FindClusters(obj, resolution = i)
    DimPlot(obj, raster = TRUE, label = TRUE)
    ggsave(file.path(out_dir, paste0(name, "_counts_obj_findvar_harmony_umap_res", i, ".jpg")), width = 5, height = 5, units = "in", limitsize = FALSE)
    saveRDS(obj, file.path(out_dir, paste0(name, "_counts_obj_findvar_harmony_umap_res", i, ".rds")), compress = FALSE)
  }
  return(obj)
}


#' Generate Main Figure Panels
#'
#' @param obj A Seurat object.
#' @param name Character. Plot output name.
#' @param genes Vector of genes for DotPlot.
#' @param lowHex Character. Low color hex.
#' @param highHex Character. High color hex.
#' @param out_dir Character. Output directory.
#' @export
makeMainPanel <- function(obj, name, genes, lowHex = '#EEF4FB', highHex = '#9FBCE1', out_dir = "sketch_v5_all/final_figures") {
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  pal <- colorRampPalette(c("#7C1E6F", "#CF597F", "#F0746E", "#EEB57A", "#FDDE9C", "#C6E174", "#9CCB85", "#7CCBA1", "#3AB185", "#029099", "#12739D", "#055275"))  

  Idents(obj) <- factor(obj$cellType, levels = rownames(table(obj$cellType)))
  
  DimPlot(obj, reduction = "umap.harmony", raster = TRUE, cols = pal(length(unique(obj$cellType)))) + theme(plot.title = element_blank())
  ggsave(file.path(out_dir, paste0(name, "_legend.jpg")), width = 4, height = 3, units = "in", limitsize = FALSE)
  
  DimPlot(obj, reduction = "umap.harmony", raster = TRUE, cols = pal(length(unique(obj$cellType)))) + NoLegend() + theme(plot.title = element_blank())
  ggsave(file.path(out_dir, paste0(name, ".jpg")), width = 6, height = 6, units = "in", limitsize = FALSE)

  Idents(obj) <- factor(Idents(obj), levels = rev(levels(obj)))
  d <- DotPlot(obj, features = genes, scale.min = 0, scale.max = 100)
  d + NoLegend() + scale_colour_gradient(low = lowHex, high = highHex) +
    theme(
      axis.title.x = element_blank(), axis.title.y = element_blank(),
      axis.text.x = element_text(angle = 90, hjust = 1),
      axis.text.y = element_text(hjust = 0), axis.line = element_blank(),
      panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.75)
    )
  ggsave(file.path(out_dir, paste0(name, "_dotplot.jpg")), width = 3.5, height = 3.5, units = "in")
}


#' Generate Supplementary Figure Panels
#'
#' @export
makeSuppPanel <- function(obj, name, genes, neg, dotColor = 'black', rasterFeaturePlot = TRUE, out_dir = "sketch_v5_all/final_figures") {
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  pal <- colorRampPalette(c("#7C1E6F", "#CF597F", "#F0746E", "#EEB57A", "#FDDE9C", "#C6E174", "#9CCB85", "#7CCBA1", "#3AB185", "#029099", "#12739D", "#055275"))
  
  lowHex <- switch(dotColor,
    'black' = "#DAE2F0", 'blue' = "#EEF4FB", 'green' = "#E2F6EC",
    'orange'= "#FAF0EA", 'pink' = "#FEEFF9", 'purple'= "#F6EFFC", 'yellow'= "#F9F7D7"
  )
  highHex <- switch(dotColor,
    'black' = "#496079", 'blue' = "#9FBCE1", 'green' = "#A5DEC9",
    'orange'= "#F5D5C0", 'pink' = "#E5BCDA", 'purple'= "#BCAADD", 'yellow'= "#F7E97D"
  )

  obj$cellType_specific <- obj$cellType
  Idents(obj) <- factor(obj$cellType_specific, levels = rownames(table(obj$cellType_specific)))
  
  DimPlot(obj, reduction = "umap.harmony", raster = TRUE, cols = pal(length(unique(obj$cellType_specific)))) + theme(plot.title = element_blank())
  ggsave(file.path(out_dir, paste0(name, "_specific_legend.jpg")), width = 6, height = 6, units = "in", limitsize = FALSE)
  
  DimPlot(obj, reduction = "umap.harmony", raster = rasterFeaturePlot, cols = pal(length(unique(obj$cellType_specific)))) + NoLegend() + theme(plot.title = element_blank())
  ggsave(file.path(out_dir, paste0(name, "_specific.jpg")), width = 6, height = 6, units = "in", limitsize = FALSE)
  
  Idents(obj) <- factor(Idents(obj), levels = rev(levels(obj)))
  d <- DotPlot(obj, features = unique(c(genes, neg)), scale.min = 0, scale.max = 100)
  d + NoLegend() + scale_colour_gradient(low = lowHex, high = highHex) +
    theme(
      axis.title.x = element_blank(), axis.title.y = element_blank(),
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
      axis.text.y = element_text(hjust = 0), axis.line = element_blank(),
      panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.75)
    )
  ggsave(file.path(out_dir, paste0(name, "_specific_dotplot.jpg")), width = 7, height = 3.5, units = "in")

  vpal <- viridis(n = 20, option = "D", direction = 1)
  f <- FeaturePlot_scCustom(obj, features = genes, raster = rasterFeaturePlot, colors_use = vpal, order = FALSE, num_columns = 3, max.cutoff = 'q95') & theme(plot.title = element_text(size = 24))
  ggsave(plot = f, file.path(out_dir, paste0(name, "_specific_featureplot.jpg")), width = 5 * 2.7, height = 5 * 4.8, units = "in")

  Idents(obj) <- factor(x = Idents(obj), levels = sort(levels(obj), decreasing = FALSE))
  cpal <- pal(length(unique(obj$cellType_specific)))
  v <- Stacked_VlnPlot(obj, features = genes, colors_use = cpal, x_lab_rotate = 45)
  ggsave(plot = v, file.path(out_dir, paste0(name, "_specific_vlnplot.jpg")), width = 5, height = 10.5, units = "in")
}


#' Generate Proportion Bar Graph
#'
#' @export
makeBarGraph <- function(obj, name, type, excludeCellType = 'NULL', height = 4, width = 8, out_dir = "sketch_v5_all/final_figures") {
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  pal <- colorRampPalette(c("#7C1E6F", "#CF597F", "#F0746E", "#EEB57A", "#FDDE9C", "#C6E174", "#9CCB85", "#7CCBA1", "#3AB185", "#029099", "#12739D", "#055275"))
  color <- pal(length(unique(obj$cellType[obj$cellType != excludeCellType])))
  
  p <- ggplot(obj@meta.data[obj$cellType != excludeCellType, ], aes(x = type[obj$cellType != excludeCellType], fill = cellType)) +
    labs(x = "Cancer Type", y = "Proportion (%)", fill = "Subtype") +
    geom_bar(position = "fill") +
    scale_fill_manual(values = color) +
    theme(
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, color = "black", size = 12),
      axis.title.x = element_blank(), axis.text.y = element_blank(),
      axis.ticks.y = element_blank(), axis.title.y = element_text(color = "black"),
      panel.background = element_blank(), plot.background = element_blank(),
      axis.line.x = element_line(color = "black")
    ) +
    scale_x_discrete(expand = c(0, 0)) +
    scale_y_continuous(labels = percent_format(scale = 100), limits = c(0, 1))

  ggsave(file.path(out_dir, paste0(name, "_combined_bargraph.jpg")), plot = p, width = width, height = height)
}


#' Generate Standardized Final Figures
#'
#' @export
makeFinalFigures <- function(obj, name, genes, dotColor = 'black', out_dir = "sketch_v5_all/final_figures") {
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  pal <- colorRampPalette(c("#7C1E6F", "#CF597F", "#F0746E", "#EEB57A", "#FDDE9C", "#C6E174", "#9CCB85", "#7CCBA1", "#3AB185", "#029099", "#12739D", "#055275")) 
  
  lowHex <- switch(dotColor,
    'black' = "#DAE2F0", 'blue' = "#EEF4FB", 'green' = "#E2F6EC",
    'orange'= "#FAF0EA", 'pink' = "#FEEFF9", 'purple'= "#F6EFFC", 'yellow'= "#F9F7D7"
  )
  highHex <- switch(dotColor,
    'black' = "#496079", 'blue' = "#9FBCE1", 'green' = "#A5DEC9",
    'orange'= "#F5D5C0", 'pink' = "#E5BCDA", 'purple'= "#BCAADD", 'yellow'= "#F7E97D"
  )

  Idents(obj) <- factor(obj$cellType, levels = rownames(table(obj$cellType)))
  
  DimPlot(obj, raster = FALSE, cols = pal(length(unique(obj$cellType)))) + theme(plot.title = element_blank())
  ggsave(file.path(out_dir, paste0(name, "_legend.jpg")), width = 4, height = 3, units = "in", limitsize = FALSE)
  
  DimPlot(obj, raster = FALSE, cols = pal(length(unique(obj$cellType)))) + NoLegend() + theme(plot.title = element_blank())
  ggsave(file.path(out_dir, paste0(name, ".jpg")), width = 6, height = 6, units = "in", limitsize = FALSE)

  Idents(obj) <- factor(Idents(obj), levels = rev(sort(unique(Idents(obj)))))
  d <- DotPlot(obj, features = genes, scale.min = 0, scale.max = 100)
  d + NoLegend() + scale_colour_gradient(low = lowHex, high = highHex) +  
    theme(
      axis.title.x = element_blank(), axis.title.y = element_blank(),
      axis.text.x = element_text(angle = 90, hjust = 1),
      axis.text.y = element_text(hjust = 0), axis.line = element_blank(),
      panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.75)
    )
  ggsave(file.path(out_dir, paste0(name, "_dotplot.jpg")), width = 3.5, height = 3.5, units = "in")
}


#' Generate Final Supplementary Figures
#'
#' @export
makeSuppFinal <- function(obj, prefix, dotColor, genes, extraGenes, out_dir = "sketch_v5_all/final_figures") {
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  dotplot_genes <- c(genes, extraGenes)
  pal <- colorRampPalette(c("#7C1E6F", "#CF597F", "#F0746E", "#EEB57A", "#FDDE9C", "#C6E174", "#9CCB85", "#7CCBA1", "#3AB185", "#029099", "#12739D", "#055275"))

  Idents(obj) <- factor(obj$cellType, levels = rownames(table(obj$cellType)))

  DimPlot(obj, reduction = "umap.harmony", raster = TRUE, cols = pal(length(unique(obj$cellType)))) + theme(plot.title = element_blank()) 
  ggsave(file.path(out_dir, paste0(prefix, "_legend.jpg")), width = 6, height = 6, units = "in", limitsize = FALSE)

  DimPlot(obj, reduction = "umap.harmony", raster = TRUE, cols = pal(length(unique(obj$cellType)))) + theme(plot.title = element_blank()) + NoLegend()
  ggsave(file.path(out_dir, paste0(prefix, ".jpg")), width = 6, height = 6, units = "in", limitsize = FALSE)

  lowHex <- switch(dotColor,
    'blue' = "#EEF4FB", 'green' = "#E2F6EC", 'orange'= "#FAF0EA",
    'pink' = "#FEEFF9", 'purple'= "#F6EFFC", 'yellow'= "#F9F7D7"
  )
  highHex <- switch(dotColor,
    'blue' = "#9FBCE1", 'green' = "#A5DEC9", 'orange'= "#F5D5C0",
    'pink' = "#E5BCDA", 'purple'= "#BCAADD", 'yellow'= "#F7E97D"
  )
  
  Idents(obj) <- factor(Idents(obj), levels = rev(levels(obj)))
  d <- DotPlot(obj, features = dotplot_genes, scale.min = 0, scale.max = 100)
  d + NoLegend() + scale_colour_gradient(low = lowHex, high = highHex) +
    theme(
      axis.title.x = element_blank(), axis.title.y = element_blank(),
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
      axis.text.y = element_text(hjust = 0), axis.line = element_blank(),
      panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.75)
    )
  ggsave(file.path(out_dir, paste0(prefix, "_full_dotplot.jpg")), width = 7, height = 3.5, units = "in")

  vpal <- viridis(n = 20, option = "D", direction = -1)
  f <- FeaturePlot_scCustom(obj, features = genes, raster = TRUE, colors_use = vpal, order = FALSE, num_columns = 3) & theme(plot.title = element_text(size = 24))
  ggsave(plot = f, file.path(out_dir, paste0(prefix, "_featureplot.jpg")), height = 15.875 * 1.5 * length(genes) / 15, width = (3.26 / 2.84) * 10 * 1.5, units = 'in')

  cpal <- pal(length(unique(obj$cellType)))
  v <- Stacked_VlnPlot(obj, features = genes, colors_use = cpal, x_lab_rotate = 45)
  ggsave(plot = v, file.path(out_dir, paste0(prefix, "_vlnplot.jpg")), height = 1.01 * .75 * 12 * length(genes) / 15 + .75 * .17 / 4.5 * 12, width = .75 * 5 * 3 / 1.88, units = "in")
}

# ==============================================================================
# Advanced Integration & Pipeline Features
# ==============================================================================

#' Combine Integrations and Align Metadata
#'
#' @export
combineIntegration <- function(target, skipStudies = c(), in_dir = ".", out_dir = "sketch_v5_all") {
  if (!dir.exists(file.path(out_dir, "BP_subset"))) dir.create(file.path(out_dir, "BP_subset"), recursive = TRUE, showWarnings = FALSE)
  if (!dir.exists(file.path(out_dir, "BP_subset_object"))) dir.create(file.path(out_dir, "BP_subset_object"), recursive = TRUE, showWarnings = FALSE)

  obj_har <- readRDS(file.path(in_dir, "sketch_v5_all", "allCells_harmony_cellType.rds"))
  obj_mnn <- readRDS(file.path(in_dir, "sketch_v5_all", "allCells_mnn_cellType.rds"))
  obj_scv <- readRDS(file.path(in_dir, "sketch_v5_all", "allCells_scvi_cellType.rds"))

  ii <- (names(obj_har) %in% names(obj_har)[(obj_har == target)]) |
        (names(obj_scv) %in% names(obj_scv)[(obj_scv == target)]) |
        (names(obj_mnn) %in% names(obj_mnn)[(obj_mnn == target)])

  keepNames <- names(obj_har)[ii]
  options(Seurat.object.assay.version = "v5")
  obj_list <- list.files(file.path(in_dir, "BP_object"), recursive = TRUE)

  my_function <- function(var1, var2) {
    if (var2 > 3) return(var1) else return(substr(var1, 1, nchar(var1) - 1 - var2))
  }

  kemp <- str_replace(keepNames, '\\.', '\\-')
  split_list <- strsplit(kemp, "_")
  last_elements_nchar <- sapply(split_list, function(x) nchar(x[length(x)]))
  results_mapply <- mapply(my_function, var1 = kemp, var2 = last_elements_nchar)
  kemp3 <- as.character(results_mapply)

  for (i in seq_along(obj_list)) {
    print(i)
    obj <- readRDS(file.path(in_dir, "BP_object", obj_list[i]))
    studyName <- strsplit(obj_list[i], '/')[[1]][1]
    obj_list[i] <- strsplit(obj_list[i], '/')[[1]][2]
    
    obj$keep <- FALSE
    temp <- str_replace(colnames(obj), '\\.', '\\-')
    split_list <- strsplit(temp, "_")
    last_elements_nchar <- sapply(split_list, function(x) nchar(x[length(x)]))
    results_mapply <- mapply(my_function, var1 = temp, var2 = last_elements_nchar)
    temp3 <- as.character(results_mapply)
    
    obj$keep[temp3 %in% kemp3] <- TRUE

    if (sum(obj$keep) < 60) {
      print('Too few cells')
      next
    }
    
    obj <- subset(obj, subset = keep == TRUE)
    counts.mat <- obj[["RNA"]]$counts
    counts.mat <- convert_matrix_type(counts.mat, type = "uint32_t")
    counts.mat <- as(counts.mat, Class = "dgCMatrix")
    colnames(counts.mat) <- paste0(studyName, '_', colnames(counts.mat))
    
    mat_dir <- file.path(out_dir, "BP_subset", paste0(target, "_it1_", substring(obj_list[i], 1, nchar(obj_list[i]) - 4)))
    write_matrix_dir(mat = counts.mat, dir = mat_dir)
    counts.mat <- open_matrix_dir(dir = mat_dir)
    
    obj <- CreateSeuratObject(counts = counts.mat, meta.data = obj@meta.data)
    saveRDS(obj, file.path(out_dir, "BP_subset_object", paste0(target, "_it1_", obj_list[i])), compress = FALSE)
  }

  options(Seurat.object.assay.version = "v5")
  obj_list <- list.files(file.path(out_dir, "BP_subset_object"), pattern = paste0(target, "_it1_[A-G]*"), recursive = FALSE)
  caf_list <- list()
  
  for (i in seq_along(obj_list)) {
    obj <- readRDS(file.path(out_dir, "BP_subset_object", obj_list[i]))
    caf_list <- append(caf_list, obj)
  }
  
  merged <- merge(caf_list[[1]], caf_list[2:length(caf_list)])
  saveRDS(merged, file.path(out_dir, paste0(target, "_it1.rds")), compress = FALSE)

  return(merged)
}


#' Iterative Filtering and Re-Integration
#'
#' @export
iterativeFilter <- function(obj, name, useSketch = FALSE, skipToScale = FALSE, out_dir = "sketch_v5_all") {
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  
  if (!skipToScale) {
    obj <- NormalizeData(obj)
    saveRDS(obj, file.path(out_dir, paste0(name, "_norm.rds")), compress = FALSE)
    if (length(Layers(obj)) < 4) obj[["RNA"]] <- split(obj[["RNA"]], f = obj$study)
    
    obj <- FindVariableFeatures(obj, verbose = TRUE)
    saveRDS(obj, file.path(out_dir, paste0(name, "_var1.rds")), compress = FALSE)
    
    if (useSketch) {
      obj <- SketchData(obj, ncells = 1000, method = "LeverageScore", sketched.assay = "sketch")
      DefaultAssay(obj) <- "sketch"
      saveRDS(obj, file.path(out_dir, paste0(name, "_sketch.rds")), compress = FALSE)
      obj <- FindVariableFeatures(obj, verbose = TRUE)
      saveRDS(obj, file.path(out_dir, paste0(name, "_var2.rds")), compress = FALSE)
    }
  }
  
  obj <- ScaleData(obj, verbose = FALSE)
  obj <- RunPCA(obj, npcs = 30, verbose = TRUE)
  saveRDS(obj, file.path(out_dir, paste0(name, "_pca.rds")), compress = FALSE)
  
  obj <- IntegrateLayers(obj, method = HarmonyIntegration, orig.reduction = "pca", new.reduction = "harmony", kmeans_init_nstart = 20, kmeans_init_iter_max = 5000, verbose = TRUE)
  saveRDS(obj, file.path(out_dir, paste0(name, "_integration.rds")), compress = FALSE)
  
  obj <- RunUMAP(obj, reduction = "harmony", dims = 1:30, reduction.name = "umap.harmony", return.model = TRUE, min.dist = 0.001)
  saveRDS(obj, file.path(out_dir, paste0(name, "_umap.rds")), compress = FALSE)
  
  options(future.globals.maxSize = 8000 * 1024^2)
  obj <- FindNeighbors(obj, reduction = "harmony", dims = 1:30)
  saveRDS(obj, file.path(out_dir, paste0(name, "_neighbors.rds")), compress = FALSE)
  
  for (i in c(0.06, 0.07, 0.08, 0.09, 0.1, 0.12, 0.14, 0.16, 0.18, 0.2)) {
    g_name <- if (useSketch) 'sketch_snn' else 'RNA_snn'
    obj <- FindClusters(obj, resolution = i, graph.name = g_name)
    DimPlot(obj, raster = TRUE, label = TRUE)
    ggsave(file.path(out_dir, paste0(name, "_clusters_res", i, ".jpg")), width = 5, height = 5, units = "in", limitsize = FALSE)
    saveRDS(obj, file.path(out_dir, paste0(name, "_clusters_res", i, ".rds")), compress = FALSE)
  }
  return(obj)
}


#' Project Sketched Data
#'
#' @export
projectSketchedData <- function(obj, name, sketchedLabels, out_dir = "sketch_v5_all") {
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  
  obj$seurat_clusters <- sketchedLabels
  Idents(obj) <- obj$seurat_clusters
  obj <- ProjectIntegration(object = obj, reduction = "harmony")
  
  options(future.globals.maxSize = 8000 * 1024^2)
  obj <- ProjectData(object = obj, sketched.reduction = "harmony.full", full.reduction = "harmony.full", umap.model = "umap.harmony", dims = 1:30, refdata = list(harmony.cluster.full = "seurat_clusters"))
  
  saveRDS(obj, file.path(out_dir, paste0("unfiltered_counts_obj_findvar_sketch1000_harmony_umap_project_", name, ".rds")), compress = FALSE)
  return(obj)
}


#' Build from Sample Split
#'
#' @export
buildFromSamples <- function(target, iteration, useSketch = FALSE, addIdents = FALSE, skipList = c('placeHolder'), in_dir = "BP_subset_object", out_dir = "sketch_v5_all") {
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  
  obj_list <- list.files(in_dir, pattern = paste0(target, "_it", iteration, "_[A-G]*"), recursive = FALSE)
  caf_list <- list()
  bad_list <- c()
  cafNumber <- 0
  
  for (i in seq_along(obj_list)) {
    print(i)
    isBad <- any(sapply(skipList, function(sl) grepl(sl, obj_list[i])))
    
    if (isBad) {
      bad_list <- c(bad_list, obj_list[i])
    } else {
      cafNumber <- cafNumber + 1
      obj <- readRDS(file.path(in_dir, obj_list[i]))
      obj <- NormalizeData(obj)
      obj <- FindVariableFeatures(obj)
      
      if (useSketch) {
        obj <- SketchData(obj, ncells = 1000, method = "LeverageScore", sketched.assay = "sketch")
        obj <- FindVariableFeatures(obj)
      }
      caf_list[[cafNumber]] <- obj
    }
  }
  
  saveRDS(caf_list, file.path(out_dir, paste0(target, '_', iteration, '_caf_list.rds')), compress = FALSE)
  print(bad_list)
  
  if (addIdents) {
    merged <- merge(caf_list[[1]], caf_list[2:length(caf_list)], add.cell.ids = 1:length(caf_list))
  } else {
    merged <- merge(caf_list[[1]], caf_list[2:length(caf_list)])
  }
  return(merged)
}


