

  getTopMarkers <- function( obj, only.pos = T, assay = "RNA" ) {
    obj <- JoinLayers(obj)
    obj[[assay]]$data <- as(object = obj[[assay]]$data, Class = "dgCMatrix")
    Idents(obj) <- factor(x = Idents(obj), levels = sort(levels(obj), decreasing=F))
    pbmc.markers = FindAllMarkers( obj, only.pos = only.pos, min.pct = -1.10, logfc.threshold = 0.25, max.cells.per.ident = 5000 )
    write.csv( pbmc.markers, '~/pbmc.markers.csv', quote=F )
    return(obj)
  }


  makeDimPlot <- function( obj, raster=T ) {
    DimPlot(obj,group.by='cellType',raster=raster,label=F); ggsave('~/featurePlot.jpg',height=8,widt=9,limitsize=F)
  }

 
 makeDotPlot <- function( obj, features, dotColor='black', w=6, h=6 ) {
    #Idents(obj) = obj$cellType
    Idents(obj) <- factor(x = Idents(obj), levels = sort(levels(obj), decreasing=T))
    lowHex = switch( dotColor,
      'black' ="#DAE2F0",
      'blue'  ="#EEF4FB",
      'green' ="#E2F6EC",
      'orange'="#FAF0EA",
      'pink'  ="#FEEFF9",
      'purple'="#F6EFFC",
      'yellow'="#F9F7D7" )
    highHex = switch( dotColor,
      'black'="#496079",
      'blue'="#9FBCE1",
      'green'="#A5DEC9",
      'orange'="#F5D5C0",
      'pink'="#E5BCDA",
      'purple'="#BCAADD",
      'yellow'="#F7E97D" )

    p = DotPlot(obj, features = features, scale.min = 0, scale.max = 100) +
      NoLegend() + scale_colour_gradient(low=lowHex, high=highHex) +
      theme(axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          axis.text.x = element_text(angle=90, hjust=1, vjust = 0.5),
          axis.text.y = element_text(hjust=0),
          axis.line=element_blank(),
          panel.border = element_rect(colour = "black", fill=NA, size=0.75))
    p
    ggsave('~/featurePlot.jpg',height=h,widt=w,limitsize=F)
    return(p)
  }


  analyzeObj <- function( obj, name, skipStudies = c() ) {
    if( length( Layers(obj) ) > 3 ) {
      obj = JoinLayers(obj)
    }
    temp = table(obj$study); badStudies = names(temp)[ temp<=9 ]; `%notin%` <- Negate(`%in%`)
    cellsToKeep = colnames(obj)[ obj$study %notin% c(badStudies,skipStudies) ]
    obj = subset( obj, cells = cellsToKeep )
    obj@reductions <- list(); DefaultAssay(obj) <- "RNA"; 
    if( length(Layers(obj)) > 1 ) {; obj[["RNA"]]$data <- NULL; }
    obj[["sketch"]] <- NULL; obj$leverage.score <- NULL; obj$umap1 <- NULL; obj$umap2 <- NULL
    counts.mat = obj[['RNA']]$counts; counts.mat = as(object = counts.mat, Class = "dgCMatrix")
    obj <- CreateSeuratObject(counts = counts.mat, meta.data = obj@meta.data)
    write_matrix_dir(mat = obj[["RNA"]]$counts, dir = paste0("sketch_v5_all/",name,"_sketch1000"), overwrite=T)
    counts.mat <- open_matrix_dir(dir = paste0("sketch_v5_all/",name,"_sketch1000"))
    counts.mat <- as(counts.mat, Class = "dgCMatrix")
    obj[["RNA"]]$counts <- counts.mat
    saveRDS(obj, paste0( "sketch_v5_all/",name,"_sketch1000.rds"), compress = F )

    options(future.globals.maxSize = 10 * 1024^3) # 10 GB
    obj <- NormalizeData(obj); obj[["RNA"]] <- split(obj[["RNA"]], f = obj$study)
    obj <- FindVariableFeatures(obj, verbose = T); obj <- ScaleData(obj, verbose = T)
    obj <- RunPCA(obj, npcs = 30, verbose = T)
    obj <- IntegrateLayers(obj, method = HarmonyIntegration, orig.reduction = "pca", new.reduction = "harmony", kmeans_init_nstart=20, kmeans_init_iter_max=5000, verbose = T)
    obj <- RunUMAP(obj, reduction = "harmony", dims = 1:30, reduction.name = "umap.harmony", return.model = T, min.dist = 0.001)
    obj <- FindNeighbors(obj, reduction = "harmony", dims = 1:30)
    for (i in c(0.06, 0.07, 0.09, 0.1, 0.12, 0.14, 0.16, 0.18, 0.2, 0.3, 0.4, 0.5 ) ) {
      obj <- FindClusters(obj, resolution = i)
      DimPlot(obj, raster = T, label=T)
      ggsave(paste0("sketch_v5_all/",name,"_counts_obj_findvar_harmony_umap_res",i,".jpg"),width=5,height=5,units="in",limitsize = F)
      saveRDS(obj, paste0("sketch_v5_all/",name,"_counts_obj_findvar_harmony_umap_res",i,".rds"), compress = F)
    }
    return(obj)
  }



makeMainPanel <- function( obj, name, genes, lowHex='#EEF4FB', highHex='#9FBCE1' ) {
  pal <- colorRampPalette(c("#7C1E6F", "#CF597F", "#F0746E", "#EEB57A", "#FDDE9C", "#C6E174", "#9CCB85", "#7CCBA1", "#3AB185", "#029099", "#12739D", "#055275"))  
 

  Idents(obj) <- factor(obj$cellType, levels = rownames(table(obj$cellType)))
  DimPlot(obj, reduction = "umap.harmony", raster = T, cols = (pal(length(unique(obj$cellType))))) + theme(plot.title=element_blank())
  ggsave(paste0("sketch_v5_all/final_figures/",name,"_legend.jpg"), width = 4, height = 3, units = "in", limitsize = FALSE)
  DimPlot(obj, reduction = "umap.harmony", raster = T, cols = (pal(length(unique(obj$cellType))))) + NoLegend() + theme(plot.title=element_blank())
  ggsave(paste0("sketch_v5_all/final_figures/",name,".jpg"), width = 6, height = 6, units = "in", limitsize = FALSE)

  Idents(obj) <- factor(Idents(obj), levels = rev(levels(obj)))
  d <- DotPlot(obj, features = genes, scale.min = 0, scale.max = 100)
  d + NoLegend() + scale_colour_gradient(low=lowHex, high=highHex) +
    theme(axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          axis.text.x = element_text(angle=90, hjust=1),
          axis.text.y = element_text(hjust=0),
          axis.line=element_blank(),
          panel.border = element_rect(colour = "black", fill=NA, size=0.75))

  ggsave(paste0("sketch_v5_all/final_figures/",name,"_dotplot.jpg"),
         width = 3.5, height = 3.5, units = "in")

}


makeSuppPanel <- function( obj, name, genes, neg, dotColor='black', rasterFeaturePlot = T ) {
  pal <- colorRampPalette(c("#7C1E6F", "#CF597F", "#F0746E", "#EEB57A", "#FDDE9C", "#C6E174", "#9CCB85", "#7CCBA1", "#3AB185", "#029099", "#12739D", "#055275"))
  lowHex = switch( dotColor,
    'black' ="#DAE2F0",
    'blue'  ="#EEF4FB",
    'green' ="#E2F6EC",
    'orange'="#FAF0EA",
    'pink'  ="#FEEFF9",
    'purple'="#F6EFFC",
    'yellow'="#F9F7D7" )
  highHex = switch( dotColor,
    'black'="#496079",
    'blue'="#9FBCE1",
    'green'="#A5DEC9",
    'orange'="#F5D5C0",
    'pink'="#E5BCDA",
    'purple'="#BCAADD",
    'yellow'="#F7E97D" )


  obj$cellType_specific = obj$cellType
  library(viridis); library(scCustomize)
  Idents(obj) <- factor(obj$cellType_specific, levels = rownames(table(obj$cellType_specific)))
  DimPlot(obj, reduction = "umap.harmony", raster = T, cols = (pal(length(unique(obj$cellType_specific))))) + theme(plot.title=element_blank())
  ggsave(paste0("sketch_v5_all/final_figures/",name,"_specific_legend.jpg"), width = 6, height = 6, units = "in", limitsize = FALSE)
  DimPlot(obj, reduction = "umap.harmony", raster = rasterFeaturePlot, cols = (pal(length(unique(obj$cellType_specific))))) + NoLegend() + theme(plot.title=element_blank())
  ggsave(paste0("sketch_v5_all/final_figures/",name,"_specific.jpg"       ), width = 6, height = 6, units = "in", limitsize = FALSE)
  Idents(obj) <- factor(Idents(obj), levels = rev(levels(obj)))
  d <- DotPlot(obj, features = unique(c(genes,neg)), scale.min = 0, scale.max = 100)
  d + NoLegend() + scale_colour_gradient(low=lowHex, high=highHex) +
    theme(axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          axis.text.x = element_text(angle=90, hjust=1, vjust = 0.5),
          axis.text.y = element_text(hjust=0),
          axis.line=element_blank(),
          panel.border = element_rect(colour = "black", fill=NA, size=0.75))

  ggsave(paste0("sketch_v5_all/final_figures/",name,"_specific_dotplot.jpg"),
         width = 7, height = 3.5, units = "in")

  # supplementary feature and vln plots
  pal <- viridis(n = 20, option = "D", direction = +1)
  f <- FeaturePlot_scCustom(obj, features = genes, raster = rasterFeaturePlot, colors_use = pal, order = F, num_columns = 3, max.cutoff = 'q95') & 
                                 theme(plot.title = element_text(size = 24))
  ggsave(plot = f, paste0("sketch_v5_all/final_figures/",name,"_specific_featureplot.jpg"),
         #width = 12*1.50*.5, height = 21*.5, units = "in")
         width = 5*2.7, height = 5*4.8, units = "in")

  Idents(obj) <- factor(x = Idents(obj), levels = sort(levels(obj), decreasing=F))
  pal <- colorRampPalette(c("#7C1E6F", "#CF597F", "#F0746E", "#EEB57A", "#FDDE9C", "#C6E174", "#9CCB85", "#7CCBA1", "#3AB185", "#029099", "#12739D", "#055275"))
  pal <- pal(length(unique(obj$cellType_specific)))

  v <- Stacked_VlnPlot(obj, features = genes, colors_use = pal, x_lab_rotate = 45)
  ggsave(plot = v, paste0("sketch_v5_all/final_figures/",name,"_specific_vlnplot.jpg"), width = 5*1.00, height = 10.5, units = "in")
  #ggsave(plot = v, paste0("sketch_v5_all/final_figures/",name,"_specific_vlnplot.jpg"), width = 2.9*1.8, height = 2.6*1.8, units = "in")
}


# makeBarGraph( obj_cal, 'Lymphoid', obj$cancer, height = 5, width=10 )
makeBarGraph <- function( obj, name, type, excludeCellType = 'NULL', height = 4, width=8 ) {
  library(patchwork)
  pal <- colorRampPalette(c("#7C1E6F", "#CF597F", "#F0746E", "#EEB57A", "#FDDE9C", "#C6E174", "#9CCB85", "#7CCBA1", "#3AB185", "#029099", "#12739D", "#055275"))
  color <- pal(length(unique(obj$cellType[ obj$cellType != excludeCellType])) )
      ggplot(obj@meta.data[obj$cellType!=excludeCellType,],aes(x = type[obj$cellType!=excludeCellType], fill = cellType)) +
      labs(x = "Cancer Type", y = "Proportion (%)", fill = "Subtype") +
      geom_bar(position = "fill") +
      scale_fill_manual(values = color) +
      theme(
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, color = "black", size = 12),
        axis.title.x = element_blank(), axis.text.y = element_blank(),
        axis.ticks.y = element_blank(), axis.title.y = element_text(color = "black"), #element_blank(),
        panel.background = element_blank(),
        plot.background = element_blank(),
        axis.line.x = element_line(color = "black")
      ) +
      scale_x_discrete(expand = c(0, 0)) +
      scale_y_continuous(labels = scales::percent_format(scale = 100),  limits = c(0, 1))

  
  ggsave(paste0("sketch_v5_all/final_figures/",name,"_combined_bargraph.jpg"), width = width, height = height)
}




makeFinalFigures = function( obj, name, genes, dotColor = 'black' ) {
  pal <- colorRampPalette(c("#7C1E6F", "#CF597F", "#F0746E", "#EEB57A", "#FDDE9C", "#C6E174", "#9CCB85", "#7CCBA1", "#3AB185", "#029099", "#12739D", "#055275")) 
  library(scales)
  lowHex = switch( dotColor,
    'black' ="#DAE2F0",
    'blue'  ="#EEF4FB",
    'green' ="#E2F6EC",
    'orange'="#FAF0EA",
    'pink'  ="#FEEFF9",
    'purple'="#F6EFFC",
    'yellow'="#F9F7D7" )
  highHex = switch( dotColor,
    'black'="#496079",
    'blue'="#9FBCE1",
    'green'="#A5DEC9",
    'orange'="#F5D5C0",
    'pink'="#E5BCDA",
    'purple'="#BCAADD",
    'yellow'="#F7E97D" )

  # main figure -- broad categories
  Idents(obj) <- factor(obj$cellType, levels = rownames(table(obj$cellType)))
  #DimPlot(obj, reduction = "umap.harmony", raster = F, cols = (pal(length(unique(obj$cellType))))) + theme(plot.title=element_blank())
  DimPlot(obj, raster = F, cols = (pal(length(unique(obj$cellType))))) + theme(plot.title=element_blank())
  ggsave(paste0("sketch_v5_all/final_figures/",name,"_legend.jpg"), width = 4, height = 3, units = "in", limitsize = FALSE)
  #DimPlot(obj, reduction = "umap.harmony", raster = F, cols = (pal(length(unique(obj$cellType))))) + NoLegend() + theme(plot.title=element_blank())
  DimPlot(obj, raster = F, cols = (pal(length(unique(obj$cellType))))) + NoLegend() + theme(plot.title=element_blank())
  ggsave(paste0("sketch_v5_all/final_figures/",name,".jpg"), width = 6, height = 6, units = "in", limitsize = FALSE)

  Idents(obj) = factor(Idents(obj), levels = sort(unique(Idents(obj))) )
  Idents(obj) <- factor(Idents(obj), levels = rev( sort(unique(Idents(obj))) ) )
  d <- DotPlot(obj, features = genes, scale.min = 0, scale.max = 100)
  d + NoLegend() + scale_colour_gradient(low=lowHex, high=highHex) +  
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_text(angle=90, hjust=1),
        axis.text.y = element_text(hjust=0),
        axis.line=element_blank(),
        #panel.border = element_rect(colour = "black", fill=NA, size=0.75))
        panel.border = element_rect(colour = "black", fill=NA, linewidth=0.75))

  ggsave(paste0("sketch_v5_all/final_figures/",name,"_dotplot.jpg"), width = 3.5, height = 3.5, units = "in")
}


makeSuppFinal <- function( obj, prefix, dotColor, genes, extraGenes ) {

  dotplot_genes <- c(genes, extraGenes)
  pal <- colorRampPalette(c("#7C1E6F", "#CF597F", "#F0746E", "#EEB57A",  "#FDDE9C", "#C6E174", "#9CCB85", "#7CCBA1", "#3AB185", "#029099", "#12739D", "#055275"))

  library(viridis); library(scCustomize)

  Idents(obj) <- factor(obj$cellType, levels = rownames(table(obj$cellType)))

  DimPlot(obj, reduction = "umap.harmony", raster = T, cols = (pal(length(unique(obj$cellType))))) + theme(plot.title=element_blank()) 
  ggsave( paste0( "sketch_v5_all/final_figures/", prefix, "_legend.jpg"), width = 6, height = 6, units = "in", limitsize = F)

  DimPlot(obj, reduction = "umap.harmony", raster = T, cols = (pal(length(unique(obj$cellType))))) + theme(plot.title=element_blank()) + NoLegend()
  ggsave( paste0( "sketch_v5_all/final_figures/", prefix,        ".jpg"), width = 6, height = 6, units = "in", limitsize = F)


  lowHex = switch( dotColor,
    'blue'  ="#EEF4FB",
    'green' ="#E2F6EC",
    'orange'="#FAF0EA",
    'pink'  ="#FEEFF9",
    'purple'="#F6EFFC",
    'yellow'="#F9F7D7" )
  highHex = switch( dotColor,
    'blue'="#9FBCE1",
    'green'="#A5DEC9",
    'orange'="#F5D5C0",
    'pink'="#E5BCDA",
    'purple'="#BCAADD",
    'yellow'="#F7E97D" )
  Idents(obj) <- factor(Idents(obj), levels = rev(levels(obj)))
  d <- DotPlot(obj, features = dotplot_genes, scale.min = 0, scale.max = 100)
  d + NoLegend() + scale_colour_gradient(low=lowHex, high=highHex) +
    theme(axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          axis.text.x = element_text(angle=90, hjust=1, vjust = 0.5),
          axis.text.y = element_text(hjust=0),
          axis.line=element_blank(),
          panel.border = element_rect(colour = "black", fill=NA, size=0.75))
  Idents(obj) <- factor(Idents(obj), levels = rev(levels(obj)))
  ggsave( paste0("sketch_v5_all/final_figures/",prefix,"_full_dotplot.jpg"), width = 7, height = 3.5, units = "in")


  pal <- viridis(n = 20, option = "D", direction = -1)
  f <- FeaturePlot_scCustom(obj, features = genes, raster = TRUE, colors_use = pal, order = F, num_columns = 3) &
    theme(plot.title = element_text(size = 24)) #& xlim(-7, 7) & ylim(-10, 7)

  ggsave(plot = f, paste0("sketch_v5_all/final_figures/",prefix,"_featureplot.jpg"), height = 15.875*1.5*length(genes)/15, width = (3.26/2.84)*10*1.5, units = 'in' )

  #Idents(obj) <- factor(Idents(obj), levels = rev(levels(obj)))
  pal <- colorRampPalette(c("#7C1E6F", "#CF597F", "#F0746E", "#EEB57A",  "#FDDE9C", "#C6E174", "#9CCB85", "#7CCBA1", "#3AB185", "#029099", "#12739D", "#055275"))
  pal <- pal(length(unique(obj$cellType)))
  v <- Stacked_VlnPlot(obj, features = genes, colors_use = pal, x_lab_rotate = 45)
  ggsave(plot = v, paste0("sketch_v5_all/final_figures/",prefix,"_vlnplot.jpg"), height = 1.01*.75*12*length(genes)/15 + .75*.17/4.5*12, width = .75*5*3/1.88, units = "in")

}




combineIntegration <- function( target, skipStudies=c() ) {

  obj_har = readRDS( "sketch_v5_all/allCells_harmony_cellType.rds")
  obj_mnn = readRDS( "sketch_v5_all/allCells_mnn_cellType.rds")
  obj_scv = readRDS( "sketch_v5_all/allCells_scvi_cellType.rds")


  ii = ( names(obj_har) %in% names(obj_har)[ (obj_har == target) ] ) |   # Will be F
       ( names(obj_scv) %in% names(obj_scv)[ (obj_scv == target) ] ) |
       ( names(obj_mnn) %in% names(obj_mnn)[ (obj_mnn == target) ] )

  keepNames <- names(obj_har)[ii]

  options(Seurat.object.assay.version = "v5")
  obj_list <- list.files(paste0("BP_object/"), recursive = T)

  library(stringr)
    my_function <- function(var1, var2) {;
      if( var2 > 3 ) {
        return( var1 )
      } else {
        return( substr( var1, 1, nchar(var1)-1-var2 ) );
      }
    }

  kemp = keepNames
  kemp = str_replace( kemp, '\\.', '\\-' )
  split_list <- strsplit( kemp, "_")
  last_elements_nchar <- sapply(split_list, function(x) nchar( x[length(x)] ) )

  results_mapply <- mapply(my_function, var1 = kemp, var2 = last_elements_nchar )
  kemp3 = as.character(results_mapply)


  start <- 1; end <- length(obj_list)
  for (i in start:end) {; print(i)
    obj <- readRDS(paste0("BP_object/",obj_list[i]))
    studyName = strsplit( obj_list[i], '/' )[[1]][1]
    obj_list[i] = strsplit( obj_list[i], '/' )[[1]][2]
    obj$keep <- F
      temp = colnames(obj)
      temp = str_replace( temp, '\\.', '\\-' )
      split_list <- strsplit( temp, "_")
      last_elements_nchar <- sapply(split_list, function(x) nchar( x[length(x)] ) )
      results_mapply <- mapply(my_function, var1 = temp, var2 = last_elements_nchar )
      temp3 = as.character(results_mapply)
    obj$keep[ temp3 %in% kemp3 ] <- T

    if( sum(obj$keep) < 60 ) {; print('Too few cells'); next; }
    obj <- subset(obj, subset = keep == TRUE)
    # re-perform BP cells 
    counts.mat <- obj[["RNA"]]$counts     ### Implement pasge0's for the target
    counts.mat <- convert_matrix_type(counts.mat, type = "uint32_t" )
    counts.mat <- as(counts.mat, Class = "dgCMatrix")
    colnames(counts.mat) = paste0( studyName, '_', colnames(counts.mat) )
    write_matrix_dir(mat = counts.mat, dir = paste0("BP_subset/",target,"_it1_",substring(obj_list[i],0,nchar(obj_list[i])-4)))
    counts.mat <- open_matrix_dir(dir = paste0("BP_subset/",target,"_it1_",substring(obj_list[i],0,nchar(obj_list[i])-4)))
    obj <- CreateSeuratObject(counts = counts.mat, meta.data = obj@meta.data)
    print(head(colnames(obj)))
    saveRDS(obj, paste0("BP_subset_object/",target,"_it1_",obj_list[i]), compress = F)
  }

  ### merge as a new layer in seurat v5
  options(Seurat.object.assay.version = "v5")
  obj_list <- list.files(paste0("BP_subset_object/"), pattern = paste0(target,"_it1_[A-G]*"), recursive = F)
  caf_list <- list()
  cell_count <- 0
  for (i in 1:length(obj_list)) {
    obj <- readRDS(paste0("BP_subset_object/",obj_list[i]))
    caf_list <- append(caf_list, obj)
    cell_count <- cell_count + dim(obj)[2]
  }
  cell_count
  merged <- merge(caf_list[[1]], caf_list[2:length(caf_list)])
  saveRDS(merged, paste0("sketch_v5_all/",target,"_it1.rds"), compress = F)

  return( merged )

}


########################### Iterative filtering ########################### 
iterativeFilter <- function( obj, name, useSketch=F, skipToScale=F ) {
  if( !skipToScale ) {
    obj <- NormalizeData(obj)
    saveRDS(obj, paste0( "sketch_v5_all/", name, "_norm.rds" ), compress = F )
    if( length( Layers(obj) ) < 4 ) {
      obj[["RNA"]] <- split(obj[["RNA"]], f = obj$study ); 
    }
    obj <- FindVariableFeatures(obj, verbose = T)
    saveRDS(obj, paste0( "sketch_v5_all/", name, "_var1.rds" ), compress = F )
    if( useSketch ) {
      obj <- SketchData(obj, ncells = 1000, method = "LeverageScore", sketched.assay = "sketch")
      DefaultAssay(obj) <- "sketch"
      saveRDS(obj, paste0( "sketch_v5_all/", name, "_sketch.rds" ), compress = F )
      obj <- FindVariableFeatures(obj, verbose = T)
      saveRDS(obj, paste0( "sketch_v5_all/", name, "_var2.rds" ), compress = F ); 
    }
  }
  obj <- ScaleData(obj, verbose = F)
  obj <- RunPCA(obj, npcs = 30, verbose = T)
  saveRDS(obj, paste0( "sketch_v5_all/", name, "_pca.rds" ), compress = F )
  obj <- IntegrateLayers(obj, method = HarmonyIntegration, orig.reduction = "pca", new.reduction="harmony", kmeans_init_nstart=20, kmeans_init_iter_max=5000, verbose = T)
  saveRDS(obj, paste0( "sketch_v5_all/", name, "_integration.rds" ), compress = F )
  obj <- RunUMAP(obj, reduction = "harmony", dims = 1:30, reduction.name = "umap.harmony", return.model = T, min.dist = 0.001)
  saveRDS(obj, paste0( "sketch_v5_all/", name, "_umap.rds" ), compress = F )
  options(future.globals.maxSize = 8000 * 1024^2)
  obj <- FindNeighbors(obj, reduction = "harmony", dims = 1:30)
  saveRDS(obj, paste0( "sketch_v5_all/", name, "_neighbors.rds" ), compress = F )
  for (i in c(0.06, 0.07, 0.08, 0.09, 0.1, 0.12, 0.14, 0.16, 0.18, 0.2)) {
    if( useSketch ) {
      obj <- FindClusters(obj, resolution = i, graph.name = 'sketch_snn')
    } else {
      obj <- FindClusters(obj, resolution = i, graph.name = 'RNA_snn')
    }
    DimPlot(obj, raster = T, label=T)
    ggsave( paste0( "sketch_v5_all/", name, "_clusters_res",i,".jpg" ), width = 5, height = 5, units = "in", limitsize = F)
    saveRDS(obj, paste0( "sketch_v5_all/", name, "_clusters_res",i,".rds" ), compress = F )
  }
  return(obj)
}


################# Project Sketched Data ################# 
projectSketchedData <- function( obj, name, sketchedLabels ) {

  obj$seurat_clusters <- sketchedLabels
  Idents(obj) <- obj$seurat_clusters
  # project
  obj <- ProjectIntegration(object = obj, reduction = "harmony")
  options(future.globals.maxSize = 8000 * 1024^2)
  obj <- ProjectData(object = obj, sketched.reduction = "harmony.full", full.reduction = "harmony.full", umap.model = "umap.harmony", dims = 1:30, refdata = list(harmony.cluster.full = "seurat_clusters"))
  saveRDS(obj, paste0( "sketch_v5_all/unfiltered_counts_obj_findvar_sketch1000_harmony_umap_project_", name, ".rds" ),compress = F)
  return(obj)
}


################### Rebuild, Norm/Var +/- Sketch, Merge ###################
buildFromSamples <- function( target, iteration, useSketch = F, addIdents = F, skipList = c('placeHolder') ) {

  obj_list <- list.files(paste0("BP_subset_object/"), pattern = paste0(target,"_it",iteration,"_[A-G]*" ), recursive = F)
  caf_list <- list()
  bad_list <- c()
  cafNumber = 0
  for (i in 1:length(obj_list)) {; print(i);
      isBad = F
      for( j in 1:length(skipList) ) {;
        if( length( grep(skipList[j], obj_list[i] ) ) ) {
          isBad = T; bad_list = c( bad_list, obj_list[i] )
        }
      }
      if( !isBad ) {
        cafNumber = cafNumber + 1
        obj <- readRDS(paste0("BP_subset_object/",obj_list[i]))
        obj <- NormalizeData(obj)
        obj <- FindVariableFeatures(obj)
        if( useSketch ) {
          obj <- SketchData(obj, ncells = 1000, method = "LeverageScore", sketched.assay = "sketch")
          obj <- FindVariableFeatures(obj)
        }
        caf_list[cafNumber] = obj
      }
  }
  saveRDS( caf_list, paste0( 'sketch_v5_all/', target, '_', iteration, '_caf_list.rds'), compress = F )
  print( bad_list )
  if( addIdents ) {
    merged <- merge( caf_list[[1]], caf_list[2:length(caf_list)], add.cell.ids=1:length(caf_list) )
  } else {
    merged <- merge( caf_list[[1]], caf_list[2:length(caf_list)] )
  }
  return(merged)
}












############################################################################################################
############################################################################################################
####################################### TME Anscillary Functions ###########################################
############################################################################################################
############################################################################################################



runScMetabolism <- function(
  obj, 
  subset = NULL, 
  name = "allCells", 
  input.pathway = c("Glycolysis / Gluconeogenesis", "Oxidative phosphorylation", "Citrate cycle (TCA cycle)")
) {

  # Check if subset is provided
  if( !is.null(subset) ) {
    obj = subset( obj, cellType_fine_abbr %in% subset )
  }

  # Convert to v3 Seurat object
  obj_v3 = obj
  obj_v3[['RNA']] = as(obj[['RNA']], 'Assay')

  # Run scMetabolism
  countexp.Seurat       <- sc.metabolism.Seurat(obj = obj_v3, method = "AUCell", imputation = F, ncores = 2, metabolism.type = "KEGG")
  countexp.Seurat_VISION <- sc.metabolism.Seurat(obj = obj_v3, method = "VISION", imputation = F, ncores = 2, metabolism.type = "KEGG")
  metabolism.matrix <- countexp.Seurat@assays$METABOLISM$score

  # Generate plot using pathways provided in the argument input.pathway
  p <- DotPlot.metabolism(obj = countexp.Seurat, pathway = input.pathway, phenotype = "cellType_fine", norm = "y")

  # Save plot, post-pending the 'name' argument for distinct filenames
  if (!dir.exists("figures")) {; dir.create("figures"); }
  ggsave(paste0("figures/metabolic_fine_", name, ".jpg"), plot = p, width = 12, height = 6)

}




makeScenicPlots <- function( objSce, objAll, name ) {

  obj = subset( objAll, cells = colnames(objAll)[ colnames(objAll) %in% colnames(objSce) ] )

  objSce[[ 'regulon' ]] = objSce[['AUC']]
  obj[[    'regulon' ]] = objSce[[ 'regulon']]
  DefaultAssay(obj) = 'regulon'

  if (!dir.exists('SCENIC')) {
    dir.create('SCENIC')
  }

  library(ggplot2)
  for( i in 1:nrow( obj[['regulon']] ) ) {;
    regulon = substr( rownames(obj[['regulon']])[i], 1, nchar( rownames(obj[['regulon']])[i] ) - 0 )
    print(regulon)
    FeaturePlot( obj, regulon, min.cutoff='q05', max.cutoff='q95', alpha=0.5 ) & scale_color_viridis_c(option = "viridis", direction = 1)
    #scale_color_gradientn(colors = c("lightgrey", "yellow", "green", "darkgreen"))
    ggsave( paste0( 'SCENIC/', name, '_', substr(regulon, 1, nchar(regulon)-3), '.jpg'),  width=8, height=8, units='in', dpi=300)
  }

}




runScenicAnalysis <- function(obj, org = "hgnc", db_dir = "cisTarget_databases", num_cores = 4) {
  
  exprMat <- as.matrix(GetAssayData(obj, assay = "RNA", slot = "counts"))
  
  scenicOptions <- initializeScenic(org = org, dbDir = db_dir, nCores = num_cores)
  
  genesKept <- geneFiltering(exprMat, scenicOptions = scenicOptions, 
                             minCountsPerGene = 3 * .01 * ncol(exprMat), 
                             minSamples = ncol(exprMat) * .01)
  exprMat_filtered <- exprMat[genesKept, ]
  
  runCorrelation(exprMat_filtered, scenicOptions)
  exprMat_filtered_log <- log2(exprMat_filtered + 1)
  runGenie3(exprMat_filtered_log, scenicOptions)
  
  scenicOptions@settings$verbose <- TRUE
  scenicOptions <- runSCENIC_1_coexNetwork2modules(scenicOptions)
  scenicOptions <- runSCENIC_2_createRegulons(scenicOptions)
  scenicOptions <- runSCENIC_3_scoreCells(scenicOptions, exprMat_filtered_log)
  
  auc_matrix <- readRDS(getIntName(scenicOptions, "aucell_regulonAUC"))
  
  # Add the AUC matrix to the Seurat object as a new assay to match downstream plotting
  obj[['AUC']] <- CreateAssayObject(data = getAUC(auc_matrix))
  
  return(obj)
}




runCellChat <- function(
  obj, 
  group.by = "ident", 
  org = "human", 
  db.use = "Secreted Signaling", 
  out_dir = "CellChat_Output", 
  min.cells = 10,
  pathways_to_show = NULL
) {
  
  # Create output directory if it doesn't exist
  if (!dir.exists(out_dir)) {
    dir.create(out_dir, recursive = TRUE)
  }
  
  # If group.by is 'ident', it uses active idents; otherwise it uses the specified metadata column
  if (group.by == "ident") {
    obj$labels <- Idents(obj)
  } else {
    obj$labels <- obj[[group.by]][, 1]
  }
  
  # Drop unused levels to prevent CellChat errors
  obj$labels <- droplevels(as.factor(obj$labels))
  
  cellchat <- createCellChat(object = obj, group.by = "labels")
  cellchat <- setIdent(cellchat, ident.use = "labels")
  groupSize <- as.numeric(table(cellchat@idents))
  
  # Configure database based on the organism
  if (org == "human") {
    CellChatDB <- CellChatDB.human
    ppi_data <- PPI.human
  } else if (org == "mouse") {
    CellChatDB <- CellChatDB.mouse
    ppi_data <- PPI.mouse
  } else {
    stop("Unsupported organism. Please specify 'human' or 'mouse'.")
  }
  
  # Subset the database if a specific category is requested
  if (db.use != "all") {
    CellChatDB.use <- subsetDB(CellChatDB, search = db.use)
  } else {
    CellChatDB.use <- CellChatDB
  }
  
  cellchat@DB <- CellChatDB.use
  
  # Preprocessing and Computing Communication Probability
  cellchat <- subsetData(cellchat)
  cellchat <- identifyOverExpressedGenes(cellchat)
  cellchat <- identifyOverExpressedInteractions(cellchat)
  cellchat <- projectData(cellchat, ppi_data)
  
  # Compute probabilities and filter small groups
  cellchat <- computeCommunProb(cellchat, raw.use = TRUE)
  cellchat <- filterCommunication(cellchat, min.cells = min.cells)
  
  # Compute pathway-level probabilities and aggregate the network
  cellchat <- computeCommunProbPathway(cellchat)
  cellchat <- aggregateNet(cellchat)
  
  # Save the finalized object
  saveRDS(cellchat, file = file.path(out_dir, "cellchat_analyzed.rds"), compress = FALSE)
  
  # Global Visualization Outputs (Number of interactions & Weights)
  jpeg(file.path(out_dir, 'netVisual_circle_all_counts.jpg'), width = 4600, height = 2875, res = 300)
  netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = TRUE, 
                   label.edge = FALSE, title.name = "Number of interactions")
  dev.off()
  
  jpeg(file.path(out_dir, 'netVisual_circle_all_weights.jpg'), width = 4600, height = 2875, res = 300)
  netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = TRUE, 
                   label.edge = FALSE, title.name = "Interaction weights/strength")
  dev.off()
  
  # Pathway-Specific Visualizations (Optional)
  if (!is.null(pathways_to_show)) {
    
    # Calculate network centrality first so role networks can be plotted
    cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")
    
    available_pathways <- cellchat@netP$pathways
    
    for (pathway in pathways_to_show) {
      if (pathway %in% available_pathways) {
        
        # Circle plot for the pathway
        jpeg(file.path(out_dir, paste0('Circle_', pathway, '.jpg')), width = 4600, height = 2875, res = 300)
        netVisual_aggregate(cellchat, signaling = pathway, layout = "circle")
        dev.off()
        
        # Chord diagram for the pathway
        jpeg(file.path(out_dir, paste0('Chord_', pathway, '.jpg')), width = 4600, height = 2875, res = 300)
        netVisual_aggregate(cellchat, signaling = pathway, layout = "chord")
        dev.off()
        
        # Contribution of each Ligand-Receptor pair
        gg <- netAnalysis_contribution(cellchat, signaling = pathway)
        ggsave(filename = file.path(out_dir, paste0("Contribution_", pathway, "_LR.pdf")), 
               plot = gg, width = 5, height = 3, units = 'in', dpi = 300)
        
        # Signaling Role Network (Sender vs Receiver)
        jpeg(file.path(out_dir, paste0('RoleNetwork_', pathway, '.jpg')), width = 4600, height = 2875, res = 300)
        netAnalysis_signalingRole_network(cellchat, signaling = pathway, width = 12, height = 4, font.size = 10)
        dev.off()
        
      } else {
        warning(paste("Pathway", pathway, "was not found as significant in this dataset. Skipping plots for this pathway."))
      }
    }
  }
  
  return(cellchat)
}




runInferCnv <- function(
  obj, 
  group.by = "cellType", 
  ref_group_names = c("Normal"), 
  gene_order_file = "hg38_gencode_v27.txt", 
  out_dir = "infercnv_output", 
  cutoff = 0.1,
  num_threads = 4
) {
  
  library(infercnv)
  library(Seurat)
  
  # Ensure output directory exists
  if (!dir.exists(out_dir)) {
    dir.create(out_dir, recursive = TRUE)
  }
  
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
    warning(paste("The following ref_group_names were not found in the object annotations:", 
                  paste(missing_refs, collapse = ", ")))
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



runCombineIntegrationBetweenSpecies <- function(
  obj_human, 
  obj_mouse, 
  target_celltype = NULL, 
  group.by = "cellType"
) {
  
  library(Seurat)
  library(nichenetr) # Required for convert_mouse_to_human_symbols
  
  # Subset objects if a target cell type is specified
  if (!is.null(target_celltype)) {
    message(paste("Subsetting data to target cell type:", target_celltype))
    
    # Extract metadata columns for subsetting
    human_labels <- obj_human[[group.by]][, 1]
    mouse_labels <- obj_mouse[[group.by]][, 1]
    
    obj_human <- subset(obj_human, cells = colnames(obj_human)[which(human_labels == target_celltype)])
    obj_mouse <- subset(obj_mouse, cells = colnames(obj_mouse)[which(mouse_labels == target_celltype)])
  }
  
  if (ncol(obj_human) == 0 | ncol(obj_mouse) == 0) {
    stop("One or both objects have 0 cells after subsetting")
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
  message(paste("Found", length(common_genes), "common genes."))
  
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
  
  message("Integration complete")
  return(merged_obj)
}



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
  
  # Load required libraries
  library(nichenetr)
  library(Seurat)
  library(tidyverse)
  library(ggplot2)
  
  # Set cell type annotations
  message(paste("Setting cell type annotations to:", cell_type_col))
  Idents(seurat_obj) <- cell_type_col
  
  # Load NicheNet prior models dynamically
  if (!file.exists(ligand_target_matrix_path) | !file.exists(lr_network_path)) {
    stop("NicheNet model files not found.")
  }
  
  message("Loading NicheNet prior models...")
  ligand_target_matrix <- readRDS(ligand_target_matrix_path)
  lr_network <- readRDS(lr_network_path)
  
  # Calculate expressed genes for receivers
  message(paste("Calculating expressed genes for receiver(s):", paste(receivers, collapse = ", ")))
  expressed_genes_receiver <- get_expressed_genes(receivers, seurat_obj, pct = 0.10)
  background_expressed_genes <- expressed_genes_receiver %>% .[. %in% rownames(ligand_target_matrix)]
  
  # Calculate expressed genes for senders
  message(paste("Calculating expressed genes for sender(s):", paste(senders, collapse = ", ")))
  expressed_genes_sender <- get_expressed_genes(senders, seurat_obj, pct = 0.10)
  
  # Filter potential ligands and receptors
  message("Filtering potential ligands and receptors...")
  ligands <- lr_network %>% pull(from) %>% unique()
  receptors <- lr_network %>% pull(to) %>% unique()
  
  expressed_ligands <- intersect(ligands, expressed_genes_sender)
  expressed_receptors <- intersect(receptors, expressed_genes_receiver)
  
  potential_ligands <- lr_network %>% 
    filter(from %in% expressed_ligands & to %in% expressed_receptors) %>% 
    pull(from) %>% unique()
  
  # Filter out blacklisted genes
  if (!is.null(blacklist)) {
    message("Removing blacklisted ambient RNA genes...")
    potential_ligands <- setdiff(potential_ligands, blacklist)
  }
  
  # Predict ligand activities
  message("Predicting ligand activities based on provided DEGs...")
  ligand_activities <- predict_ligand_activities(geneset = degs, 
                                                 background_expressed_genes = background_expressed_genes, 
                                                 ligand_target_matrix = ligand_target_matrix, 
                                                 potential_ligands = potential_ligands)
  
  ligand_activities <- ligand_activities %>% arrange(-pearson) %>% mutate(rank = rank(desc(pearson)))
  top_ligands <- ligand_activities %>% top_n(20, pearson) %>% pull(test_ligand) %>% unique()
  
  # Infer specific ligand-target links
  message("Inferring specific ligand-target links...")
  active_ligand_target_links_df <- top_ligands %>% 
    lapply(get_weighted_ligand_target_links, 
           geneset = degs, 
           ligand_target_matrix = ligand_target_matrix,
           n = 200) %>% 
    bind_rows() %>% drop_na()
  
  active_ligand_target_links <- prepare_ligand_target_visualization(
    ligand_target_df = active_ligand_target_links_df, 
    ligand_target_matrix = ligand_target_matrix, 
    cutoff = 0.33)
  
  # Generate heatmap visualization
  message("Generating heatmap visualization...")
  order_ligands <- intersect(top_ligands, colnames(active_ligand_target_links)) %>% rev()
  order_targets <- active_ligand_target_links_df$target %>% unique() %>% intersect(rownames(active_ligand_target_links))
  vis_ligand_target <- active_ligand_target_links[order_targets, order_ligands] %>% t()
  
  p_ligand_target_network <- vis_ligand_target %>% 
    make_heatmap_ggplot(y_name = "Prioritized Sender Ligands", 
                        x_name = "Predicted Target Genes", 
                        color = "purple", legend_position = "top", 
                        x_axis_position = "top", legend_title = "Regulatory potential") +
    theme(axis.text.x = element_text(angle = 45, hjust = 0, size = 6))
  
  # Save heatmap
  message("Saving heatmap to output directory...")
  
  safe_senders <- paste(gsub("[^A-Za-z0-9]", "_", senders), collapse = "_")
  safe_receivers <- paste(gsub("[^A-Za-z0-9]", "_", receivers), collapse = "_")
  file_name <- paste0(safe_senders, "__", safe_receivers, "_.jpg")
  
  if (!dir.exists(out_dir)) {
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  }
  
  file_path <- file.path(out_dir, file_name)
  
  ggsave(filename = file_path, plot = p_ligand_target_network, width = 14, height = 8, dpi = 300)
  message(paste("Saved plot successfully to:", file_path))
  
  return(list(top_ligands = top_ligands, ligand_activities = ligand_activities, heatmap = p_ligand_target_network))
}



