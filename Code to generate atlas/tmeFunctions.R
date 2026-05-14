

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





