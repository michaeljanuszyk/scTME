### HUMAN META-ANALYSIS CODE
### Written by John Lu and Michael Januszyk
### Date: Feb 28, 2024

### Quality control, pre-processing, and merging data ----
library(Seurat)
library(ggplot2)
library(dplyr)
options(Seurat.object.assay.version = 'v3')

# QC per sample
obj_list <- list.files(paste0("readyForSeurat/"),"*object.rds", recursive = TRUE)
obj_list <- substring(obj_list,0,nchar(obj_list)-11)
write.csv(obj_list, "obj_list_01_preQC.csv")
for (obj in obj_list) {
  x <- readRDS(paste0("readyForSeurat/",obj,".object.rds"))
  x <- subset(x, subset = nFeature_RNA > 200)
  x[["percent.mt"]] <- PercentageFeatureSet(x, pattern = "^mt-")
  VlnPlot(x, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
  x <- subset(x, subset = nFeature_RNA > 200 & nFeature_RNA < quantile(x$nFeature_RNA, 0.95) & percent.mt < 20)
  VlnPlot(x, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
  dir.create(paste0("processed/",strsplit(obj,"/")[[1]][1]),suppressWarnings(FALSE))
  saveRDS(x, paste0("processed/",obj, "_01_qc.rds"), compress = FALSE)
}

# format individual post-QC objects (streamline metadata, rename gene names)
obj_list <- as.vector(read.csv("obj_list_01_preQC.csv", colClasses=c("NULL",NA))$x)
for (i in 1:length(obj_list)) {
  obj <- obj_list[i]
  x <- readRDS(paste0("processed/",obj,"_01_qc.rds"))
  x <- DietSeurat(x)
  x@meta.data <- x@meta.data %>%
    select(any_of(c("orig.ident", "nCount_RNA","nFeature_RNA","study",
                    "organ","cancer","model","sorting","site","percent.mt")))
  x <- RenameCells(x, new.names = paste(x$orig.ident, Cells(x), sep = "_"))
  # deal with duplicate gene names and non-standard gene names
  counts.mat <- GetAssayData( x, slot = 'counts' )
  geneNames <- rownames(counts.mat)
  geneNames <- sub("\\.ps", "-ps", geneNames, ignore.case = T) # replace .ps for -ps
  geneNames <- sub("-ENSMUSG.*$", "", geneNames) # replace eliminate ENSMUSG suffix
  # download human gene names, symbols for GRCh38.p14
  synToName <- read.csv("GRCh38.p14_geneNames_geneSynonyms.csv", header = T)
  for (i in 1:length(geneNames)) {
    name <- geneNames[i]
    if (name %in% synToName$Gene.name) {
      next
    } else if (name %in% synToName$Gene.Synonym){
      geneNames[i] <- synToName$Gene.name[which(synToName$Gene.Synonym == name)]
    } 
  }
  rownames(counts.mat) <- geneNames
  # counts mat now has duplicated row names. sum dup rows
  geneNamesDup <- unique(geneNames[which(duplicated(geneNames))])
  if (length(geneNamesDup) > 0) {
    counts.mat <- convert_matrix_type( counts.mat, type = "uint32_t" )
    dup.counts.mat <- as(matrix(0,nrow = length(geneNamesDup), ncol = ncol(counts.mat)), Class = "dgCMatrix")
    colnames(dup.counts.mat) <- colnames(counts.mat)
    rownames(dup.counts.mat) <- geneNamesDup
    for (i in 1:length(geneNamesDup)) {
      gene <- geneNamesDup[i]
      rows <- which(rownames(counts.mat) == gene)
      temp <- as(counts.mat[rows,], Class = "dgCMatrix")
      temp <- as(t(colSums(temp)), Class = "dgCMatrix")
      dup.counts.mat[i,] <- temp
    }
    duprows <- (which(duplicated(geneNames) | duplicated(geneNames, fromLast = T)))
    counts.mat <- rbind(
      as(counts.mat[-duprows,], Class = "dgCMatrix"),
      dup.counts.mat
    )
  }
  x <- CreateSeuratObject(counts = counts.mat, meta.data = x@meta.data)
  saveRDS(x, paste0("processed/",obj, "_01b_correctGeneName.rds"), compress = FALSE)
}

### generate list of seurat objects and merge
obj_list <- as.vector(read.csv("obj_list_01_preQC.csv", colClasses=c("NULL",NA))$x)
caf.list <- list()
for (obj in obj_list) {
  x <- readRDS(paste0("processed/",obj,"_01b_correctGeneName.rds"))
  caf.list <- append(caf.list, x)
}
obj <- merge(caf.list[[1]],caf.list[2:length(caf.list)])
saveRDS(obj, "sketch_v5_all/merged_newNames.rds", compress = F)
any(GetAssayData( obj, slot = 'counts' )%%1!=0) # confirm no non-integer values

### remove genes not included in standard nomenclature
counts.mat <- GetAssayData( obj, slot = 'counts' )
counts.mat <- convert_matrix_type( counts.mat, type = "uint32_t" )
# remove genes in less than 0.1% of cells
nCount_Feature <- rowSums( counts.mat > 0)
counts.mat <- counts.mat[nCount_Feature >=0.001*ncol(obj), ]
synToName <- read.csv("GRCh38.p14_geneNames_geneSynonyms.csv", header = T)
counts.mat <- counts.mat[rownames(counts.mat) %in% synToName$Gene.name,]
counts.mat <- as(counts.mat, Class = "dgCMatrix")

### BPCell-ize object
library(BPCells)
options(Seurat.object.assay.version = "v5")

obj <- CreateSeuratObject(counts = counts.mat, meta.data = obj@meta.data)
write_matrix_dir(mat = obj[["RNA"]]$counts, dir = "sketch_v5_all/unfiltered_counts")
counts.mat <- open_matrix_dir(dir = "sketch_v5_all/unfiltered_counts")
obj[["RNA"]]$counts <- counts.mat
saveRDS(obj, "sketch_v5_all/unfiltered_counts_obj.rds", compress = F)

### Seurat workflow with sketching 1000 cells per study
obj <- NormalizeData(obj)
obj[["RNA"]] <- split(obj[["RNA"]], f = obj$study)
obj <- FindVariableFeatures(obj, verbose = FALSE)
obj <- SketchData(obj, ncells = 1000, method = "LeverageScore", sketched.assay = "sketch")
DefaultAssay(obj) <- "sketch"
obj <- FindVariableFeatures(obj, verbose = F)
obj <- ScaleData(obj, verbose = F)
obj <- RunPCA(obj, npcs = 30, verbose = F)
saveRDS(obj, "sketch_v5_all/unfiltered_counts_obj_findvar_sketch1000_pca.rds", compress = F)

### Integrate using scVI----
library(Seurat); library(ggplot2); library( dplyr )
library(BPCells); library(Matrix); library(SeuratWrappers)
library(scales)
options(Seurat.object.assay.version = "v5")

Sys.setenv( RETICULATE_PYTHON = "/home/johnlu/.conda/envs/scvi-env/bin/python" )
library(reticulate)

obj <- readRDS("n0c/humanMerged_norm_findvar_sketch1000_pca.rds")
obj <- IntegrateLayers(
  object = obj, method = scVIIntegration,
  orig.reduction = "pca", new.reduction = "integrated.scvi", conda_env = '~/.conda/envs/scvi-env/', verbose = T
)
obj <- RunUMAP(obj, reduction = "integrated.scvi", dims = 1:30, reduction.name = "umap.scvi", return.model = T, min.dist = 0.001)
saveRDS(obj, "n0d/humanMerged_norm_findvar_sketch1000_scvi_umap.rds",compress = F)

obj <- FindNeighbors(obj, reduction = "integrated.scvi", dims = 1:30)
for (i in seq(0.03, 0.1, 0.01)) {  #choose 0.09 
  DefaultAssay(obj) = 'RNA'
  obj <- FindClusters(obj, resolution = i)
  saveRDS(obj, paste0("n0d/humanMerged_norm_findvar_sketch1000_scvi_umap_res",i,".rds"), compress = FALSE)
}

# find markers
obj <- readRDS("n0d/humanMerged_norm_findvar_sketch1000_scvi_umap_res0.09.rds")
DefaultAssay(obj) <- "sketch"
obj <- JoinLayers(obj)
markers <- FindAllMarkers(obj, max.cells.per.ident = 1000, only.pos = T)
write.csv(markers, "n0d/humanMerged_norm_findvar_sketch1000_scvi_umap_res0.09_markers.csv")
markers$rank <- markers$pct.1 / (markers$pct.2 + 0.001) * markers$avg_log2FC
topmarkers <- markers %>%
  group_by(cluster) %>%
  slice_max(n = 100, order_by = rank)
write.csv(topmarkers, "n0d/humanMerged_norm_findvar_sketch1000_scvi_umap_res0.09_topmarkers.csv")

# project to full dataset
obj <- readRDS("n0d/humanMerged_norm_findvar_sketch1000_scvi_umap_res0.09.rds")
obj <- ProjectIntegration(object = obj, reduction = "integrated.scvi")
options(future.globals.maxSize = 8000 * 1024^2)
obj <- ProjectData(object = obj, sketched.reduction = "integrated.scvi.full", full.reduction = "integrated.scvi.full", umap.model = "umap.scvi", dims = 1:30, refdata = list(scvi.cluster.full = "seurat_clusters"))
saveRDS(obj, "n0d/humanMerged_norm_findvar_sketch1000_scvi_umap_res0.09_project.rds",compress = F)
saveRDS(obj$scvi.cluster.full, "n0d/humanMerged_norm_findvar_sketch1000_scvi_umap_res0.09_clusterList.rds")

# join layers
DefaultAssay(obj) <- "sketch"
obj <- JoinLayers(obj)
DefaultAssay(obj) <- "RNA"
obj <- JoinLayers(obj)
saveRDS(obj, "n0d/humanMerged_norm_findvar_sketch1000_scvi_umap_res0.09_project_joined.rds", compress = F)


### Integrate using Harmony----
obj <- readRDS("n0c/humanMerged_norm_findvar_sketch1000_pca.rds")
obj <- IntegrateLayers(obj, method = HarmonyIntegration, orig.reduction = "pca", new.reduction = "harmony", kmeans_init_nstart=20, kmeans_init_iter_max=5000, verbose = T)
obj <- RunUMAP(obj, reduction = "harmony", dims = 1:30, reduction.name = "umap.harmony", return.model = T, min.dist = 0.001)

# clustering
obj <- FindNeighbors(obj, reduction = "harmony", dims = 1:30)
for (i in seq(0.03, 0.1, 0.01) ) { # choose 0.06
  obj <- FindClusters(obj, resolution = i)
  saveRDS(obj, paste0("n0d/humanMerged_norm_findvar_sketch1000_harmony_umap_res",i,".rds"), compress = FALSE)
}

# find clusters
obj <- readRDS("n0d/humanMerged_norm_findvar_sketch1000_harmony_umap_res0.06.rds")
DefaultAssay(obj) <- "sketch"
obj <- JoinLayers(obj)
markers <- FindAllMarkers(obj, max.cells.per.ident = 1000, only.pos = T)
write.csv(markers, "n0d/humanMerged_norm_findvar_sketch1000_harmony_umap_res0.06_markers.csv")
markers$rank <- markers$pct.1 / (markers$pct.2 + 0.001) * markers$avg_log2FC
topmarkers <- markers %>%
  group_by(cluster) %>%
  slice_max(n = 100, order_by = rank)
write.csv(topmarkers, "n0d/humanMerged_norm_findvar_sketch1000_harmony_umap_res0.06_topmarkers.csv")

# project to full dataset
obj <- readRDS("n0d/humanMerged_norm_findvar_sketch1000_harmony_umap_res0.06.rds")
obj <- ProjectIntegration(object = obj, reduction = "harmony")
options(future.globals.maxSize = 8000 * 1024^2)
obj <- ProjectData(object = obj, sketched.reduction = "harmony.full", full.reduction = "harmony.full", umap.model = "umap.harmony", dims = 1:30, refdata = list(harmony.cluster.full = "seurat_clusters"))
saveRDS(obj, "n0d/humanMerged_norm_findvar_sketch1000_harmony_umap_res0.06_project.rds",compress = F)
saveRDS(obj$harmony.cluster.full, "n0d/humanMerged_norm_findvar_sketch1000_harmony_umap_res0.06_clusterList.rds")

# join layers
DefaultAssay(obj) <- "sketch"
obj <- JoinLayers(obj)
DefaultAssay(obj) <- "RNA"
obj <- JoinLayers(obj)
saveRDS(obj, "n0d/humanMerged_norm_findvar_sketch1000_harmony_umap_res0.06_project_joined.rds", compress = F)

### Integrate using MNN----
obj <- readRDS("n0c/humanMerged_norm_findvar_sketch1000_pca.rds")
obj <- IntegrateLayers(
  object = obj, method = FastMNNIntegration,
  orig.reduction = "pca", new.reduction = "integrated.mnn", verbose = T
)
obj <- RunUMAP(obj, reduction = "integrated.mnn", dims = 1:30, reduction.name = "umap.mnn", return.model = T, min.dist = 0.001)
saveRDS(obj, "n0d/humanMerged_norm_findvar_sketch1000_mnn_umap.rds",compress = F)

# find clusters
obj <- FindNeighbors(obj, reduction = "integrated.mnn", dims = 1:30)
DefaultAssay(obj) <- "RNA"
for (i in seq(0.03, 0.1, 0.01)) { # choose 0.08
  obj <- FindClusters(obj, resolution = i)
  saveRDS(obj, paste0("n0d/humanMerged_norm_findvar_sketch1000_mnn_umap_res",i,".rds"), compress = FALSE)
}


# find markers
obj <- readRDS("n0d/humanMerged_norm_findvar_sketch1000_mnn_umap_res0.08.rds")
DefaultAssay(obj) <- "sketch"
obj <- JoinLayers(obj)
markers <- FindAllMarkers(obj, max.cells.per.ident = 1000, only.pos = T)
write.csv(markers, "n0d/humanMerged_norm_findvar_sketch1000_mnn_umap_res0.08_markers.csv")
markers$rank <- markers$pct.1 / (markers$pct.2 + 0.001) * markers$avg_log2FC
topmarkers <- markers %>%
  group_by(cluster) %>%
  slice_max(n = 100, order_by = rank)
write.csv(topmarkers, "n0d/humanMerged_norm_findvar_sketch1000_mnn_umap_res0.08_topmarkers.csv")


# project
obj <- readRDS("n0d/humanMerged_norm_findvar_sketch1000_mnn_umap_res0.08.rds")
obj <- ProjectIntegration(object = obj, reduction = "integrated.mnn")
options(future.globals.maxSize = 8000 * 1024^2)
obj <- ProjectData(object = obj, sketched.reduction = "integrated.mnn.full", full.reduction = "integrated.mnn.full", umap.model = "umap.mnn", dims = 1:30, refdata = list(mnn.cluster.full = "seurat_clusters"))
saveRDS(obj, "n0d/humanMerged_norm_findvar_sketch1000_mnn_umap_res0.08_project.rds",compress = F)
saveRDS(obj$mnn.cluster.full, "n0d/humanMerged_norm_findvar_sketch1000_mnn_umap_res0.08_clusterList.rds")

# join layers
DefaultAssay(obj) <- "sketch"
obj <- JoinLayers(obj)
DefaultAssay(obj) <- "RNA"
obj <- JoinLayers(obj)
saveRDS(obj, "n0d/humanMerged_norm_findvar_sketch1000_mnn_umap_res0.08_project_joined.rds", compress = F)


### Fibroblast Iteration 1 ----
obj     = readRDS( "/oak/stanford/groups/longaker/CAF_META/newHuman/n0d/humanMerged_norm_findvar_sketch1000_harmony_umap_res0.06_project_joined.rds" )
obj_har = readRDS( "/oak/stanford/groups/longaker/CAF_META/newHuman/n0d/humanMerged_norm_findvar_sketch1000_harmony_umap_res0.06_clusterList.rds" )
obj_scv = readRDS( "/oak/stanford/groups/longaker/CAF_META/newHuman/n0d/humanMerged_norm_findvar_sketch1000_scvi_umap_res0.09_clusterList.rds"    )
obj_mnn = readRDS( "/oak/stanford/groups/longaker/CAF_META/newHuman/n0d/humanMerged_norm_findvar_sketch1000_mnn_umap_res0.08_clusterList.rds"    )

# Fibroblasts
ii = ( names(obj_har) %in% names(obj_har)[ (obj_har == 3) ] ) |
     ( names(obj_scv) %in% names(obj_scv)[ (obj_scv == 3) | (obj_scv ==19) ] ) |
     ( names(obj_mnn) %in% names(obj_mnn)[ (obj_mnn == 1) ] )

DefaultAssay(obj) = 'RNA'
obj$keep = FALSE
obj$keep[ ii ] = TRUE
obj <- subset(obj, subset = keep == TRUE )

# re-generate seurat object and remove extraneous data
options(Seurat.object.assay.version = "v5")
obj <- CreateSeuratObject(counts = obj[["RNA"]]$counts, meta.data = obj@meta.data)
saveRDS( obj, 'o0c/union_fibro.rds', compress = F )

# re-perform BPcells
counts.mat = obj[['RNA']]$counts
counts.mat = convert_matrix_type( counts.mat, type="uint32_t" )
counts.mat = as( counts.mat, Class='dgCMatrix' )
write_matrix_dir(mat = counts.mat, dir = "o0c/union_fibro_counts", overwrite = T )
counts.mat <- open_matrix_dir(     dir = "o0c/union_fibro_counts" )
obj = CreateSeuratObject( counts.mat, meta.data = obj@meta.data )
saveRDS( obj, 'o0c/union_fibro_counts.rds', compress = F )

# run downstream analysis
obj <- NormalizeData(obj)
obj[["RNA"]] <- split(obj[["RNA"]], f = obj$study)
obj <- FindVariableFeatures(obj, verbose = T)
obj <- SketchData(object = obj, ncells = 1000, method = "LeverageScore", sketched.assay = "sketch")
DefaultAssay(obj) <- "sketch"
obj <- FindVariableFeatures(obj, verbose = T)
obj <- ScaleData(obj, verbose = F)
obj <- RunPCA(obj, npcs = 30, verbose = F)
saveRDS(obj, "o0c/fibro_it1_var_sketch_pca.rds", compress = F)

obj <- IntegrateLayers(obj, method = HarmonyIntegration, orig.reduction = "pca", new.reduction = "harmony", kmeans_init_nstart=20, kmeans_init_iter_max=5000, verbose = T)
obj <- RunUMAP(obj, reduction = "harmony", dims = 1:30, reduction.name = "umap.harmony", return.model = T, min.dist = 0.001)
saveRDS(obj, "o0c/fibro_it1_var_sketch_pca_harmony_umap.rds", compress = F)

obj <- FindNeighbors(obj, reduction = "harmony", dims = 1:30)
for (i in c(seq(0.01, 0.1, 0.01), 0.025)) { # 0.025
  obj <- FindClusters(obj, resolution = i)
  saveRDS(obj, paste0("o0c/fibro_it1_var_sketch_pca_harmony_umap_res",i,".rds"), compress = FALSE)
}

# find markers
obj <- readRDS("o0c/fibro_it1_var_sketch_pca_harmony_umap_res0.025.rds")
DefaultAssay(obj) <- "sketch"
obj <- JoinLayers(obj)
markers <- FindAllMarkers(obj, max.cells.per.ident = 1000, only.pos = T)
write.csv(markers, "o0c/fibro_it1_var_sketch_pca_harmony_umap_res0.025_markers.csv")
markers$rank <- markers$pct.1 / (markers$pct.2 + 0.001) * markers$avg_log2FC
topmarkers <- markers %>%
  group_by(cluster) %>%
  slice_max(n = 100, order_by = rank)
write.csv(topmarkers, "o0c/fibro_it1_var_sketch_pca_harmony_umap_res0.025_topmarkers.csv")

# project to full dataset
obj <- readRDS("o0c/fibro_it1_var_sketch_pca_harmony_umap_res0.025.rds") # needs to be split (non-joined)
obj <- ProjectIntegration(object = obj, reduction = "harmony")
options(future.globals.maxSize = 8000 * 1024^2)
obj <- ProjectData(object = obj, sketched.reduction = "harmony.full", full.reduction = "harmony.full", umap.model = "umap.harmony", dims = 1:30, refdata = list(harmony.cluster.full = "seurat_clusters"))
saveRDS(obj, "o0c/fibro_it1_var_sketch_pca_harmony_umap_res0.025_project.rds",compress = F)
saveRDS(obj$harmony.cluster.full, "o0c/fibro_it1_var_sketch_pca_harmony_umap_res0.025_project_clusterList.rds")

# join layers
DefaultAssay(obj) <- "sketch"
obj <- JoinLayers(obj)
DefaultAssay(obj) <- "RNA"
obj <- JoinLayers(obj)
saveRDS(obj, "o0c/fibro_it1_var_sketch_pca_harmony_umap_res0.025_project_joined.rds", compress = F)


### Fibroblast Iteration 2 ----

obj <- readRDS( "o0c/fibro_it1_var_sketch_pca_harmony_umap_res0.025_project_joined.rds" )

# Fibroblast
DefaultAssay(obj) = 'RNA'
obj$keep = FALSE
obj$keep[ obj$seurat_clusters %in% c(0,6) ] = TRUE
obj <- subset(obj, subset = keep == TRUE )

# re-generate seurat object and remove extraneous data
options(Seurat.object.assay.version = "v5")
obj <- CreateSeuratObject(counts = obj[["RNA"]]$counts, meta.data = obj@meta.data)
saveRDS( obj, 'n0f/it2_union_fibro.rds', compress = F )

# re-perform BPcells
counts.mat = obj[['RNA']]$counts
counts.mat = convert_matrix_type( counts.mat, type="uint32_t" )
counts.mat = as( counts.mat, Class='dgCMatrix' )
write_matrix_dir(mat = counts.mat, dir = "n0f/it2_union_fibro_counts", overwrite = T )
counts.mat <- open_matrix_dir(     dir = "n0f/it2_union_fibro_counts" )
obj = CreateSeuratObject( counts.mat, meta.data = obj@meta.data )
saveRDS( obj, 'n0f/it2_union_fibro_counts.rds', compress = F )

# run downstream analysis
obj <- NormalizeData(obj)
obj[["RNA"]] <- split(obj[["RNA"]], f = obj$study)
obj <- FindVariableFeatures(obj, verbose = T )
obj <- ScaleData(obj, verbose = T)
obj <- RunPCA(obj, npcs = 30, verbose = T)
obj <- IntegrateLayers(obj, method = HarmonyIntegration, orig.reduction = "pca", new.reduction = "harmony" )
obj <- RunUMAP(obj, reduction = "harmony", dims = 1:30, reduction.name = "umap.harmony", return.model = T, min.dist = 0.001)
saveRDS( obj, 'n0f/it2_union_fibro_harmony_umap.rds', compress = F ) 

# find clusters
obj <- FindNeighbors(obj, reduction = "harmony", dims = 1:30)
for( res in seq( 0.01, 0.1, 0.01 ) ) { # choose 0.04
  obj <- FindClusters(obj, resolution = res )
  saveRDS( obj, paste0( 'n0f/it2_union_fibro_harmony_clusters_', res, '.rds' ), compress = F ) 
}

# join layers and find clusters
obj <- readRDS("n0f/it2_union_fibro_harmony_clusters_0.04.rds")
DefaultAssay(obj) <- "RNA"
obj <- JoinLayers(obj)
saveRDS(obj, "n0f/it2_union_fibro_harmony_clusters_0.04_joined.rds", compress = F)
obj[["RNA"]]$data <- as(object = obj[["RNA"]]$data, Class = "dgCMatrix")
markers <- FindAllMarkers(obj, max.cells.per.ident = 1000, only.pos = T)
write.csv(markers, "n0f/it2_union_fibro_harmony_clusters_0.04_markers.csv")
markers$rank <- markers$pct.1 / (markers$pct.2 + 0.001) * markers$avg_log2FC
topmarkers <- markers %>%
  group_by(cluster) %>%
  slice_max(n = 100, order_by = rank)
write.csv(topmarkers, "n0f/it2_union_fibro_harmony_clusters_0.04_topmarkers.csv")

### Fibroblast Iteration 3 ----

obj <- readRDS( "n0f/it2_union_fibro_harmony_clusters_0.04_joined.rds" )

# Fibroblast
DefaultAssay(obj) = 'RNA'
obj$keep = FALSE
obj$keep[ obj$seurat_clusters %in% c(0, 2, 3, 4) ] = TRUE
obj <- subset(obj, subset = keep == TRUE )

# re-generate seurat object and remove extraneous data
options(Seurat.object.assay.version = "v5")
obj <- CreateSeuratObject(counts = obj[["RNA"]]$counts, meta.data = obj@meta.data)
saveRDS( obj, 'n0e/it3_union_fibro.rds', compress = F )

# re-perform BPcells
counts.mat = obj[['RNA']]$counts
counts.mat = convert_matrix_type( counts.mat, type="uint32_t" )
counts.mat = as( counts.mat, Class='dgCMatrix' )
write_matrix_dir(mat = counts.mat, dir = "n0e/it3_union_fibro_counts", overwrite = T )
counts.mat <- open_matrix_dir(     dir = "n0e/it3_union_fibro_counts" )
obj = CreateSeuratObject( counts.mat, meta.data = obj@meta.data )
saveRDS( obj, 'n0e/it3_union_fibro_counts.rds', compress = F )

# run downstream analysis
obj <- NormalizeData(obj)
obj[["RNA"]] <- split(obj[["RNA"]], f = obj$study)
obj <- FindVariableFeatures(obj, verbose = T )
obj <- ScaleData(obj, verbose = T)
obj <- RunPCA(obj, npcs = 30, verbose = T)
obj <- IntegrateLayers(obj, method = HarmonyIntegration, orig.reduction = "pca", new.reduction = "harmony" )
obj <- RunUMAP(obj, reduction = "harmony", dims = 1:30, reduction.name = "umap.harmony", return.model = T, min.dist = 0.001)
saveRDS( obj, 'n0e/it3_union_fibro_harmony_umap.rds', compress = F ) 

# find clusters
for( res in seq(0.03, 0.25, 0.01) ) {
  obj <- FindClusters(obj, resolution = res )
  saveRDS( obj, paste0( 'n0e/union_fibro3_harmony_clusters_', res, '.rds' ), compress = F ) 
}

# join layers and find clusters
obj <- readRDS("n0e/union_fibro3_harmony_clusters_0.07.rds")
DefaultAssay(obj) <- "RNA"
obj <- JoinLayers(obj)
obj[["RNA"]]$data <- as(object = obj[["RNA"]]$data, Class = "dgCMatrix")
saveRDS(obj, "n0e/union_fibro3_harmony_clusters_0.07_joined.rds", compress = F)
markers <- FindAllMarkers(obj, max.cells.per.ident = 1000, only.pos = T)
write.csv(markers, "n0e/union_fibro3_harmony_clusters_0.07_markers.csv")
markers$rank <- markers$pct.1 / (markers$pct.2 + 0.001) * markers$avg_log2FC
topmarkers <- markers %>%
  group_by(cluster) %>%
  slice_max(n = 100, order_by = rank)
write.csv(topmarkers, "n0e/union_fibro3_harmony_clusters_0.07_topmarkers.csv")

### Fibroblast Iteration 4 ----
obj <- readRDS("n0e/union_fibro3_harmony_clusters_0.07_joined.rds")

# Fibroblasts
DefaultAssay(obj) = 'RNA'
obj$keep = FALSE
obj$keep[ obj$seurat_clusters %in% c(0,1,3,4) ] = TRUE
obj <- subset(obj, subset = keep == TRUE )

# re-generate seurat object and remove extraneous data
options(Seurat.object.assay.version = "v5")
obj <- CreateSeuratObject(counts = obj[["RNA"]]$counts, meta.data = obj@meta.data)
saveRDS( obj, 'n0h/union_fibro4.rds', compress = F )

# re-perform BPcells
counts.mat = obj[['RNA']]$counts
counts.mat = convert_matrix_type( counts.mat, type="uint32_t" )
counts.mat = as( counts.mat, Class='dgCMatrix' )
write_matrix_dir(mat = counts.mat, dir = "n0h/union_fibro4_counts", overwrite = T )
counts.mat <- open_matrix_dir(     dir = "n0h/union_fibro4_counts" )
obj = CreateSeuratObject( counts.mat, meta.data = obj@meta.data )
saveRDS( obj, 'n0h/union_fibro4_counts.rds', compress = F )

# run downstream analysis
obj <- NormalizeData(obj)
obj[["RNA"]] <- split(obj[["RNA"]], f = obj$study)
saveRDS( obj, 'n0h/union_fibro4_split.rds', compress = F )
obj <- FindVariableFeatures(obj, verbose = T )
obj <- ScaleData(obj, verbose = T)
obj <- RunPCA(obj, npcs = 30, verbose = T)
obj <- IntegrateLayers(obj, method = HarmonyIntegration, orig.reduction = "pca", new.reduction = "harmony" )
obj <- RunUMAP(obj, reduction = "harmony", dims = 1:30, reduction.name = "umap.harmony", return.model = T, min.dist = 0.001)
obj <- JoinLayers(obj)
saveRDS( obj, 'n0h/union_fibro4_harmony_umap_joined.rds', compress = F ) 

obj <- FindNeighbors(obj, reduction = "harmony", dims = 1:30)
for( res in seq(0.04, 0.3, 0.01) ) { # choose 0.16
  obj <- FindClusters(obj, resolution = res )
  saveRDS( obj, paste0( 'n0h/union_fibro4_harmony_umap_clusters_', res, '.rds' ), compress = F ) 
}

# join layers and find clusters
obj <- readRDS("n0h/union_fibro4_harmony_umap_clusters_0.16.rds")
DefaultAssay(obj) <- "RNA"
obj[["RNA"]]$data <- as(object = obj[["RNA"]]$data, Class = "dgCMatrix")
markers <- FindAllMarkers(obj, max.cells.per.ident = 1000, only.pos = T)
write.csv(markers, "n0h/union_fibro4_harmony_umap_clusters_0.16_markers.csv")
markers$rank <- markers$pct.1 / (markers$pct.2 + 0.001) * markers$avg_log2FC
topmarkers <- markers %>%
  group_by(cluster) %>%
  slice_max(n = 100, order_by = rank)
write.csv(topmarkers, "n0h/union_fibro4_harmony_umap_clusters_0.16_topmarkers.csv")

obj <- readRDS("n0h/union_fibro4_harmony_umap_clusters_0.16.rds")
# already joined
obj[["RNA"]]$counts <- as(object = obj[["RNA"]]$counts, Class = "dgCMatrix")
obj[["RNA"]]$data <- as(object = obj[["RNA"]]$data, Class = "dgCMatrix")
obj[["RNA"]]$scale.data <- as(object = obj[["RNA"]]$scale.data, Class = "dgCMatrix")

saveRDS(obj, "n0h/union_fibro4_harmony_umap_clusters_0.16_sparse.rds", compress = F)

# manual annotations for entire dataset
obj <- subset(obj, subset = seurat_clusters %in% c(0:7))
obj$manual.annotation <- "tbd"
obj$manual.annotation[obj$seurat_clusters == 0] <-  "ssCAF PI16+"
obj$manual.annotation[obj$seurat_clusters == 1] <-  "mCAF LRRC15+"
obj$manual.annotation[obj$seurat_clusters == 2] <-  "iCAF CXCL8+"
obj$manual.annotation[obj$seurat_clusters == 3] <-  "mCAF LRRC15+"
obj$manual.annotation[obj$seurat_clusters == 4] <-  "ssCAF CXCL14+"
obj$manual.annotation[obj$seurat_clusters == 5] <-  "mCAF COL4A1+"
obj$manual.annotation[obj$seurat_clusters == 6] <-  "apCAF CD74+"
obj$manual.annotation[obj$seurat_clusters == 7] <-  "iCAF ISG15+"

saveRDS(obj, "n0h/final_figures/human_fibro4.rds", compress = F)

### Endothelial Iteration 1  ----

obj     = readRDS( "/oak/stanford/groups/longaker/CAF_META/newHuman/n0d/humanMerged_norm_findvar_sketch1000_harmony_umap_res0.06_project_joined.rds" )
obj_har = readRDS( "/oak/stanford/groups/longaker/CAF_META/newHuman/n0d/humanMerged_norm_findvar_sketch1000_harmony_umap_res0.06_clusterList.rds" )
obj_scv = readRDS( "/oak/stanford/groups/longaker/CAF_META/newHuman/n0d/humanMerged_norm_findvar_sketch1000_scvi_umap_res0.09_clusterList.rds"    )
obj_mnn = readRDS( "/oak/stanford/groups/longaker/CAF_META/newHuman/n0d/humanMerged_norm_findvar_sketch1000_mnn_umap_res0.08_clusterList.rds"    )

ii = ( colnames(obj) %in% colnames(obj)[ (obj_har == 4) ] ) |
  ( colnames(obj) %in% colnames(obj)[ (obj_scv == 5) ] ) |
  ( colnames(obj) %in% colnames(obj)[ (obj_mnn == 4) ] )

DefaultAssay(obj) = 'RNA'
obj$keep = FALSE
obj$keep[ ii ] = TRUE
obj <- subset(obj, subset = keep == TRUE )

# re-generate seurat object and remove extraneous data
options(Seurat.object.assay.version = "v5")
obj <- CreateSeuratObject(counts = obj[["RNA"]]$counts, meta.data = obj@meta.data)
saveRDS( obj, 'n0f/union_endo.rds', compress = F )

# re-perform BPcells
counts.mat = obj[['RNA']]$counts
counts.mat = convert_matrix_type( counts.mat, type="uint32_t" )
counts.mat = as( counts.mat, Class='dgCMatrix' )
write_matrix_dir(mat = counts.mat, dir = "n0f/union_endo_counts", overwrite = T )
counts.mat <- open_matrix_dir(     dir = "n0f/union_endo_counts" )
obj = CreateSeuratObject( counts.mat, meta.data = obj@meta.data )
saveRDS( obj, 'n0f/union_endo_counts.rds', compress = F )

# run downstream analysis
obj <- NormalizeData(obj)
obj[["RNA"]] <- split(obj[["RNA"]], f = obj$study)
obj <- FindVariableFeatures(obj, verbose = T )
obj <- ScaleData(obj, verbose = T)
obj <- RunPCA(obj, npcs = 30, verbose = T)
obj <- IntegrateLayers(obj, method = HarmonyIntegration, orig.reduction = "pca", new.reduction = "harmony" )
obj <- RunUMAP(obj, reduction = "harmony", dims = 1:30, reduction.name = "umap.harmony", return.model = T, min.dist = 0.001)
saveRDS( obj, 'n0f/union_endo_harmony_umap.rds', compress = F ) 

# find clusters
obj <- FindNeighbors(obj, reduction = "harmony", dims = 1:30)
for( res in seq( 0.03, 0.1, 0.01) ) { # choose 0.08
  obj <- FindClusters(obj, resolution = res )
  saveRDS( obj, paste0( 'n0f/union_endo_harmony_clusters_', res, '.rds' ), compress = F ) 
}

# join layers and find markers
obj <- readRDS("n0f/union_endo_harmony_clusters_0.08.rds")
DefaultAssay(obj) <- "RNA"
obj <- JoinLayers(obj)
saveRDS(obj, "n0f/union_endo_harmony_clusters_0.08_joined.rds", compress = F)
obj[["RNA"]]$data <- as(object = obj[["RNA"]]$data, Class = "dgCMatrix")
markers <- FindAllMarkers(obj, max.cells.per.ident = 1000, only.pos = T)
write.csv(markers, "n0f/union_endo_harmony_clusters_0.08_markers.csv")
markers$rank <- markers$pct.1 / (markers$pct.2 + 0.001) * markers$avg_log2FC
topmarkers <- markers %>%
  group_by(cluster) %>%
  slice_max(n = 100, order_by = rank)
write.csv(topmarkers, "n0f/union_endo_harmony_clusters_0.08_topmarkers.csv")

### Endothelial Iteration 2 ----

obj     = readRDS( "n0f/union_endo_harmony_clusters_0.08_joined.rds" )

# Endothelial
DefaultAssay(obj) = 'RNA'
obj$keep = FALSE
obj$keep[ obj$seurat_clusters %in% c(0, 1, 7) ] = TRUE
obj <- subset(obj, subset = keep == TRUE )

# re-generate seurat object and remove extraneous data
options(Seurat.object.assay.version = "v5")
obj <- CreateSeuratObject(counts = obj[["RNA"]]$counts, meta.data = obj@meta.data)
saveRDS( obj, 'n0f/it2_union_endo.rds', compress = F )

# re-perform BPcells
counts.mat = obj[['RNA']]$counts
counts.mat = convert_matrix_type( counts.mat, type="uint32_t" )
counts.mat = as( counts.mat, Class='dgCMatrix' )
write_matrix_dir(mat = counts.mat, dir = "n0f/it2_union_endo_counts", overwrite = T )
counts.mat <- open_matrix_dir(     dir = "n0f/it2_union_endo_counts" )
obj = CreateSeuratObject( counts.mat, meta.data = obj@meta.data )
saveRDS( obj, 'n0f/it2_union_endo_counts.rds', compress = F )

# run downstream analysis
obj <- NormalizeData(obj)
obj[["RNA"]] <- split(obj[["RNA"]], f = obj$study)
obj <- FindVariableFeatures(obj, verbose = T )
obj <- ScaleData(obj, verbose = T)
obj <- RunPCA(obj, npcs = 30, verbose = T)
obj <- IntegrateLayers(obj, method = HarmonyIntegration, orig.reduction = "pca", new.reduction = "harmony" )
obj <- RunUMAP(obj, reduction = "harmony", dims = 1:30, reduction.name = "umap.harmony", return.model = T, min.dist = 0.001)
saveRDS( obj, 'n0f/it2_union_endo_harmony_umap.rds', compress = F ) 

# find clusters
obj <- FindNeighbors(obj, reduction = "harmony", dims = 1:30)
for( res in seq( 0.04, 0.1, 0.01 ) ) { # choose 0.07
  obj <- FindClusters(obj, resolution = res )
  saveRDS( obj, paste0( 'n0f/it2_union_endo_harmony_clusters_', res, '.rds' ), compress = F ) 
}

# join layers and find clusters
obj <- readRDS("n0f/it2_union_endo_harmony_clusters_0.07.rds")
DefaultAssay(obj) <- "RNA"
obj <- JoinLayers(obj)
saveRDS(obj, "n0f/it2_union_endo_harmony_clusters_0.07_joined.rds", compress = F)
obj[["RNA"]]$data <- as(object = obj[["RNA"]]$data, Class = "dgCMatrix")
markers <- FindAllMarkers(obj, max.cells.per.ident = 1000, only.pos = T)
write.csv(markers, "n0f/it2_union_endo_harmony_clusters_0.07_markers.csv")
markers$rank <- markers$pct.1 / (markers$pct.2 + 0.001) * markers$avg_log2FC
topmarkers <- markers %>%
  group_by(cluster) %>%
  slice_max(n = 100, order_by = rank)
write.csv(topmarkers, "n0f/it2_union_endo_harmony_clusters_0.07_topmarkers.csv")

# convert to sparse matrices
obj <- readRDS("n0f/it2_union_endo_harmony_clusters_0.07_joined.rds")
obj[['RNA']]$counts <- as( obj[['RNA']]$counts, Class='dgCMatrix' )
obj[['RNA']]$data <- as( obj[['RNA']]$data, Class='dgCMatrix' )
obj[['RNA']]$scale.data <- as( obj[['RNA']]$scale.data, Class='dgCMatrix' )
saveRDS(obj, "n0f/it2_union_endo_harmony_clusters_0.07_joined_sparse.rds", compress = F)

# manual annotations for entire dataset
obj <- subset(obj, subset = seurat_clusters %in% (0:5))
obj$manual.annotation <- "tbd"
obj$manual.annotation[obj$seurat_clusters == 0] <- "Vein ACKR1+"
obj$manual.annotation[obj$seurat_clusters == 1] <- "Tip ESM1+"
obj$manual.annotation[obj$seurat_clusters == 2] <- "Capillary CA4+"
obj$manual.annotation[obj$seurat_clusters == 3] <- "Artery GJA5+"
obj$manual.annotation[obj$seurat_clusters == 4] <- "Lymphatic PROX1+"
obj$manual.annotation[obj$seurat_clusters == 5] <- "EC ISG15+"

saveRDS(obj, "n0f/final_figures/human_endo.rds", compress = F)

### Mural Iteration 1  ----

obj     = readRDS( "/oak/stanford/groups/longaker/CAF_META/newHuman/n0d/humanMerged_norm_findvar_sketch1000_harmony_umap_res0.06_project_joined.rds" )
obj_har = readRDS( "/oak/stanford/groups/longaker/CAF_META/newHuman/n0d/humanMerged_norm_findvar_sketch1000_harmony_umap_res0.06_clusterList.rds" )
obj_scv = readRDS( "/oak/stanford/groups/longaker/CAF_META/newHuman/n0d/humanMerged_norm_findvar_sketch1000_scvi_umap_res0.09_clusterList.rds"    )
obj_mnn = readRDS( "/oak/stanford/groups/longaker/CAF_META/newHuman/n0d/humanMerged_norm_findvar_sketch1000_mnn_umap_res0.08_clusterList.rds"    )

ii = ( colnames(obj) %in% colnames(obj)[ (obj_har == 5) ] ) |
  ( colnames(obj) %in% colnames(obj)[ (obj_scv == 6) ] ) |
  ( colnames(obj) %in% colnames(obj)[ (obj_mnn == 5) ] )

DefaultAssay(obj) = 'RNA'
obj$keep = FALSE
obj$keep[ ii ] = TRUE
obj <- subset(obj, subset = keep == TRUE )

# re-generate seurat object and remove extraneous data
options(Seurat.object.assay.version = "v5")
obj <- CreateSeuratObject(counts = obj[["RNA"]]$counts, meta.data = obj@meta.data)
saveRDS( obj, 'n0f/union_mural.rds', compress = F )

# re-perform BPcells
counts.mat = obj[['RNA']]$counts
counts.mat = convert_matrix_type( counts.mat, type="uint32_t" )
counts.mat = as( counts.mat, Class='dgCMatrix' )
write_matrix_dir(mat = counts.mat, dir = "n0f/union_mural_counts", overwrite = T )
counts.mat <- open_matrix_dir(     dir = "n0f/union_mural_counts" )
obj = CreateSeuratObject( counts.mat, meta.data = obj@meta.data )
saveRDS( obj, 'n0f/union_mural_counts.rds', compress = F )

# run downstream analysis
obj <- NormalizeData(obj)
obj[["RNA"]] <- split(obj[["RNA"]], f = obj$study)
obj <- FindVariableFeatures(obj, verbose = T )
obj <- ScaleData(obj, verbose = T)
obj <- RunPCA(obj, npcs = 30, verbose = T)
obj <- IntegrateLayers(obj, method = HarmonyIntegration, orig.reduction = "pca", new.reduction = "harmony" )
obj <- RunUMAP(obj, reduction = "harmony", dims = 1:30, reduction.name = "umap.harmony", return.model = T, min.dist = 0.001)
saveRDS( obj, 'n0f/union_mural_harmony_umap.rds', compress = F ) 

# find clusters
obj <- FindNeighbors(obj, reduction = "harmony", dims = 1:30)
for( res in seq( 0.04, 0.1, 0.01 ) ) {
  obj <- FindClusters(obj, resolution = res )
  saveRDS( obj, paste0( 'n0f/union_mural_harmony_clusters_', res, '.rds' ), compress = F ) 
}

# join and find markers
obj <- readRDS("n0f/union_mural_harmony_clusters_0.05.rds")
DefaultAssay(obj) <- "RNA"
obj <- JoinLayers(obj)
saveRDS(obj, "n0f/union_mural_harmony_clusters_0.05_joined.rds", compress = F)
obj[["RNA"]]$data <- as(object = obj[["RNA"]]$data, Class = "dgCMatrix")
markers <- FindAllMarkers(obj, max.cells.per.ident = 1000, only.pos = T)
write.csv(markers, "n0f/union_mural_harmony_clusters_0.05_markers.csv")
markers$rank <- markers$pct.1 / (markers$pct.2 + 0.001) * markers$avg_log2FC
topmarkers <- markers %>%
  group_by(cluster) %>%
  slice_max(n = 100, order_by = rank)
write.csv(topmarkers, "n0f/union_mural_harmony_clusters_0.05_topmarkers.csv")

### Mural Iteration 2  ----

obj     = readRDS( "n0f/union_mural_harmony_clusters_0.05_joined.rds" )

# Pericytes
DefaultAssay(obj) = 'RNA'
obj$keep = FALSE
obj$keep[ obj$seurat_clusters %in% c(0, 1) ] = TRUE
obj <- subset(obj, subset = keep == TRUE )

# re-generate seurat object and remove extraneous data
options(Seurat.object.assay.version = "v5")
obj <- CreateSeuratObject(counts = obj[["RNA"]]$counts, meta.data = obj@meta.data)
saveRDS( obj, 'n0f/it2_union_mural.rds', compress = F )

# re-perform BPcells
counts.mat = obj[['RNA']]$counts
counts.mat = convert_matrix_type( counts.mat, type="uint32_t" )
counts.mat = as( counts.mat, Class='dgCMatrix' )
write_matrix_dir(mat = counts.mat, dir = "n0f/it2_union_mural_counts", overwrite = T )
counts.mat <- open_matrix_dir(     dir = "n0f/it2_union_mural_counts" )
obj = CreateSeuratObject( counts.mat, meta.data = obj@meta.data )
saveRDS( obj, 'n0f/it2_union_mural_counts.rds', compress = F )

# run downstream analysis
obj <- NormalizeData(obj)
obj[["RNA"]] <- split(obj[["RNA"]], f = obj$study)
obj <- FindVariableFeatures(obj, verbose = T )
obj <- ScaleData(obj, verbose = T)
obj <- RunPCA(obj, npcs = 30, verbose = T)
obj <- IntegrateLayers(obj, method = HarmonyIntegration, orig.reduction = "pca", new.reduction = "harmony" )
obj <- RunUMAP(obj, reduction = "harmony", dims = 1:30, reduction.name = "umap.harmony", return.model = T, min.dist = 0.001)
saveRDS( obj, 'n0f/it2_union_mural_harmony_umap.rds', compress = F ) 

# find clusters
obj <- FindNeighbors(obj, reduction = "harmony", dims = 1:30)
for( res in seq( 0.04, 0.1, 0.01 ) ) {
  obj <- FindClusters(obj, resolution = res )
  saveRDS( obj, paste0( 'n0f/it2_union_mural_harmony_clusters_', res, '.rds' ), compress = F ) 
}

# join and find markers
obj <- readRDS("n0f/it2_union_mural_harmony_clusters_0.08.rds")
DefaultAssay(obj) <- "RNA"
obj <- JoinLayers(obj)
saveRDS(obj, "n0f/it2_union_mural_harmony_clusters_0.08_joined.rds", compress = F)
obj[["RNA"]]$data <- as(object = obj[["RNA"]]$data, Class = "dgCMatrix")
markers <- FindAllMarkers(obj, max.cells.per.ident = 10000, only.pos = T)
write.csv(markers, "n0f/it2_union_mural_harmony_clusters_0.08_markers.csv")
markers$rank <- markers$pct.1 / (markers$pct.2 + 0.001) * markers$avg_log2FC
topmarkers <- markers %>%
  group_by(cluster) %>%
  slice_max(n = 100, order_by = rank)
write.csv(topmarkers, "n0f/it2_union_mural_harmony_clusters_0.08_topmarkers.csv")

### Mural Iteration 3 ----

obj     = readRDS( "n0f/it2_union_mural_harmony_clusters_0.08_joined.rds" )

DefaultAssay(obj) = 'RNA'
obj$keep = FALSE
obj$keep[ obj$seurat_clusters %in% c(0:4) ] = TRUE
obj <- subset(obj, subset = keep == TRUE )

# re-generate seurat object and remove extraneous data
options(Seurat.object.assay.version = "v5")
obj <- CreateSeuratObject(counts = obj[["RNA"]]$counts, meta.data = obj@meta.data)
saveRDS( obj, 'n0f/it3_union_mural.rds', compress = F )

# re-perform BPcells
counts.mat = obj[['RNA']]$counts
counts.mat = convert_matrix_type( counts.mat, type="uint32_t" )
counts.mat = as( counts.mat, Class='dgCMatrix' )
write_matrix_dir(mat = counts.mat, dir = "n0f/it3_union_mural_counts", overwrite = T )
counts.mat <- open_matrix_dir(     dir = "n0f/it3_union_mural_counts" )
obj = CreateSeuratObject( counts.mat, meta.data = obj@meta.data )
saveRDS( obj, 'n0f/it3_union_mural_counts.rds', compress = F )

# run downstream analysis
obj <- NormalizeData(obj)
obj[["RNA"]] <- split(obj[["RNA"]], f = obj$study)
obj <- FindVariableFeatures(obj, verbose = T )
obj <- ScaleData(obj, verbose = T)
obj <- RunPCA(obj, npcs = 30, verbose = T)
obj <- IntegrateLayers(obj, method = HarmonyIntegration, orig.reduction = "pca", new.reduction = "harmony" )
obj <- RunUMAP(obj, reduction = "harmony", dims = 1:30, reduction.name = "umap.harmony", return.model = T, min.dist = 0.001)
saveRDS( obj, 'n0f/it3_union_mural_harmony_umap.rds', compress = F ) 

# find clusters
obj <- FindNeighbors(obj, reduction = "harmony", dims = 1:30)
for( res in seq(0.03, 0.25, 0.01) ) { # choose 0.17
  obj <- FindClusters(obj, resolution = res )
  saveRDS( obj, paste0( 'n0f/it3_union_mural_harmony_clusters_', res, '.rds' ), compress = F ) 
}

# join and find markers
obj <- readRDS("n0f/it3_union_mural_harmony_clusters_0.17.rds")
DefaultAssay(obj) <- "RNA"
obj <- JoinLayers(obj)
saveRDS(obj, "n0f/it3_union_mural_harmony_clusters_0.17_joined.rds", compress = F)
obj[["RNA"]]$data <- as(object = obj[["RNA"]]$data, Class = "dgCMatrix")
obj <- subset(obj, subset = seurat_clusters %in% c(0:12))
markers <- FindAllMarkers(obj, max.cells.per.ident = 1000, only.pos = T)
write.csv(markers, "n0f/it3_union_mural_harmony_clusters_0.17_markers.csv")
markers$rank <- markers$pct.1 / (markers$pct.2 + 0.001) * markers$avg_log2FC
topmarkers <- markers %>%
  group_by(cluster) %>%
  slice_max(n = 100, order_by = rank)
write.csv(topmarkers, "n0f/it3_union_mural_harmony_clusters_0.17_topmarkers.csv")

# un-BP-cellize for further distribution
obj <- readRDS("n0f/it3_union_mural_harmony_clusters_0.17_joined.rds")
obj[['RNA']]$counts <- as( obj[['RNA']]$counts, Class='dgCMatrix' )
obj[['RNA']]$data <- as( obj[['RNA']]$data, Class='dgCMatrix' )
obj[['RNA']]$scale.data <- as( obj[['RNA']]$scale.data, Class='dgCMatrix' )
saveRDS(obj, "n0f/it3_union_mural_harmony_clusters_0.17_joined_sparse.rds", compress = F)

# manual annotations for entire dataset
obj <- subset(obj, subset = seurat_clusters %in% (0:6))
obj$manual.annotation <- "tbd"
obj$manual.annotation[obj$seurat_clusters == 0] <- "Pericyte CD248+"
obj$manual.annotation[obj$seurat_clusters == 1] <- "SMC vascular RERGL+"
obj$manual.annotation[obj$seurat_clusters == 2] <- "Pericyte CCL2+"
obj$manual.annotation[obj$seurat_clusters == 3] <- "SMC vascular WFDC1+"
obj$manual.annotation[obj$seurat_clusters == 4] <- "SMC visceral DES+"
obj$manual.annotation[obj$seurat_clusters == 5] <- "SMC vascular WFDC1+"
obj$manual.annotation[obj$seurat_clusters == 6] <- "Pericyte ISG15+"

saveRDS(obj, "n0f/final_figures/human_mural.rds", compress = F)

### Lymph Iteration 1 ----
obj     = readRDS( "/oak/stanford/groups/longaker/CAF_META/newHuman/n0d/humanMerged_norm_findvar_sketch1000_harmony_umap_res0.06_project_joined.rds" )
obj_har = readRDS( "/oak/stanford/groups/longaker/CAF_META/newHuman/n0d/humanMerged_norm_findvar_sketch1000_harmony_umap_res0.06_clusterList.rds" )
obj_scv = readRDS( "/oak/stanford/groups/longaker/CAF_META/newHuman/n0d/humanMerged_norm_findvar_sketch1000_scvi_umap_res0.09_clusterList.rds"    )
obj_mnn = readRDS( "/oak/stanford/groups/longaker/CAF_META/newHuman/n0d/humanMerged_norm_findvar_sketch1000_mnn_umap_res0.08_clusterList.rds"    )

# Lymph
ii = ( names(obj_har) %in% names(obj_har)[ (obj_har == 1) | (obj_har == 6) | (obj_har == 7)] ) |
  ( names(obj_scv) %in% names(obj_scv)[ (obj_scv == 0) | (obj_scv == 9) | (obj_scv == 10) ] ) |
  ( names(obj_mnn) %in% names(obj_mnn)[ (obj_mnn == 2) | (obj_mnn == 6) | (obj_mnn == 9) ] )

DefaultAssay(obj) = 'RNA'
obj$keep = FALSE
obj$keep[ ii ] = TRUE
obj <- subset(obj, subset = keep == TRUE )

# re-generate seurat object and remove extraneous data
options(Seurat.object.assay.version = "v5")
obj <- CreateSeuratObject(counts = obj[["RNA"]]$counts, meta.data = obj@meta.data)
saveRDS( obj, 'o0c/union_lymph.rds', compress = F )

# re-perform BPcells
counts.mat = obj[['RNA']]$counts
counts.mat = convert_matrix_type( counts.mat, type="uint32_t" )
counts.mat = as( counts.mat, Class='dgCMatrix' )
write_matrix_dir(mat = counts.mat, dir = "o0c/union_lymph_counts", overwrite = T )
counts.mat <- open_matrix_dir(     dir = "o0c/union_lymph_counts" )
obj = CreateSeuratObject( counts.mat, meta.data = obj@meta.data )
saveRDS( obj, 'o0c/union_lymph_counts.rds', compress = F )

# run downstream analysis
obj <- NormalizeData(obj)
obj[["RNA"]] <- split(obj[["RNA"]], f = obj$study)
obj <- FindVariableFeatures(obj, verbose = T)
obj <- SketchData(object = obj, ncells = 1000, method = "LeverageScore", sketched.assay = "sketch")
DefaultAssay(obj) <- "sketch"
obj <- FindVariableFeatures(obj, verbose = T)
obj <- ScaleData(obj, verbose = F)
obj <- RunPCA(obj, npcs = 30, verbose = F)
saveRDS(obj, "o0c/lymph_it1_var_sketch_pca.rds", compress = F)

obj <- IntegrateLayers(obj, method = HarmonyIntegration, orig.reduction = "pca", new.reduction = "harmony", kmeans_init_nstart=20, kmeans_init_iter_max=5000, verbose = T)
obj <- RunUMAP(obj, reduction = "harmony", dims = 1:30, reduction.name = "umap.harmony", return.model = T, min.dist = 0.001)
saveRDS(obj, "o0c/lymph_it1_var_sketch_pca_harmony_umap.rds", compress = F)

obj <- FindNeighbors(obj, reduction = "harmony", dims = 1:30)
for (i in seq(0.03, 0.15, 0.01)) { # 0.04
  obj <- FindClusters(obj, resolution = i)
  saveRDS(obj, paste0("o0c/lymph_it1_var_sketch_pca_harmony_umap_res",i,".rds"), compress = FALSE)
}

# find markers
obj <- readRDS("o0c/lymph_it1_var_sketch_pca_harmony_umap_res0.04.rds")
DefaultAssay(obj) <- "sketch"
obj <- JoinLayers(obj)
markers <- FindAllMarkers(obj, max.cells.per.ident = 1000, only.pos = T)
write.csv(markers, "o0c/lymph_it1_var_sketch_pca_harmony_umap_res0.04_markers.csv")
markers$rank <- markers$pct.1 / (markers$pct.2 + 0.001) * markers$avg_log2FC
topmarkers <- markers %>%
  group_by(cluster) %>%
  slice_max(n = 100, order_by = rank)
write.csv(topmarkers, "o0c/lymph_it1_var_sketch_pca_harmony_umap_res0.04_topmarkers.csv")

# project to full dataset
obj <- readRDS("o0c/lymph_it1_var_sketch_pca_harmony_umap_res0.04.rds") # needs to be split (non-joined)
obj <- ProjectIntegration(object = obj, reduction = "harmony")
options(future.globals.maxSize = 8000 * 1024^2)
obj <- ProjectData(object = obj, sketched.reduction = "harmony.full", full.reduction = "harmony.full", umap.model = "umap.harmony", dims = 1:30, refdata = list(harmony.cluster.full = "seurat_clusters"))
saveRDS(obj, "o0c/lymph_it1_var_sketch_pca_harmony_umap_res0.04_project.rds",compress = F)
saveRDS(obj$harmony.cluster.full, "o0c/lymph_it1_var_sketch_pca_harmony_umap_res0.04_project_clusterList.rds")

# join layers
DefaultAssay(obj) <- "sketch"
obj <- JoinLayers(obj)
DefaultAssay(obj) <- "RNA"
obj <- JoinLayers(obj)
saveRDS(obj, "o0c/lymph_it1_var_sketch_pca_harmony_umap_res0.04_project_joined.rds", compress = F)

### Lymph Iteration 2 ----

obj <- readRDS( "o0c/lymph_it1_var_sketch_pca_harmony_umap_res0.04_project_joined.rds" )

# lymphoid
DefaultAssay(obj) = 'RNA'
obj$keep = FALSE
obj$keep[ obj$harmony.cluster.full %in% c(0, 2, 5) ] = TRUE
obj <- subset(obj, subset = keep == TRUE )

# re-generate seurat object and remove extraneous data
options(Seurat.object.assay.version = "v5")
obj <- CreateSeuratObject(counts = obj[["RNA"]]$counts, meta.data = obj@meta.data)
saveRDS( obj, 'o0c/it2_union_lymph.rds', compress = F )

# re-perform BPcells
counts.mat = obj[['RNA']]$counts
counts.mat = convert_matrix_type( counts.mat, type="uint32_t" )
counts.mat = as( counts.mat, Class='dgCMatrix' )
write_matrix_dir(mat = counts.mat, dir = "o0c/it2_union_lymph_counts", overwrite = T )
counts.mat <- open_matrix_dir(     dir = "o0c/it2_union_lymph_counts" )
obj = CreateSeuratObject( counts.mat, meta.data = obj@meta.data )
saveRDS( obj, 'o0c/it2_union_lymph_counts.rds', compress = F )

# run downstream analysis
obj <- NormalizeData(obj)
obj[["RNA"]] <- split(obj[["RNA"]], f = obj$study)
obj <- FindVariableFeatures(obj, verbose = T)
obj <- SketchData(object = obj, ncells = 1000, method = "LeverageScore", sketched.assay = "sketch")
DefaultAssay(obj) <- "sketch"
obj <- FindVariableFeatures(obj, verbose = T)
obj <- ScaleData(obj, verbose = F)
obj <- RunPCA(obj, npcs = 30, verbose = F)
saveRDS(obj, "o0c/lymph_it2_var_sketch_pca.rds", compress = F)


obj <- IntegrateLayers(obj, method = HarmonyIntegration, orig.reduction = "pca", new.reduction = "harmony", kmeans_init_nstart=20, kmeans_init_iter_max=5000, verbose = T)
obj <- RunUMAP(obj, reduction = "harmony", dims = 1:30, reduction.name = "umap.harmony", return.model = T, min.dist = 0.001)
saveRDS(obj, "o0c/lymph_it2_var_sketch_pca_harmony_umap.rds", compress = F)

obj <- FindNeighbors(obj, reduction = "harmony", dims = 1:30)
for (i in c(0.1, 0.03, 0.2, 0.3, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.11, 0.12, 0.13, 0.14, 0.15)) { #choose 0.3
  obj <- FindClusters(obj, resolution = i)
  saveRDS(obj, paste0("o0c/lymph_it2_var_sketch_pca_harmony_umap_res",i,".rds"), compress = FALSE)
}

# find markers
obj <- readRDS("o0c/lymph_it2_var_sketch_pca_harmony_umap_res0.3.rds")
DefaultAssay(obj) <- "sketch"
obj <- JoinLayers(obj)
markers <- FindAllMarkers(obj, max.cells.per.ident = 1000, only.pos = T)
write.csv(markers, "o0c/lymph_it2_var_sketch_pca_harmony_umap_res0.3_markers.csv")
markers$rank <- markers$pct.1 / (markers$pct.2 + 0.001) * markers$avg_log2FC
topmarkers <- markers %>%
  group_by(cluster) %>%
  slice_max(n = 100, order_by = rank)
write.csv(topmarkers, "o0c/lymph_it2_var_sketch_pca_harmony_umap_res0.3_topmarkers.csv")

# project to full dataset
obj <- readRDS("o0c/lymph_it2_var_sketch_pca_harmony_umap_res0.3.rds") # needs to be split (non-joined)
obj <- ProjectIntegration(object = obj, reduction = "harmony")
options(future.globals.maxSize = 8000 * 1024^2)
obj <- ProjectData(object = obj, sketched.reduction = "harmony.full", full.reduction = "harmony.full", umap.model = "umap.harmony", dims = 1:30, refdata = list(harmony.cluster.full = "seurat_clusters"))
saveRDS(obj, "o0c/lymph_it2_var_sketch_pca_harmony_umap_res0.3_project.rds",compress = F)
saveRDS(obj$harmony.cluster.full, "o0c/lymph_it2_var_sketch_pca_harmony_umap_res0.3_project_clusterList.rds")

# join layers
DefaultAssay(obj) <- "sketch"
obj <- JoinLayers(obj)
DefaultAssay(obj) <- "RNA"
obj <- JoinLayers(obj)
saveRDS(obj, "o0c/lymph_it2_var_sketch_pca_harmony_umap_res0.3_project_joined.rds", compress = F)

### Lymph Iteration 3 ----

obj <- readRDS( "o0c/lymph_it2_var_sketch_pca_harmony_umap_res0.3_project_joined.rds" )

# lymphoid
DefaultAssay(obj) = 'RNA'
obj$keep = FALSE
obj$keep[ obj$harmony.cluster.full %in% c(0:5, 7:8, 11) ] = TRUE
obj <- subset(obj, subset = keep == TRUE )

# re-generate seurat object and remove extraneous data
options(Seurat.object.assay.version = "v5")
obj <- CreateSeuratObject(counts = obj[["RNA"]]$counts, meta.data = obj@meta.data)
saveRDS( obj, 'o0c/it3_union_lymph.rds', compress = F )

# re-perform BPcells
counts.mat = obj[['RNA']]$counts
counts.mat = convert_matrix_type( counts.mat, type="uint32_t" )
counts.mat = as( counts.mat, Class='dgCMatrix' )
write_matrix_dir(mat = counts.mat, dir = "o0c/it3_union_lymph_counts", overwrite = T )
counts.mat <- open_matrix_dir(     dir = "o0c/it3_union_lymph_counts" )
obj = CreateSeuratObject( counts.mat, meta.data = obj@meta.data )
saveRDS( obj, 'o0c/it3_union_lymph_counts.rds', compress = F )

# run downstream analysis
obj <- NormalizeData(obj)
obj[["RNA"]] <- split(obj[["RNA"]], f = obj$study)
obj <- FindVariableFeatures(obj, verbose = T)
obj <- SketchData(object = obj, ncells = 1000, method = "LeverageScore", sketched.assay = "sketch")
DefaultAssay(obj) <- "sketch"
obj <- FindVariableFeatures(obj, verbose = T)
obj <- ScaleData(obj, verbose = F)
obj <- RunPCA(obj, npcs = 30, verbose = F)
saveRDS(obj, "o0c/lymph_it3_var_sketch_pca.rds", compress = F)

obj <- IntegrateLayers(obj, group.by = "study", method = HarmonyIntegration, orig.reduction = "pca", new.reduction = "harmony", kmeans_init_nstart=20, kmeans_init_iter_max=5000, verbose = T)
obj <- RunUMAP(obj, reduction = "harmony", dims = 1:30, reduction.name = "umap.harmony", return.model = T, min.dist = 0.001)
saveRDS(obj, "o0c/lymph_it3_var_sketch_pca_harmony_umap.rds", compress = F)

obj <- FindNeighbors(obj, reduction = "harmony", dims = 1:30)
for (i in seq(0.1, 1.0, 0.1) ) { # choose 0.3
  obj <- FindClusters(obj, resolution = i)
  saveRDS(obj, paste0("o0c/lymph_it3_var_sketch_pca_harmony_umap_res",i,".rds"), compress = FALSE)
}

# find markers
obj <- readRDS("o0c/lymph_it3_var_sketch_pca_harmony_umap_res0.3.rds")
DefaultAssay(obj) <- "sketch"
obj <- JoinLayers(obj)
markers <- FindAllMarkers(obj, max.cells.per.ident = 1000, only.pos = T)
write.csv(markers, "o0c/lymph_it3_var_sketch_pca_harmony_umap_res0.3_markers.csv")
markers$rank <- markers$pct.1 / (markers$pct.2 + 0.001) * markers$avg_log2FC
topmarkers <- markers %>%
  group_by(cluster) %>%
  slice_max(n = 100, order_by = rank)
write.csv(topmarkers, "o0c/lymph_it3_var_sketch_pca_harmony_umap_res0.3_topmarkers.csv")

# project to full dataset
obj <- readRDS("o0c/lymph_it3_var_sketch_pca_harmony_umap_res0.3.rds") # needs to be split (non-joined)
obj <- ProjectIntegration(object = obj, reduction = "harmony")
options(future.globals.maxSize = 8000 * 1024^2)
obj <- ProjectData(object = obj, sketched.reduction = "harmony.full", full.reduction = "harmony.full", umap.model = "umap.harmony", dims = 1:30, refdata = list(harmony.cluster.full = "seurat_clusters"))
saveRDS(obj, "o0c/lymph_it3_var_sketch_pca_harmony_umap_res0.3_project.rds",compress = F)
saveRDS(obj$harmony.cluster.full, "o0c/lymph_it3_var_sketch_pca_harmony_umap_res0.3_project_clusterList.rds")

# join layers
DefaultAssay(obj) <- "sketch"
obj <- JoinLayers(obj)
DefaultAssay(obj) <- "RNA"
obj <- JoinLayers(obj)
saveRDS(obj, "o0c/lymph_it3_var_sketch_pca_harmony_umap_res0.3_joined.rds", compress = F)

### Lymph Iteration 4 ----
obj <- readRDS( "o0c/lymph_it3_var_sketch_pca_harmony_umap_res0.3_project.rds" )

# lymphoid
DefaultAssay(obj) = 'RNA'
obj$keep = F
obj$keep[ obj$harmony.cluster.full %in% c(0:4, 6:7, 9:12, 16) ] = T
obj <- subset(obj, subset = keep == TRUE )

# re-generate seurat object and remove extraneous data
options(Seurat.object.assay.version = "v5")
obj <- CreateSeuratObject(counts = obj[["RNA"]]$counts, meta.data = obj@meta.data)
saveRDS( obj, 'o0c/it4_union_lymph.rds', compress = F )

# re-perform BPcells
counts.mat = obj[['RNA']]$counts
counts.mat = convert_matrix_type( counts.mat, type="uint32_t" )
counts.mat = as( counts.mat, Class='dgCMatrix' )
write_matrix_dir(mat = counts.mat, dir = "o0c/it4_union_lymph_counts", overwrite = T )
counts.mat <- open_matrix_dir(     dir = "o0c/it4_union_lymph_counts" )
obj = CreateSeuratObject( counts.mat, meta.data = obj@meta.data )
saveRDS( obj, 'o0c/it4_union_lymph_counts.rds', compress = F )

# run downstream analysis
obj <- NormalizeData(obj)
obj[["RNA"]] <- split(obj[["RNA"]], f = obj$study)
obj <- FindVariableFeatures(obj, verbose = T)
obj <- SketchData(object = obj, ncells = 1000, method = "LeverageScore", sketched.assay = "sketch")
DefaultAssay(obj) <- "sketch"
obj <- FindVariableFeatures(obj, verbose = T)
obj <- ScaleData(obj, verbose = F)
obj <- RunPCA(obj, npcs = 30, verbose = F)
saveRDS(obj, "o0c/lymph_it4_var_sketch_pca.rds", compress = F)

# need to specify group.by = 'study'
obj <- IntegrateLayers(obj, group.by = "study", method = HarmonyIntegration, orig.reduction = "pca", new.reduction = "harmony", kmeans_init_nstart=20, kmeans_init_iter_max=5000, verbose = T)
obj <- RunUMAP(obj, reduction = "harmony", dims = 1:30, reduction.name = "umap.harmony", return.model = T, min.dist = 0.001)
saveRDS(obj, "o0c/lymph_it4_var_sketch_pca_harmony_umap.rds", compress = F)

obj <- FindNeighbors(obj, reduction = "harmony", dims = 1:30)
for (i in seq(0.4, 1.0, 0.1) ) { # choose 0.8
  obj <- FindClusters(obj, resolution = i)
  saveRDS(obj, paste0("o0c/lymph_it4_var_sketch_pca_harmony_umap_res",i,".rds"), compress = FALSE)
}

# find markers
obj <- readRDS("o0c/lymph_it4_var_sketch_pca_harmony_umap_res0.8.rds")
DefaultAssay(obj) <- "sketch"
obj <- JoinLayers(obj)
markers <- FindAllMarkers(obj, max.cells.per.ident = 1000, only.pos = T)
write.csv(markers, "o0c/lymph_it4_var_sketch_pca_harmony_umap_res0.8_markers.csv")
markers$rank <- markers$pct.1 / (markers$pct.2 + 0.001) * markers$avg_log2FC
topmarkers <- markers %>%
  group_by(cluster) %>%
  slice_max(n = 100, order_by = rank)
write.csv(topmarkers, "o0c/lymph_it4_var_sketch_pca_harmony_umap_res0.8_topmarkers.csv")

# project to full dataset
obj <- readRDS("o0c/lymph_it4_var_sketch_pca_harmony_umap_res0.8.rds") # needs to be split (non-joined)
obj <- ProjectIntegration(object = obj, reduction = "harmony")
options(future.globals.maxSize = 8000 * 1024^2)
obj <- ProjectData(object = obj, sketched.reduction = "harmony.full", full.reduction = "harmony.full", umap.model = "umap.harmony", dims = 1:30, refdata = list(harmony.cluster.full = "seurat_clusters"))
saveRDS(obj, "o0c/lymph_it4_var_sketch_pca_harmony_umap_res0.8_project.rds",compress = F)
saveRDS(obj$harmony.cluster.full, "o0c/lymph_it4_var_sketch_pca_harmony_umap_res0.8_project_clusterList.rds")

# join layers
obj <- readRDS("o0c/lymph_it4_var_sketch_pca_harmony_umap_res0.8.rds") 
DefaultAssay(obj) <- "sketch"
obj <- JoinLayers(obj)
DefaultAssay(obj) <- "RNA"
obj <- JoinLayers(obj)
saveRDS(obj, "o0c/lymph_it4_var_sketch_pca_harmony_umap_res0.8_joined.rds", compress = F)


### Lymph Iteration 5 ----

obj <- readRDS( "o0c/lymph_it4_var_sketch_pca_harmony_umap_res0.8_project.rds" )

# lymphoid
DefaultAssay(obj) = 'RNA'
obj$keep = F
obj$keep[ obj$harmony.cluster.full %in% c(0:1, 3:6, 8:13, 15:17, 19, 21) ] = T
obj <- subset(obj, subset = keep == TRUE )

# re-generate seurat object and remove extraneous data
options(Seurat.object.assay.version = "v5")
obj <- CreateSeuratObject(counts = obj[["RNA"]]$counts, meta.data = obj@meta.data)
saveRDS( obj, 'o0c/it5_union_lymph.rds', compress = F )

# re-perform BPcells
counts.mat = obj[['RNA']]$counts
counts.mat = convert_matrix_type( counts.mat, type="uint32_t" )
counts.mat = as( counts.mat, Class='dgCMatrix' )
write_matrix_dir(mat = counts.mat, dir = "o0c/it5_union_lymph_counts", overwrite = T )
counts.mat <- open_matrix_dir(     dir = "o0c/it5_union_lymph_counts" )
obj = CreateSeuratObject( counts.mat, meta.data = obj@meta.data )
saveRDS( obj, 'o0c/it5_union_lymph_counts.rds', compress = F )

# run downstream analysis
obj <- NormalizeData(obj)
obj[["RNA"]] <- split(obj[["RNA"]], f = obj$study)
obj <- FindVariableFeatures(obj, verbose = T)
obj <- SketchData(object = obj, ncells = 1000, method = "LeverageScore", sketched.assay = "sketch")
DefaultAssay(obj) <- "sketch"
obj <- FindVariableFeatures(obj, verbose = T)
obj <- ScaleData(obj, verbose = F)
obj <- RunPCA(obj, npcs = 30, verbose = F)
saveRDS(obj, "o0c/lymph_it5_var_sketch_pca.rds", compress = F)

# need to specify group.by = 'study'
obj <- IntegrateLayers(obj, group.by = "study", method = HarmonyIntegration, orig.reduction = "pca", new.reduction = "harmony", kmeans_init_nstart=20, kmeans_init_iter_max=5000, verbose = T)
obj <- RunUMAP(obj, reduction = "harmony", dims = 1:30, reduction.name = "umap.harmony", return.model = T, min.dist = 0.001)
saveRDS(obj, "o0c/lymph_it5_var_sketch_pca_harmony_umap.rds", compress = F)

obj <- FindNeighbors(obj, reduction = "harmony", dims = 1:30)
for (i in seq(0.4, 1.0, 0.1)) { # choose 1.0
  obj <- FindClusters(obj, resolution = i)
  saveRDS(obj, paste0("o0c/lymph_it5_var_sketch_pca_harmony_umap_res",i,".rds"), compress = FALSE)
}

# find markers
obj <- readRDS("o0c/lymph_it5_var_sketch_pca_harmony_umap_res1.rds")
DefaultAssay(obj) <- "sketch"
obj <- JoinLayers(obj)
markers <- FindAllMarkers(obj, max.cells.per.ident = 1000, only.pos = T)
write.csv(markers, "o0c/lymph_it5_var_sketch_pca_harmony_umap_res1_markers.csv")
markers$rank <- markers$pct.1 / (markers$pct.2 + 0.001) * markers$avg_log2FC
topmarkers <- markers %>%
  group_by(cluster) %>%
  slice_max(n = 100, order_by = rank)
write.csv(topmarkers, "o0c/lymph_it5_var_sketch_pca_harmony_umap_res1_topmarkers.csv")

# project to full dataset
obj <- readRDS("o0c/lymph_it5_var_sketch_pca_harmony_umap_res1.rds") # needs to be split (non-joined)
obj <- ProjectIntegration(object = obj, reduction = "harmony")
options(future.globals.maxSize = 8000 * 1024^2)
obj <- ProjectData(object = obj, sketched.reduction = "harmony.full", full.reduction = "harmony.full", umap.model = "umap.harmony", dims = 1:30, refdata = list(harmony.cluster.full = "seurat_clusters"))
saveRDS(obj, "o0c/lymph_it5_var_sketch_pca_harmony_umap_res1_project.rds",compress = F)
saveRDS(obj$harmony.cluster.full, "o0c/lymph_it5_var_sketch_pca_harmony_umap_res1_project_clusterList.rds")

# join layers
obj <- readRDS("o0c/lymph_it5_var_sketch_pca_harmony_umap_res1.rds") 
DefaultAssay(obj) <- "sketch"
obj <- JoinLayers(obj)
DefaultAssay(obj) <- "RNA"
obj <- JoinLayers(obj)
saveRDS(obj, "o0c/lymph_it5_var_sketch_pca_harmony_umap_res1_joined.rds", compress = F)

### Lymph Iteration 6 ----
obj <- readRDS("o0c/lymph_it5_var_sketch_pca_harmony_umap_res1_joined.rds")
temp <- readRDS("o0c/lymph_it5_var_sketch_pca_harmony_umap_res1_project_clusterList.rds")
obj$harmony.cluster.full <- temp

# lymphoid
DefaultAssay(obj) = 'RNA'
obj$keep = FALSE
obj$keep[ obj$harmony.cluster.full %in% c(0:8, 10:13, 15, 17, 20:22) ] = TRUE
obj <- subset(obj, subset = keep == TRUE )

# re-generate seurat object and remove extraneous data
options(Seurat.object.assay.version = "v5")
obj <- CreateSeuratObject(counts = obj[["RNA"]]$counts, meta.data = obj@meta.data)
saveRDS( obj, 'o0c/lymph_it6.rds', compress = F )

# re-perform BPcells
counts.mat = obj[['RNA']]$counts
counts.mat = convert_matrix_type( counts.mat, type="uint32_t" )
counts.mat = as( counts.mat, Class='dgCMatrix' )
write_matrix_dir(mat = counts.mat, dir = "o0c/lymph_it6_counts", overwrite = T )
counts.mat <- open_matrix_dir(     dir = "o0c/lymph_it6_counts" )
obj = CreateSeuratObject( counts.mat, meta.data = obj@meta.data )
saveRDS( obj, 'o0c/lymph_it6_counts.rds', compress = F )

# run downstream analysis
obj <- NormalizeData(obj)
obj[["RNA"]] <- split(obj[["RNA"]], f = obj$study)
obj <- FindVariableFeatures(obj, verbose = T )
obj <- SketchData(object = obj, ncells = 20000, method = "LeverageScore", sketched.assay = "sketch")
DefaultAssay(obj) <- "sketch"
obj <- FindVariableFeatures(obj, verbose = T)
obj <- ScaleData(obj, verbose = F)
obj <- RunPCA(obj, npcs = 30, verbose = F)

obj <- IntegrateLayers(obj, group.by = "study", method = HarmonyIntegration, orig.reduction = "pca", new.reduction = "harmony", kmeans_init_nstart=20, kmeans_init_iter_max=5000, verbose = T)
obj <- RunUMAP(obj, reduction = "harmony", dims = 1:30, reduction.name = "umap.harmony", return.model = T, min.dist = 0.001)
saveRDS(obj, "o0c/lymph_it6_var_sketch_pca_harmony_umap.rds", compress = F)

obj <- FindNeighbors(obj, reduction = "harmony", dims = 1:30)
for (i in seq(0.4, 1.0, 0.1)) { # choose 0.6
  obj <- FindClusters(obj, resolution = i)
  saveRDS(obj, paste0("o0c/lymph_it6_var_sketch_pca_harmony_umap_res",i,".rds"), compress = FALSE)
}

# find markers
obj <- readRDS("o0c/lymph_it6_var_sketch_pca_harmony_umap_res0.6.rds")
DefaultAssay(obj) <- "sketch"
obj <- JoinLayers(obj)
markers <- FindAllMarkers(obj, max.cells.per.ident = 1000, only.pos = T)
write.csv(markers, "o0c/lymph_it6_var_sketch_pca_harmony_umap_res0.6_markers.csv")
markers$rank <- markers$pct.1 / (markers$pct.2 + 0.001) * markers$avg_log2FC
topmarkers <- markers %>%
  group_by(cluster) %>%
  slice_max(n = 100, order_by = rank)
write.csv(topmarkers, "o0c/lymph_it6_var_sketch_pca_harmony_umap_res0.6_topmarkers.csv")

# project to full dataset
obj <- readRDS("o0c/lymph_it6_var_sketch_pca_harmony_umap_res0.6.rds") # needs to be split (non-joined)
obj <- ProjectIntegration(object = obj, reduction = "harmony")
options(future.globals.maxSize = 8000 * 1024^2)
obj <- ProjectData(object = obj, sketched.reduction = "harmony.full", full.reduction = "harmony.full", umap.model = "umap.harmony", dims = 1:30, refdata = list(harmony.cluster.full = "seurat_clusters"))
saveRDS(obj, "o0c/lymph_it6_var_sketch_pca_harmony_umap_res0.6_project.rds",compress = F)
saveRDS(obj$harmony.cluster.full, "o0c/lymph_it6_var_sketch_pca_harmony_umap_res0.6_project_clusterList.rds")

# join layers
obj <- readRDS("o0c/lymph_it6_var_sketch_pca_harmony_umap_res0.6_project.rds") 
DefaultAssay(obj) <- "sketch"
obj <- JoinLayers(obj)
DefaultAssay(obj) <- "RNA"
obj <- JoinLayers(obj)
saveRDS(obj, "o0c/lymph_it6_var_sketch_pca_harmony_umap_res0.6_project_joined.rds", compress = F)

# un-BP-cellize for further distribution
obj <- readRDS("n0f/lymph_it6_var_sketch_pca_harmony_umap_res0.6_project_joined.rds")
DefaultAssay(obj) <- "RNA"
obj[['RNA']]$counts <- as( obj[['RNA']]$counts, Class='dgCMatrix' )
obj[['RNA']]$data <- as( obj[['RNA']]$data, Class='dgCMatrix' )
obj[['sketch']]$counts <- as( obj[['RNA']]$counts, Class='dgCMatrix' )
obj[['sketch']]$data <- as( obj[['RNA']]$data, Class='dgCMatrix' )
# obj[['RNA']]$scale.data <- as( obj[['RNA']]$scale.data, Class='dgCMatrix' )
saveRDS(obj, "n0f/lymph_it6_var_sketch_pca_harmony_umap_res0.6_project_joined_sparse.rds", compress = F)

# manual annotations for entire dataset
obj <- subset(obj, subset = harmony.cluster.full %in% c(0:12, 14:15, 17:20))
obj$manual.annotation <- "tbd"
obj$manual.annotation[obj$harmony.cluster.full == 0] <- "T CD8+/GZMK+"
obj$manual.annotation[obj$harmony.cluster.full == 1] <- "T CD4+/CD40LG+"
obj$manual.annotation[obj$harmony.cluster.full == 2] <- "B mem BANK1+"
obj$manual.annotation[obj$harmony.cluster.full == 3] <- "Plasma JCHAIN+"
obj$manual.annotation[obj$harmony.cluster.full == 4] <- "Treg FOXP3+"
obj$manual.annotation[obj$harmony.cluster.full == 5] <- "NK CD16+"
obj$manual.annotation[obj$harmony.cluster.full == 6] <- "NK CD56+"
obj$manual.annotation[obj$harmony.cluster.full == 7] <- "Tex CD8+/PDCD1+"
obj$manual.annotation[obj$harmony.cluster.full == 8] <- "B naive AFF3+"
obj$manual.annotation[obj$harmony.cluster.full == 9] <- "T CD8+/ISG15+"
obj$manual.annotation[obj$harmony.cluster.full == 10] <- "Tfh NR3C1+"
obj$manual.annotation[obj$harmony.cluster.full == 11] <- "Plasma JCHAIN+"
obj$manual.annotation[obj$harmony.cluster.full == 12] <- "T CD8+/MT+"
obj$manual.annotation[obj$harmony.cluster.full == 14] <- "B germinal RGS13+"
obj$manual.annotation[obj$harmony.cluster.full == 15] <- "B mem BANK1+"
obj$manual.annotation[obj$harmony.cluster.full == 17] <- "Plasma JCHAIN+"
obj$manual.annotation[obj$harmony.cluster.full == 18] <- "Plasma JCHAIN+"
obj$manual.annotation[obj$harmony.cluster.full == 19] <- "B naive AFF3+"
obj$manual.annotation[obj$harmony.cluster.full == 20] <- "Plasma JCHAIN+"

saveRDS(obj, "n0f/final_figures/human_lymph.rds", compress = F)

### Myeloid Iteration 1 ----

obj     = readRDS( "/oak/stanford/groups/longaker/CAF_META/newHuman/n0d/humanMerged_norm_findvar_sketch1000_harmony_umap_res0.06_project_joined.rds" )
obj_har = readRDS( "/oak/stanford/groups/longaker/CAF_META/newHuman/n0d/humanMerged_norm_findvar_sketch1000_harmony_umap_res0.06_clusterList.rds" )
obj_scv = readRDS( "/oak/stanford/groups/longaker/CAF_META/newHuman/n0d/humanMerged_norm_findvar_sketch1000_scvi_umap_res0.09_clusterList.rds"    )
obj_mnn = readRDS( "/oak/stanford/groups/longaker/CAF_META/newHuman/n0d/humanMerged_norm_findvar_sketch1000_mnn_umap_res0.08_clusterList.rds"    )

# Myeloid
ii = ( colnames(obj) %in% colnames(obj)[ (obj_har %in% c(2,9)) ] ) |
  ( colnames(obj) %in% colnames(obj)[ (obj_scv %in% c(2,14,15,18)) ] ) |
  ( colnames(obj) %in% colnames(obj)[ (obj_mnn %in% c(3,11)) ] )

DefaultAssay(obj) = 'RNA'
obj$keep = FALSE
obj$keep[ ii ] = TRUE
obj <- subset(obj, subset = keep == TRUE )

# re-generate seurat object and remove extraneous data
options(Seurat.object.assay.version = "v5")
obj <- CreateSeuratObject(counts = obj[["RNA"]]$counts, meta.data = obj@meta.data)
saveRDS( obj, 'n0f/union_myel.rds', compress = F )

# re-perform BPcells
counts.mat = obj[['RNA']]$counts
counts.mat = convert_matrix_type( counts.mat, type="uint32_t" )
counts.mat = as( counts.mat, Class='dgCMatrix' )
write_matrix_dir(mat = counts.mat, dir = "n0f/union_myel_counts", overwrite = T )
counts.mat <- open_matrix_dir(     dir = "n0f/union_myel_counts" )
obj = CreateSeuratObject( counts.mat, meta.data = obj@meta.data )
saveRDS( obj, 'n0f/union_myel_counts.rds', compress = F )

# run downstream analysis
obj <- NormalizeData(obj)
obj[["RNA"]] <- split(obj[["RNA"]], f = obj$study)
obj <- FindVariableFeatures(obj, verbose = T )
obj <- SketchData(object = obj, ncells = 1000, method = "LeverageScore", sketched.assay = "sketch")
DefaultAssay(obj) <- "sketch"
obj <- FindVariableFeatures(obj, verbose = T)
obj <- ScaleData(obj, verbose = F)
obj <- RunPCA(obj, npcs = 30, verbose = F)
saveRDS(obj, "n0f/union_myel_var_sketch_pca.rds", compress = F)

obj <- IntegrateLayers(obj, method = HarmonyIntegration, orig.reduction = "pca", new.reduction = "harmony", kmeans_init_nstart=20, kmeans_init_iter_max=5000, verbose = T)
obj <- RunUMAP(obj, reduction = "harmony", dims = 1:30, reduction.name = "umap.harmony", return.model = T, min.dist = 0.001)
saveRDS(obj, "n0f/union_myel_var_sketch_pca_harmony_umap.rds", compress = F)

obj <- FindNeighbors(obj, reduction = "harmony", dims = 1:30)
for (i in seq(0.03, 0.15, 0.01)) { # choose 0.04
  obj <- FindClusters(obj, resolution = i)
  saveRDS(obj, paste0("n0f/union_myel_var_sketch_pca_harmony_umap_res",i,".rds"), compress = FALSE)
}

# find markers
obj <- readRDS("n0f/union_myel_var_sketch_pca_harmony_umap_res0.04.rds")
DefaultAssay(obj) <- "sketch"
obj <- JoinLayers(obj)
markers <- FindAllMarkers(obj, max.cells.per.ident = 1000, only.pos = T)
write.csv(markers, "n0f/union_myel_var_sketch_pca_harmony_umap_res0.04_markers.csv")
markers$rank <- markers$pct.1 / (markers$pct.2 + 0.001) * markers$avg_log2FC
topmarkers <- markers %>%
  group_by(cluster) %>%
  slice_max(n = 100, order_by = rank)
write.csv(topmarkers, "n0f/union_myel_var_sketch_pca_harmony_umap_res0.04_topmarkers.csv")

# project to full dataset
obj <- readRDS("n0f/union_myel_var_sketch_pca_harmony_umap_res0.04.rds") # needs to be split (non-joined)
obj <- ProjectIntegration(object = obj, reduction = "harmony")
options(future.globals.maxSize = 8000 * 1024^2)
obj <- ProjectData(object = obj, sketched.reduction = "harmony.full", full.reduction = "harmony.full", umap.model = "umap.harmony", dims = 1:30, refdata = list(harmony.cluster.full = "seurat_clusters"))
saveRDS(obj, "n0f/union_myel_var_sketch_pca_harmony_umap_res0.04_project.rds",compress = F)
saveRDS(obj$harmony.cluster.full, "n0f/union_myel_var_sketch_pca_harmony_umap_res0.04_project_clusterList.rds")

# join layers
DefaultAssay(obj) <- "sketch"
obj <- JoinLayers(obj)
DefaultAssay(obj) <- "RNA"
obj <- JoinLayers(obj)
saveRDS(obj, "n0f/union_myel_var_sketch_pca_harmony_umap_res0.04_project_joined.rds", compress = F)


### Myeloid Iteration 2 ----

obj <- readRDS( "n0f/union_myel_var_sketch_pca_harmony_umap_res0.04_project_joined.rds" )

# myeloid
DefaultAssay(obj) = 'RNA'
obj$keep = FALSE
obj$keep[ obj$harmony.cluster.full %in% c(0, 4, 7, 8) ] = TRUE
obj <- subset(obj, subset = keep == TRUE )

# re-generate seurat object and remove extraneous data
options(Seurat.object.assay.version = "v5")
obj <- CreateSeuratObject(counts = obj[["RNA"]]$counts, meta.data = obj@meta.data)
saveRDS( obj, 'n0f/it2_union_myel.rds', compress = F )

# re-perform BPcells
counts.mat = obj[['RNA']]$counts
counts.mat = convert_matrix_type( counts.mat, type="uint32_t" )
counts.mat = as( counts.mat, Class='dgCMatrix' )
write_matrix_dir(mat = counts.mat, dir = "n0f/it2_union_myel_counts", overwrite = T )
counts.mat <- open_matrix_dir(     dir = "n0f/it2_union_myel_counts" )
obj = CreateSeuratObject( counts.mat, meta.data = obj@meta.data )
saveRDS( obj, 'n0f/it2_union_myel_counts.rds', compress = F )

# run downstream analysis
obj <- NormalizeData(obj)
obj[["RNA"]] <- split(obj[["RNA"]], f = obj$study)
obj <- FindVariableFeatures(obj, verbose = T )
obj <- ScaleData(obj, verbose = T)
obj <- RunPCA(obj, npcs = 30, verbose = T)
obj <- IntegrateLayers(obj, method = HarmonyIntegration, orig.reduction = "pca", new.reduction = "harmony" )
obj <- RunUMAP(obj, reduction = "harmony", dims = 1:30, reduction.name = "umap.harmony", return.model = T, min.dist = 0.001)
saveRDS( obj, 'n0f/it2_union_myel_harmony_umap.rds', compress = F ) 

obj <- FindNeighbors(obj, reduction = "harmony", dims = 1:30)
for( res in seq(0.03, 0.25, 0.01) ) {
  obj <- FindClusters(obj, resolution = res )
  saveRDS( obj, paste0( 'n0f/it2_union_myel_harmony_clusters_', res, '.rds' ), compress = F ) 
}

# join layers and find clusters
obj <- readRDS("n0f/it2_union_myel_harmony_clusters_0.06.rds")
DefaultAssay(obj) <- "RNA"
obj <- JoinLayers(obj)
saveRDS(obj, "n0f/it2_union_myel_harmony_clusters_0.06_joined.rds", compress = F)
obj[["RNA"]]$data <- as(object = obj[["RNA"]]$data, Class = "dgCMatrix")
markers <- FindAllMarkers(obj, max.cells.per.ident = 1000, only.pos = T)
write.csv(markers, "n0f/it2_union_myel_harmony_clusters_0.06_markers.csv")
markers$rank <- markers$pct.1 / (markers$pct.2 + 0.001) * markers$avg_log2FC
topmarkers <- markers %>%
  group_by(cluster) %>%
  slice_max(n = 100, order_by = rank)
write.csv(topmarkers, "n0f/it2_union_myel_harmony_clusters_0.06_topmarkers.csv")


### Myeloid Iteration 3 ----

obj <- readRDS( "n0f/it2_union_myel_harmony_clusters_0.06_joined.rds" )

# myeloid
DefaultAssay(obj) = 'RNA'
obj$keep = FALSE
obj$keep[ obj$seurat_clusters %in% c(0:3, 5:7) ] = TRUE
obj <- subset(obj, subset = keep == TRUE )

# re-generate seurat object and remove extraneous data
options(Seurat.object.assay.version = "v5")
obj <- CreateSeuratObject(counts = obj[["RNA"]]$counts, meta.data = obj@meta.data)
saveRDS( obj, 'n0f/it3_union_myel.rds', compress = F )

# re-perform BPcells
counts.mat = obj[['RNA']]$counts
counts.mat = convert_matrix_type( counts.mat, type="uint32_t" )
counts.mat = as( counts.mat, Class='dgCMatrix' )
write_matrix_dir(mat = counts.mat, dir = "n0f/it3_union_myel_counts", overwrite = T )
counts.mat <- open_matrix_dir(     dir = "n0f/it3_union_myel_counts" )
obj = CreateSeuratObject( counts.mat, meta.data = obj@meta.data )
saveRDS( obj, 'n0f/it3_union_myel_counts.rds', compress = F )

# run downstream analysis
obj <- NormalizeData(obj)
obj[["RNA"]] <- split(obj[["RNA"]], f = obj$study)
obj <- FindVariableFeatures(obj, verbose = T )
obj <- ScaleData(obj, verbose = T)
obj <- RunPCA(obj, npcs = 30, verbose = T)
obj <- IntegrateLayers(obj, method = HarmonyIntegration, orig.reduction = "pca", new.reduction = "harmony" )
obj <- RunUMAP(obj, reduction = "harmony", dims = 1:30, reduction.name = "umap.harmony", return.model = T, min.dist = 0.001)
obj <- JoinLayers(obj)
saveRDS( obj, 'n0f/it3_union_myel_harmony_umap_joined.rds', compress = F ) 

obj <- FindNeighbors(obj, reduction = "harmony", dims = 1:30)
for( res in seq( 0.04, 0.3, 0.02) ) { # choose 0.2
  obj <- FindClusters(obj, resolution = res )
  saveRDS( obj, paste0( 'n0f/it3_union_myel_harmony_clusters_', res, '.rds' ), compress = F ) 
}

# join layers and find clusters
obj <- readRDS("n0f/it3_union_myel_harmony_clusters_0.2.rds")
obj[["RNA"]]$data <- as(object = obj[["RNA"]]$data, Class = "dgCMatrix")
markers <- FindAllMarkers(obj, max.cells.per.ident = 1000, only.pos = T)
write.csv(markers, "n0f/it3_union_myel_harmony_clusters_0.2_markers.csv")
markers$rank <- markers$pct.1 / (markers$pct.2 + 0.001) * markers$avg_log2FC
topmarkers <- markers %>%
  group_by(cluster) %>%
  slice_max(n = 100, order_by = rank)
write.csv(topmarkers, "n0f/it3_union_myel_harmony_clusters_0.2_topmarkers.csv")

### Myeloid Iteration 4 ----

obj <- readRDS( "n0f/it3_union_myel_harmony_clusters_0.2.rds" )

# myeloid
DefaultAssay(obj) = 'RNA'
obj$keep = FALSE
obj$keep[ obj$seurat_clusters %in% c( 0:8, 10:13, 18, 20:21 ) ] = TRUE
obj <- subset(obj, subset = keep == TRUE )

# re-generate seurat object and remove extraneous data
options(Seurat.object.assay.version = "v5")
obj <- CreateSeuratObject(counts = obj[["RNA"]]$counts, meta.data = obj@meta.data)
saveRDS( obj, 'n0f/it4_union_myel.rds', compress = F )

# re-perform BPcells
counts.mat = obj[['RNA']]$counts
counts.mat = convert_matrix_type( counts.mat, type="uint32_t" )
counts.mat = as( counts.mat, Class='dgCMatrix' )
write_matrix_dir(mat = counts.mat, dir = "n0f/it4_union_myel_counts", overwrite = T )
counts.mat <- open_matrix_dir(     dir = "n0f/it4_union_myel_counts" )
obj = CreateSeuratObject( counts.mat, meta.data = obj@meta.data )
saveRDS( obj, 'n0f/it4_union_myel_counts.rds', compress = F )

# run downstream analysis
obj <- NormalizeData(obj)
obj[["RNA"]] <- split(obj[["RNA"]], f = obj$study)
obj <- FindVariableFeatures(obj, verbose = T )
obj <- ScaleData(obj, verbose = T)
obj <- RunPCA(obj, npcs = 30, verbose = T)
obj <- IntegrateLayers(obj, method = HarmonyIntegration, orig.reduction = "pca", new.reduction = "harmony" )
obj <- RunUMAP(obj, reduction = "harmony", dims = 1:30, reduction.name = "umap.harmony", return.model = T, min.dist = 0.001)
obj <- JoinLayers(obj)
saveRDS( obj, 'n0f/it4_union_myel_harmony_umap_joined.rds', compress = F ) 

obj <- FindNeighbors(obj, reduction = "harmony", dims = 1:30)
for( res in seq( 0.04, 0.3, 0.02) ) { # choose 0.26
  obj <- FindClusters(obj, resolution = res )
  saveRDS( obj, paste0( 'n0f/it4_union_myel_harmony_clusters_', res, '.rds' ), compress = F ) 
}

# join layers and find clusters
obj <- readRDS("n0f/it4_union_myel_harmony_clusters_0.26.rds")
obj[["RNA"]]$data <- as(object = obj[["RNA"]]$data, Class = "dgCMatrix")
markers <- FindAllMarkers(obj, max.cells.per.ident = 1000, only.pos = T)
write.csv(markers, "n0f/it4_union_myel_harmony_clusters_0.26_markers.csv")
markers$rank <- markers$pct.1 / (markers$pct.2 + 0.001) * markers$avg_log2FC
topmarkers <- markers %>%
  group_by(cluster) %>%
  slice_max(n = 100, order_by = rank)
write.csv(topmarkers, "n0f/it4_union_myel_harmony_clusters_0.26_topmarkers.csv")

### Myeloid Iteration 5 ----

obj <- readRDS("n0f/it4_union_myel_harmony_clusters_0.26.rds")

# myeloid
DefaultAssay(obj) = 'RNA'
obj$fibro = FALSE
obj$fibro[ obj$seurat_clusters %in% c(0:3, 5:8, 10:14, 17:22) ] = TRUE
obj <- subset(obj, subset = fibro == TRUE )

# re-generate seurat object and remove extraneous data
options(Seurat.object.assay.version = "v5")
obj <- CreateSeuratObject(counts = obj[["RNA"]]$counts, meta.data = obj@meta.data)
saveRDS( obj, 'n0f/it5_union_myel.rds', compress = F )

# re-perform BPcells
counts.mat = obj[['RNA']]$counts
counts.mat = convert_matrix_type( counts.mat, type="uint32_t" )
counts.mat = as( counts.mat, Class='dgCMatrix' )
write_matrix_dir(mat = counts.mat, dir = "n0f/it5_union_myel_counts", overwrite = T )
counts.mat <- open_matrix_dir(     dir = "n0f/it5_union_myel_counts" )
obj = CreateSeuratObject( counts.mat, meta.data = obj@meta.data )
saveRDS( obj, 'n0f/it5_union_myel_counts.rds', compress = F )

# run downstream analysis
obj <- NormalizeData(obj)
obj[["RNA"]] <- split(obj[["RNA"]], f = obj$study)
obj <- FindVariableFeatures(obj, verbose = T )
obj <- ScaleData(obj, verbose = T)
obj <- RunPCA(obj, npcs = 30, verbose = T)
obj <- IntegrateLayers(obj, method = HarmonyIntegration, orig.reduction = "pca", new.reduction = "harmony" )
obj <- RunUMAP(obj, reduction = "harmony", dims = 1:30, reduction.name = "umap.harmony", return.model = T, min.dist = 0.001)
obj <- JoinLayers(obj)
saveRDS( obj, 'n0f/it5_union_myel_harmony_umap_joined.rds', compress = F ) 

obj <- FindNeighbors(obj, reduction = "harmony", dims = 1:30)
for( res in seq(0.04, 0.3, 0.02) ) { # choose 0.3
  obj <- FindClusters(obj, resolution = res )
  saveRDS( obj, paste0( 'n0f/it5_union_myel_harmony_clusters_', res, '.rds' ), compress = F ) 
}

# join layers and find clusters
obj <- readRDS("n0f/it5_union_myel_harmony_clusters_0.3.rds")
obj <- subset(obj, subset = seurat_clusters %in% c(0:30))
obj[["RNA"]]$data <- as(object = obj[["RNA"]]$data, Class = "dgCMatrix")
markers <- FindAllMarkers(obj, max.cells.per.ident = 1000, only.pos = T)
write.csv(markers, "n0f/it5_union_myel_harmony_clusters_0.3_markers.csv")
markers$rank <- markers$pct.1 / (markers$pct.2 + 0.001) * markers$avg_log2FC
topmarkers <- markers %>%
  group_by(cluster) %>%
  slice_max(n = 100, order_by = rank)
write.csv(topmarkers, "n0f/it5_union_myel_harmony_clusters_0.3_topmarkers.csv")

# un-BP-cellize for further distribution
obj <- readRDS("n0f/it5_union_myel_harmony_clusters_0.3.rds")
obj[['RNA']]$counts <- as( obj[['RNA']]$counts, Class='dgCMatrix' )
obj[['RNA']]$data <- as( obj[['RNA']]$data, Class='dgCMatrix' )
obj[['RNA']]$scale.data <- as( obj[['RNA']]$scale.data, Class='dgCMatrix' )
saveRDS(obj, "n0f/it5_union_myel_harmony_clusters_0.3_sparse.rds", compress = F)

# manual annotations for entire dataset
obj <- subset(obj, subset = seurat_clusters %in% (0:12))
obj$manual.annotation <- "tbd"
obj$manual.annotation[obj$seurat_clusters == 0] <- "Macrophage C1QC+"
obj$manual.annotation[obj$seurat_clusters == 1] <- "DC cDC2 CD1C+"
obj$manual.annotation[obj$seurat_clusters == 2] <- "Monocyte CD14+"
obj$manual.annotation[obj$seurat_clusters == 3] <- "Neutrophil CSF3R+"
obj$manual.annotation[obj$seurat_clusters == 4] <- "Mast CPA3+"
obj$manual.annotation[obj$seurat_clusters == 5] <- "Macrophage SPP1+"
obj$manual.annotation[obj$seurat_clusters == 6] <- "Macrophage ISG15+"
obj$manual.annotation[obj$seurat_clusters == 7] <- "Macrophage MARCO+"
obj$manual.annotation[obj$seurat_clusters == 8] <- "Monocyte CD16+"
obj$manual.annotation[obj$seurat_clusters == 9] <- "DC pDC TCF4+"
obj$manual.annotation[obj$seurat_clusters == 10] <- "DC mregDC LAMP3+"
obj$manual.annotation[obj$seurat_clusters == 11] <- "DC cDC1 CLEC9A+"
obj$manual.annotation[obj$seurat_clusters == 12] <- "Osteoclast CTSK+"

saveRDS(obj, "n0f/final_figures/human_myel.rds", compress = F)

### Nerve Iteration 1 (no neural crest tumors) ----
# remove melanoma, MPNST, NB, RB, teratoma, GIST because challenging to differentiate tumor vs non tumor 

obj     = readRDS( "/oak/stanford/groups/longaker/CAF_META/newHuman/n0d/humanMerged_norm_findvar_sketch1000_harmony_umap_res0.06_project_joined.rds" )
obj_har = readRDS( "/oak/stanford/groups/longaker/CAF_META/newHuman/n0d/humanMerged_norm_findvar_sketch1000_harmony_umap_res0.06_clusterList.rds" )
obj_scv = readRDS( "/oak/stanford/groups/longaker/CAF_META/newHuman/n0d/humanMerged_norm_findvar_sketch1000_scvi_umap_res0.09_clusterList.rds"    )
obj_mnn = readRDS( "/oak/stanford/groups/longaker/CAF_META/newHuman/n0d/humanMerged_norm_findvar_sketch1000_mnn_umap_res0.08_clusterList.rds"    )

# Nerve
ii = ( colnames(obj) %in% colnames(obj)[ (obj_har %in% c(8)) ] ) |
  ( colnames(obj) %in% colnames(obj)[ (obj_scv %in% c(17)) ] ) |
  ( colnames(obj) %in% colnames(obj)[ (obj_mnn %in% c(8,10)) ] )

DefaultAssay(obj) = 'RNA'
obj$keep = FALSE
obj$keep[ ii ] = TRUE
obj$keep[obj$cancer %in% c("cMelanoma", "GIST", "MPNST", "NB", "RB", "Teratoma", "uMelanoma")] <- F
DimPlot(obj, group.by = "keep", raster = T, reduction = "full.umap.harmony", alpha = 0.05)
obj <- subset(obj, subset = fibro == TRUE )

# re-generate seurat object and remove extraneous data
options(Seurat.object.assay.version = "v5")
obj <- CreateSeuratObject(counts = obj[["RNA"]]$counts, meta.data = obj@meta.data)
saveRDS( obj, 'n0h/union_nerve.rds', compress = F )

# re-perform BPcells
counts.mat = obj[['RNA']]$counts
counts.mat = convert_matrix_type( counts.mat, type="uint32_t" )
counts.mat = as( counts.mat, Class='dgCMatrix' )
write_matrix_dir(mat = counts.mat, dir = "n0h/union_nerve_counts", overwrite = T )
counts.mat <- open_matrix_dir(     dir = "n0h/union_nerve_counts" )
obj = CreateSeuratObject( counts.mat, meta.data = obj@meta.data )
saveRDS( obj, 'n0h/union_nerve_counts.rds', compress = F )

# run downstream analysis
obj <- NormalizeData(obj)
obj[["RNA"]] <- split(obj[["RNA"]], f = obj$study)
obj <- FindVariableFeatures(obj, verbose = T )
obj <- ScaleData(obj, verbose = T)
obj <- RunPCA(obj, npcs = 30, verbose = T)
obj <- IntegrateLayers(obj, method = HarmonyIntegration, orig.reduction = "pca", new.reduction = "harmony" )
obj <- RunUMAP(obj, reduction = "harmony", dims = 1:30, reduction.name = "umap.harmony", return.model = T, min.dist = 0.001)
obj <- JoinLayers(obj)
saveRDS( obj, 'n0h/union_nerve_harmony_umap_joined.rds', compress = F ) 

obj <- FindNeighbors(obj, reduction = "harmony", dims = 1:30)
for( res in seq(0.03, 0.1, 0.01) ) { # choose 0.03
  obj <- FindClusters(obj, resolution = res )
  saveRDS( obj, paste0( 'n0h/union_nerve_harmony_clusters_', res, '.rds' ), compress = F ) 
}

# join layers and find clusters
obj <- readRDS("n0h/union_nerve_harmony_clusters_0.03.rds")
obj[["RNA"]]$data <- as(object = obj[["RNA"]]$data, Class = "dgCMatrix")
markers <- FindAllMarkers(obj, max.cells.per.ident = 1000, only.pos = T)
write.csv(markers, "n0h/union_nerve_harmony_clusters_0.03_markers.csv")
markers$rank <- markers$pct.1 / (markers$pct.2 + 0.001) * markers$avg_log2FC
topmarkers <- markers %>%
  group_by(cluster) %>%
  slice_max(n = 100, order_by = rank)
write.csv(topmarkers, "n0h/union_nerve_harmony_clusters_0.03_topmarkers.csv")


### Nerve Iteration 2 (no neural crest tumors) ----

obj <- readRDS( "n0h/union_nerve_harmony_clusters_0.03.rds" )

# nerve
DefaultAssay(obj) = 'RNA'
obj$keep = FALSE
obj$keep[ obj$seurat_clusters %in% c(1,2) ] = TRUE
obj <- subset(obj, subset = keep == TRUE )

# re-generate seurat object and remove extraneous data
options(Seurat.object.assay.version = "v5")
obj <- CreateSeuratObject(counts = obj[["RNA"]]$counts, meta.data = obj@meta.data)
saveRDS( obj, 'n0h/it2_union_nerve.rds', compress = F )

# re-perform BPcells
counts.mat = obj[['RNA']]$counts
counts.mat = convert_matrix_type( counts.mat, type="uint32_t" )
counts.mat = as( counts.mat, Class='dgCMatrix' )
write_matrix_dir(mat = counts.mat, dir = "n0h/it2_union_nerve_counts", overwrite = T )
counts.mat <- open_matrix_dir(     dir = "n0h/it2_union_nerve_counts" )
obj = CreateSeuratObject( counts.mat, meta.data = obj@meta.data )
saveRDS( obj, 'n0h/it2_union_nerve_counts.rds', compress = F )

# run downstream analysis
obj <- NormalizeData(obj)
obj[["RNA"]] <- split(obj[["RNA"]], f = obj$study)
obj <- FindVariableFeatures(obj, verbose = T )
obj <- ScaleData(obj, verbose = T)
obj <- RunPCA(obj, npcs = 30, verbose = T)
obj <- IntegrateLayers(obj, method = HarmonyIntegration, orig.reduction = "pca", new.reduction = "harmony" )
obj <- RunUMAP(obj, reduction = "harmony", dims = 1:30, reduction.name = "umap.harmony", return.model = T, min.dist = 0.001)
obj <- JoinLayers(obj)
saveRDS( obj, 'n0h/it2_union_nerve_harmony_umap_joined.rds', compress = F ) 

obj <- FindNeighbors(obj, reduction = "harmony", dims = 1:30)
for( res in seq( 0.04, 0.09, 0.01 )) { # choose 0.08
  obj <- FindClusters(obj, resolution = res )
  saveRDS( obj, paste0( 'n0h/it2_union_nerve_harmony_clusters_', res, '.rds' ), compress = F ) 
}

# join layers and find clusters
obj <- readRDS("n0h/it2_union_nerve_harmony_clusters_0.08.rds")
DefaultAssay(obj) <- "RNA"
obj[["RNA"]]$data <- as(object = obj[["RNA"]]$data, Class = "dgCMatrix")
markers <- FindAllMarkers(obj, max.cells.per.ident = 1000, only.pos = T)
write.csv(markers, "n0h/it2_union_nerve_harmony_clusters_0.08_markers.csv")
markers$rank <- markers$pct.1 / (markers$pct.2 + 0.001) * markers$avg_log2FC
topmarkers <- markers %>%
  group_by(cluster) %>%
  slice_max(n = 100, order_by = rank)
write.csv(topmarkers, "n0h/it2_union_nerve_harmony_clusters_0.08_topmarkers.csv")

### Nerve Iteration 3 (no neural crest tumors) ----

obj <- readRDS( "n0h/it2_union_nerve_harmony_clusters_0.08.rds" )

# nerve
DefaultAssay(obj) = 'RNA'
obj$keep = FALSE
obj$keep[ obj$seurat_clusters %in% c(4) ] = TRUE
obj <- subset(obj, subset = keep == TRUE )

# re-generate seurat object and remove extraneous data
options(Seurat.object.assay.version = "v5")
obj <- CreateSeuratObject(counts = obj[["RNA"]]$counts, meta.data = obj@meta.data)
saveRDS( obj, 'n0h/it3_union_nerve.rds', compress = F )

# re-perform BPcells
counts.mat = obj[['RNA']]$counts
counts.mat = convert_matrix_type( counts.mat, type="uint32_t" )
counts.mat = as( counts.mat, Class='dgCMatrix' )
write_matrix_dir(mat = counts.mat, dir = "n0h/it3_union_nerve_counts", overwrite = T )
counts.mat <- open_matrix_dir(     dir = "n0h/it3_union_nerve_counts" )
obj = CreateSeuratObject( counts.mat, meta.data = obj@meta.data )
saveRDS( obj, 'n0h/it3_union_nerve_counts.rds', compress = F )

# run downstream analysis
obj <- NormalizeData(obj)
obj[["RNA"]] <- split(obj[["RNA"]], f = obj$study)
obj <- FindVariableFeatures(obj, verbose = T )
obj <- ScaleData(obj, verbose = T)
obj <- RunPCA(obj, npcs = 30, verbose = T)
obj <- IntegrateLayers(obj, method = HarmonyIntegration, orig.reduction = "pca", new.reduction = "harmony" )
obj <- RunUMAP(obj, reduction = "harmony", dims = 1:30, reduction.name = "umap.harmony", return.model = T, min.dist = 0.001)
saveRDS( obj, 'n0h/it3_union_nerve_harmony_umap.rds', compress = F ) 

obj <- FindNeighbors(obj, reduction = "harmony", dims = 1:30)
for( res in seq( 0.03, 0.2, 0.01)) { # choose 0.03
  obj <- FindClusters(obj, resolution = res )
  saveRDS( obj, paste0( 'n0h/it3_union_nerve_harmony_clusters_', res, '.rds' ), compress = F ) 
}

# join layers and find clusters
obj <- readRDS("n0h/it3_union_nerve_harmony_clusters_0.03.rds")
DefaultAssay(obj) <- "RNA"
obj <- JoinLayers(obj)
saveRDS(obj, "n0h/it3_union_nerve_harmony_clusters_0.03_joined.rds", compress = F)
obj[["RNA"]]$data <- as(object = obj[["RNA"]]$data, Class = "dgCMatrix")
markers <- FindAllMarkers(obj, max.cells.per.ident = 1000, only.pos = T)
write.csv(markers, "n0h/it3_union_nerve_harmony_clusters_0.03_markers.csv")
markers$rank <- markers$pct.1 / (markers$pct.2 + 0.001) * markers$avg_log2FC
topmarkers <- markers %>%
  group_by(cluster) %>%
  slice_max(n = 100, order_by = rank)
write.csv(topmarkers, "n0h/it3_union_nerve_harmony_clusters_0.03_topmarkers.csv")


### Nerve Iteration 4 (no neural crest tumors) ----

obj <- readRDS( "n0h/it3_union_nerve_harmony_clusters_0.03_joined.rds" )

# nerve
DefaultAssay(obj) = 'RNA'
obj$keep = FALSE
obj$keep[ obj$seurat_clusters %in% c(0,2) ] = TRUE
obj <- subset(obj, subset = keep == TRUE )

# re-generate seurat object and remove extraneous data
options(Seurat.object.assay.version = "v5")
obj <- CreateSeuratObject(counts = obj[["RNA"]]$counts, meta.data = obj@meta.data)
saveRDS( obj, 'n0h/it4_union_nerve.rds', compress = F )

# re-perform BPcells
counts.mat = obj[['RNA']]$counts
counts.mat = convert_matrix_type( counts.mat, type="uint32_t" )
counts.mat = as( counts.mat, Class='dgCMatrix' )
write_matrix_dir(mat = counts.mat, dir = "n0h/it4_union_nerve_counts", overwrite = T )
counts.mat <- open_matrix_dir(     dir = "n0h/it4_union_nerve_counts" )
obj = CreateSeuratObject( counts.mat, meta.data = obj@meta.data )
saveRDS( obj, 'n0h/it4_union_nerve_counts.rds', compress = F )

# run downstream analysis
obj <- NormalizeData(obj)
obj[["RNA"]] <- split(obj[["RNA"]], f = obj$study)
obj <- FindVariableFeatures(obj, verbose = T )
obj <- ScaleData(obj, verbose = T)
obj <- RunPCA(obj, npcs = 30, verbose = T)
obj <- IntegrateLayers(obj, method = HarmonyIntegration, orig.reduction = "pca", new.reduction = "harmony" )
obj <- RunUMAP(obj, reduction = "harmony", dims = 1:30, reduction.name = "umap.harmony", return.model = T, min.dist = 0.001)
saveRDS( obj, 'n0h/it4_union_nerve_harmony_umap.rds', compress = F ) 

obj <- FindNeighbors(obj, reduction = "harmony", dims = 1:30)
for( res in seq( 0.03, 0.2, 0.01) ) { # choose 0.03
  obj <- FindClusters(obj, resolution = res )
  saveRDS( obj, paste0( 'n0h/it4_union_nerve_harmony_clusters_', res, '.rds' ), compress = F ) 
}

# join layers and find clusters
obj <- readRDS("n0h/it4_union_nerve_harmony_clusters_0.03.rds")
DefaultAssay(obj) <- "RNA"
obj <- JoinLayers(obj)
saveRDS(obj, "n0h/it4_union_nerve_harmony_clusters_0.03_joined.rds", compress = F)
obj[["RNA"]]$data <- as(object = obj[["RNA"]]$data, Class = "dgCMatrix")
markers <- FindAllMarkers(obj, max.cells.per.ident = 1000, only.pos = T)
write.csv(markers, "n0h/it4_union_nerve_harmony_clusters_0.03_markers.csv")
markers$rank <- markers$pct.1 / (markers$pct.2 + 0.001) * markers$avg_log2FC
topmarkers <- markers %>%
  group_by(cluster) %>%
  slice_max(n = 100, order_by = rank)
write.csv(topmarkers, "n0h/it4_union_nerve_harmony_clusters_0.03_topmarkers.csv")

### Nerve Iteration 5 (no neural crest tumors) ----

obj <- readRDS( "n0h/it4_union_nerve_harmony_clusters_0.03_joined.rds" )

# nerve
DefaultAssay(obj) = 'RNA'
obj$keep = FALSE
obj$keep[ obj$seurat_clusters %in% c(1,2) ] = TRUE
obj <- subset(obj, subset = keep == TRUE )

# re-generate seurat object and remove extraneous data
options(Seurat.object.assay.version = "v5")
obj <- CreateSeuratObject(counts = obj[["RNA"]]$counts, meta.data = obj@meta.data)
saveRDS( obj, 'n0h/it5_union_nerve.rds', compress = F )

# re-perform BPcells
counts.mat = obj[['RNA']]$counts
counts.mat = convert_matrix_type( counts.mat, type="uint32_t" )
counts.mat = as( counts.mat, Class='dgCMatrix' )
write_matrix_dir(mat = counts.mat, dir = "n0h/it5_union_nerve_counts", overwrite = T )
counts.mat <- open_matrix_dir(     dir = "n0h/it5_union_nerve_counts" )
obj = CreateSeuratObject( counts.mat, meta.data = obj@meta.data )
saveRDS( obj, 'n0h/it5_union_nerve_counts.rds', compress = F )

# run downstream analysis
obj <- NormalizeData(obj)
obj[["RNA"]] <- split(obj[["RNA"]], f = obj$study)
obj <- FindVariableFeatures(obj, verbose = T )
obj <- ScaleData(obj, verbose = T)
obj <- RunPCA(obj, npcs = 30, verbose = T)
obj <- IntegrateLayers(obj, method = HarmonyIntegration, orig.reduction = "pca", new.reduction = "harmony" )
obj <- RunUMAP(obj, reduction = "harmony", dims = 1:30, reduction.name = "umap.harmony", return.model = T, min.dist = 0.001)
saveRDS( obj, 'n0h/it5_union_nerve_harmony_umap.rds', compress = F ) 

obj <- FindNeighbors(obj, reduction = "harmony", dims = 1:30)
for( res in seq( 0.03, 0.1, 0.01) ) { # choose 0.06
  obj <- FindClusters(obj, resolution = res )
  saveRDS( obj, paste0( 'n0h/it5_union_nerve_harmony_clusters_', res, '.rds' ), compress = F ) 
}

# join layers and find clusters
obj <- readRDS("n0h/it5_union_nerve_harmony_clusters_0.06.rds")
DefaultAssay(obj) <- "RNA"
obj <- JoinLayers(obj)
saveRDS(obj, "n0h/it5_union_nerve_harmony_clusters_0.06_joined.rds", compress = F)
obj[["RNA"]]$data <- as(object = obj[["RNA"]]$data, Class = "dgCMatrix")
markers <- FindAllMarkers(obj, max.cells.per.ident = 1000, only.pos = T)
write.csv(markers, "n0h/it5_union_nerve_harmony_clusters_0.06_markers.csv")
markers$rank <- markers$pct.1 / (markers$pct.2 + 0.001) * markers$avg_log2FC
topmarkers <- markers %>%
  group_by(cluster) %>%
  slice_max(n = 100, order_by = rank)
write.csv(topmarkers, "n0h/it5_union_nerve_harmony_clusters_0.06_topmarkers.csv")

# un-BP-cellize for further distribution
obj <- readRDS("n0h/it5_union_nerve_harmony_clusters_0.06_joined.rds")
obj[['RNA']]$counts <- as( obj[['RNA']]$counts, Class='dgCMatrix' )
obj[['RNA']]$data <- as( obj[['RNA']]$data, Class='dgCMatrix' )
obj[['RNA']]$scale.data <- as( obj[['RNA']]$scale.data, Class='dgCMatrix' )
saveRDS(obj, "n0h/it5_union_nerve_harmony_clusters_0.06_joined_sparse.rds", compress = F)

# manual annotations for entire dataset
obj$manual.annotation <- "tbd"
obj$manual.annotation[obj$seurat_clusters == 0] <- "Schwann myel MPZ+"
obj$manual.annotation[obj$seurat_clusters == 1] <- "Melanocyte MLANA+"
obj$manual.annotation[obj$seurat_clusters == 2] <- "Schwann non-myel NGFR+"

saveRDS(obj, "n0h/final_figures/human_nerve.rds", compress = F)

### Human-Mouse Integration Fibroblast ----
obj_mouse <- readRDS("../newMouse/sketch_v5_all/final_figures/mouse_union_fibro_it4_res0.26_joined_annotated_sparse.rds")
obj_human <- readRDS("n0h/final_figures/human_fibro4.rds")

library(nichenetr)

mouseToHumanNames = convert_mouse_to_human_symbols( rownames( obj_mouse ) )
ii = !(is.na(mouseToHumanNames))
counts.mat = GetAssayData( obj_mouse, slot = 'counts' )
counts.mat = counts.mat[ ii, ]
rownames(counts.mat) = mouseToHumanNames[ii]
counts.mat = as(counts.mat, Class = "dgCMatrix")
counts.mat = counts.mat[ !duplicated(rownames(counts.mat)), ]
options(Seurat.object.assay.version = "v5")
obj_mouse = CreateSeuratObject( counts.mat, meta.data = obj_mouse@meta.data )
obj_mouse[["RNA"]]$counts = as(obj_mouse[["RNA"]]$counts, Class = "dgCMatrix")
obj_mouse@meta.data <- obj_mouse@meta.data %>%
  select(any_of(c("orig.ident", "nCount_RNA","nFeature_RNA","study",
                  "organ","cancer","model","sorting","site","percent.mt",
                  "manual.annotation")))
write_matrix_dir(mat = obj_mouse[["RNA"]]$counts, dir = "species/mouse_toAnchor_fibro_counts", overwrite=T)
obj_mouse[["RNA"]]$counts <- open_matrix_dir(dir = "species/mouse_toAnchor_fibro_counts")
saveRDS(obj_mouse, "species/mouse_toAnchor_fibro_counts.rds",compress = F)

obj_human = CreateSeuratObject( obj_human[["RNA"]]$counts, meta.data = obj_human@meta.data )
obj_human[["RNA"]]$counts = as(obj_human[["RNA"]]$counts, Class = "dgCMatrix")
obj_human@meta.data <- obj_human@meta.data %>%
  select(any_of(c("orig.ident", "nCount_RNA","nFeature_RNA","study",
                  "organ","cancer","sorting","site","percent.mt",
                  "manual.annotation")))
write_matrix_dir(mat = obj_human[["RNA"]]$counts, dir = "species/human_toAnchor_fibro_counts", overwrite=T)
obj_human[["RNA"]]$counts <- open_matrix_dir(dir = "species/human_toAnchor_fibro_counts")
saveRDS(obj_human, "species/human_toAnchor_fibro_counts.rds",compress = F)

rm(list = ls())
obj_mouse = readRDS("species/mouse_toAnchor_fibro_counts.rds")
obj_human = readRDS("species/human_toAnchor_fibro_counts.rds")
data_mouse = GetAssayData(obj_mouse, slot = 'counts' )
data_human = GetAssayData(obj_human, slot = 'counts' )
commonGenes = intersect( rownames(data_mouse), rownames(data_human) )
data_mouse = data_mouse[ rownames(data_mouse) %in% commonGenes, ]
data_human = data_human[ rownames(data_human) %in% commonGenes, ]
obj_mouse = CreateSeuratObject( data_mouse, meta.data = obj_mouse@meta.data )
obj_human = CreateSeuratObject( data_human, meta.data = obj_human@meta.data )
saveRDS( obj_mouse, "species/mouse_toAnchor_fibro_counts_commonGenes.rds", compress = F )
saveRDS( obj_human, "species/human_toAnchor_fibro_counts_commonGenes.rds", compress = F )

obj_mouse$species = 'mouse'
obj_human$species = 'human'

obj_mouse$species.study <-  paste(obj_mouse$species, obj_mouse$study, sep = "_")
obj_human$species.study <-  paste(obj_human$species, obj_human$study, sep = "_")

obj = merge( obj_mouse, obj_human )
obj = JoinLayers( obj )
saveRDS( obj, "species/human_mouse_fibro_merged.rds", compress = F )

# run downstream analysis
obj <- NormalizeData(obj)
obj[["RNA"]] <- split(obj[["RNA"]], f = obj$species.study)
obj <- FindVariableFeatures(obj, verbose = T )
obj <- ScaleData(obj, verbose = T)
obj <- RunPCA(obj, npcs = 30, verbose = T)
obj <- IntegrateLayers(obj, method = HarmonyIntegration, orig.reduction = "pca", new.reduction = "harmony" )
obj <- RunUMAP(obj, reduction = "harmony", dims = 1:30, reduction.name = "umap.harmony", return.model = T, min.dist = 0.001)
obj <- JoinLayers(obj)
saveRDS( obj, 'species/human_mouse_fibro_harmony_umap_joined.rds', compress = F ) 

### Human-Mouse Integration Endothelial ----
obj_mouse <- readRDS("../newMouse/sketch_v5_all/final_figures/mouse_union_endo_it3_res0.22_joined_annotated_sparse.rds")
obj_human <- readRDS("n0f/final_figures/human_endo.rds")

library(nichenetr)

mouseToHumanNames = convert_mouse_to_human_symbols( rownames( obj_mouse ) )
ii = !(is.na(mouseToHumanNames))
counts.mat = GetAssayData( obj_mouse, slot = 'counts' )
counts.mat = counts.mat[ ii, ]
rownames(counts.mat) = mouseToHumanNames[ii]
counts.mat = as(counts.mat, Class = "dgCMatrix")
counts.mat = counts.mat[ !duplicated(rownames(counts.mat)), ]
options(Seurat.object.assay.version = "v5")
obj_mouse = CreateSeuratObject( counts.mat, meta.data = obj_mouse@meta.data )
obj_mouse[["RNA"]]$counts = as(obj_mouse[["RNA"]]$counts, Class = "dgCMatrix")
obj_mouse@meta.data <- obj_mouse@meta.data %>%
  select(any_of(c("orig.ident", "nCount_RNA","nFeature_RNA","study",
                  "organ","cancer","model","sorting","site","percent.mt",
                  "manual.annotation")))
write_matrix_dir(mat = obj_mouse[["RNA"]]$counts, dir = "species/mouse_toAnchor_endo_counts", overwrite=T)
obj_mouse[["RNA"]]$counts <- open_matrix_dir(dir = "species/mouse_toAnchor_endo_counts")
saveRDS(obj_mouse, "species/mouse_toAnchor_endo_counts.rds",compress = F)

obj_human = CreateSeuratObject( obj_human[["RNA"]]$counts, meta.data = obj_human@meta.data )
obj_human[["RNA"]]$counts = as(obj_human[["RNA"]]$counts, Class = "dgCMatrix")
obj_human@meta.data <- obj_human@meta.data %>%
  select(any_of(c("orig.ident", "nCount_RNA","nFeature_RNA","study",
                  "organ","cancer","sorting","site","percent.mt",
                  "manual.annotation")))
write_matrix_dir(mat = obj_human[["RNA"]]$counts, dir = "species/human_toAnchor_endo_counts", overwrite=T)
obj_human[["RNA"]]$counts <- open_matrix_dir(dir = "species/human_toAnchor_endo_counts")
saveRDS(obj_human, "species/human_toAnchor_endo_counts.rds",compress = F)

rm(list = ls())
obj_mouse = readRDS("species/mouse_toAnchor_endo_counts.rds")
obj_human = readRDS("species/human_toAnchor_endo_counts.rds")
data_mouse = GetAssayData(obj_mouse, slot = 'counts' )
data_human = GetAssayData(obj_human, slot = 'counts' )
commonGenes = intersect( rownames(data_mouse), rownames(data_human) )
data_mouse = data_mouse[ rownames(data_mouse) %in% commonGenes, ]
data_human = data_human[ rownames(data_human) %in% commonGenes, ]
obj_mouse = CreateSeuratObject( data_mouse, meta.data = obj_mouse@meta.data )
obj_human = CreateSeuratObject( data_human, meta.data = obj_human@meta.data )
saveRDS( obj_mouse, "species/mouse_toAnchor_endo_counts_commonGenes.rds", compress = F )
saveRDS( obj_human, "species/human_toAnchor_endo_counts_commonGenes.rds", compress = F )

obj_mouse$species = 'mouse'
obj_human$species = 'human'

obj_mouse$species.study <-  paste(obj_mouse$species, obj_mouse$study, sep = "_")
obj_human$species.study <-  paste(obj_human$species, obj_human$study, sep = "_")

obj = merge( obj_mouse, obj_human )
obj = JoinLayers( obj )
saveRDS( obj, "species/human_mouse_endo_merged.rds", compress = F )

# run downstream analysis
obj <- NormalizeData(obj)
obj[["RNA"]] <- split(obj[["RNA"]], f = obj$species.study)
obj <- FindVariableFeatures(obj, verbose = T )
obj <- ScaleData(obj, verbose = T)
obj <- RunPCA(obj, npcs = 30, verbose = T)
obj <- IntegrateLayers(obj, method = HarmonyIntegration, orig.reduction = "pca", new.reduction = "harmony" )
obj <- RunUMAP(obj, reduction = "harmony", dims = 1:30, reduction.name = "umap.harmony", return.model = T, min.dist = 0.001)
obj <- JoinLayers(obj)
saveRDS( obj, 'species/human_mouse_endo_harmony_umap_joined.rds', compress = F ) 

### Human-Mouse Integration Mural ----
obj_mouse <- readRDS("../newMouse/sketch_v5_all/final_figures/mouse_union_mural_it3_res0.05_joined_annotated_sparse.rds")
obj_human <- readRDS("n0f/final_figures/human_mural.rds")

library(nichenetr)

mouseToHumanNames = convert_mouse_to_human_symbols( rownames( obj_mouse ) )
ii = !(is.na(mouseToHumanNames))
counts.mat = GetAssayData( obj_mouse, slot = 'counts' )
counts.mat = counts.mat[ ii, ]
rownames(counts.mat) = mouseToHumanNames[ii]
counts.mat = as(counts.mat, Class = "dgCMatrix")
counts.mat = counts.mat[ !duplicated(rownames(counts.mat)), ]
options(Seurat.object.assay.version = "v5")
obj_mouse = CreateSeuratObject( counts.mat, meta.data = obj_mouse@meta.data )
obj_mouse[["RNA"]]$counts = as(obj_mouse[["RNA"]]$counts, Class = "dgCMatrix")
obj_mouse@meta.data <- obj_mouse@meta.data %>%
  select(any_of(c("orig.ident", "nCount_RNA","nFeature_RNA","study",
                  "organ","cancer","model","sorting","site","percent.mt",
                  "manual.annotation")))
write_matrix_dir(mat = obj_mouse[["RNA"]]$counts, dir = "species/mouse_toAnchor_mural_counts", overwrite=T)
obj_mouse[["RNA"]]$counts <- open_matrix_dir(dir = "species/mouse_toAnchor_mural_counts")
saveRDS(obj_mouse, "species/mouse_toAnchor_mural_counts.rds",compress = F)

obj_human = CreateSeuratObject( obj_human[["RNA"]]$counts, meta.data = obj_human@meta.data )
obj_human[["RNA"]]$counts = as(obj_human[["RNA"]]$counts, Class = "dgCMatrix")
obj_human@meta.data <- obj_human@meta.data %>%
  select(any_of(c("orig.ident", "nCount_RNA","nFeature_RNA","study",
                  "organ","cancer","sorting","site","percent.mt",
                  "manual.annotation")))
write_matrix_dir(mat = obj_human[["RNA"]]$counts, dir = "species/human_toAnchor_mural_counts", overwrite=T)
obj_human[["RNA"]]$counts <- open_matrix_dir(dir = "species/human_toAnchor_mural_counts")
saveRDS(obj_human, "species/human_toAnchor_mural_counts.rds",compress = F)

rm(list = ls())
obj_mouse = readRDS("species/mouse_toAnchor_mural_counts.rds")
obj_human = readRDS("species/human_toAnchor_mural_counts.rds")
data_mouse = GetAssayData(obj_mouse, slot = 'counts' )
data_human = GetAssayData(obj_human, slot = 'counts' )
commonGenes = intersect( rownames(data_mouse), rownames(data_human) )
data_mouse = data_mouse[ rownames(data_mouse) %in% commonGenes, ]
data_human = data_human[ rownames(data_human) %in% commonGenes, ]
obj_mouse = CreateSeuratObject( data_mouse, meta.data = obj_mouse@meta.data )
obj_human = CreateSeuratObject( data_human, meta.data = obj_human@meta.data )
saveRDS( obj_mouse, "species/mouse_toAnchor_mural_counts_commonGenes.rds", compress = F )
saveRDS( obj_human, "species/human_toAnchor_mural_counts_commonGenes.rds", compress = F )

obj_mouse$species = 'mouse'
obj_human$species = 'human'

obj_mouse$species.study <-  paste(obj_mouse$species, obj_mouse$study, sep = "_")
obj_human$species.study <-  paste(obj_human$species, obj_human$study, sep = "_")

obj = merge( obj_mouse, obj_human )
obj = JoinLayers( obj )
saveRDS( obj, "species/human_mouse_mural_merged.rds", compress = F )

# run downstream analysis
obj <- NormalizeData(obj)
obj[["RNA"]] <- split(obj[["RNA"]], f = obj$species.study)
obj <- FindVariableFeatures(obj, verbose = T )
obj <- ScaleData(obj, verbose = T)
obj <- RunPCA(obj, npcs = 30, verbose = T)
obj <- IntegrateLayers(obj, method = HarmonyIntegration, orig.reduction = "pca", new.reduction = "harmony" )
obj <- RunUMAP(obj, reduction = "harmony", dims = 1:30, reduction.name = "umap.harmony", return.model = T, min.dist = 0.001)
obj <- JoinLayers(obj)
saveRDS( obj, 'species/human_mouse_mural_harmony_umap_joined.rds', compress = F ) 

### Human-Mouse Integration Myeloid ----
obj_mouse <- readRDS("../newMouse/sketch_v5_all/final_figures/mouse_union_myel_it4_res0.6_joined_annotated_sparse.rds")
obj_human <- readRDS("n0f/final_figures/human_myel.rds")

library(nichenetr)

mouseToHumanNames = convert_mouse_to_human_symbols( rownames( obj_mouse ) )
ii = !(is.na(mouseToHumanNames))
counts.mat = GetAssayData( obj_mouse, slot = 'counts' )
counts.mat = counts.mat[ ii, ]
rownames(counts.mat) = mouseToHumanNames[ii]
counts.mat = as(counts.mat, Class = "dgCMatrix")
counts.mat = counts.mat[ !duplicated(rownames(counts.mat)), ]
options(Seurat.object.assay.version = "v5")
obj_mouse = CreateSeuratObject( counts.mat, meta.data = obj_mouse@meta.data )
obj_mouse[["RNA"]]$counts = as(obj_mouse[["RNA"]]$counts, Class = "dgCMatrix")
obj_mouse@meta.data <- obj_mouse@meta.data %>%
  select(any_of(c("orig.ident", "nCount_RNA","nFeature_RNA","study",
                  "organ","cancer","model","sorting","site","percent.mt",
                  "manual.annotation")))
write_matrix_dir(mat = obj_mouse[["RNA"]]$counts, dir = "species/mouse_toAnchor_myel_counts", overwrite=T)
obj_mouse[["RNA"]]$counts <- open_matrix_dir(dir = "species/mouse_toAnchor_myel_counts")
saveRDS(obj_mouse, "species/mouse_toAnchor_myel_counts.rds",compress = F)

obj_human = CreateSeuratObject( obj_human[["RNA"]]$counts, meta.data = obj_human@meta.data )
obj_human[["RNA"]]$counts = as(obj_human[["RNA"]]$counts, Class = "dgCMatrix")
obj_human@meta.data <- obj_human@meta.data %>%
  select(any_of(c("orig.ident", "nCount_RNA","nFeature_RNA","study",
                  "organ","cancer","sorting","site","percent.mt",
                  "manual.annotation")))
write_matrix_dir(mat = obj_human[["RNA"]]$counts, dir = "species/human_toAnchor_myel_counts", overwrite=T)
obj_human[["RNA"]]$counts <- open_matrix_dir(dir = "species/human_toAnchor_myel_counts")
saveRDS(obj_human, "species/human_toAnchor_myel_counts.rds",compress = F)

rm(list = ls())
obj_mouse = readRDS("species/mouse_toAnchor_myel_counts.rds")
obj_human = readRDS("species/human_toAnchor_myel_counts.rds")
data_mouse = GetAssayData(obj_mouse, slot = 'counts' )
data_human = GetAssayData(obj_human, slot = 'counts' )
commonGenes = intersect( rownames(data_mouse), rownames(data_human) )
data_mouse = data_mouse[ rownames(data_mouse) %in% commonGenes, ]
data_human = data_human[ rownames(data_human) %in% commonGenes, ]
obj_mouse = CreateSeuratObject( data_mouse, meta.data = obj_mouse@meta.data )
obj_human = CreateSeuratObject( data_human, meta.data = obj_human@meta.data )
saveRDS( obj_mouse, "species/mouse_toAnchor_myel_counts_commonGenes.rds", compress = F )
saveRDS( obj_human, "species/human_toAnchor_myel_counts_commonGenes.rds", compress = F )

obj_mouse$species = 'mouse'
obj_human$species = 'human'

obj_mouse$species.study <-  paste(obj_mouse$species, obj_mouse$study, sep = "_")
obj_human$species.study <-  paste(obj_human$species, obj_human$study, sep = "_")

obj = merge( obj_mouse, obj_human )
obj = JoinLayers( obj )
saveRDS( obj, "species/human_mouse_myel_merged.rds", compress = F )

# run downstream analysis
obj <- NormalizeData(obj)
obj[["RNA"]] <- split(obj[["RNA"]], f = obj$species.study)
obj <- FindVariableFeatures(obj, verbose = T )
obj <- ScaleData(obj, verbose = T)
obj <- RunPCA(obj, npcs = 30, verbose = T)
obj <- IntegrateLayers(obj, method = HarmonyIntegration, orig.reduction = "pca", new.reduction = "harmony" )
obj <- RunUMAP(obj, reduction = "harmony", dims = 1:30, reduction.name = "umap.harmony", return.model = T, min.dist = 0.001)
obj <- JoinLayers(obj)
saveRDS( obj, 'species/human_mouse_myel_harmony_umap_joined.rds', compress = F ) 

### Human-Mouse Integration Lymph ----
obj_mouse <- readRDS("../newMouse/sketch_v5_all/final_figures/mouse_union_lymph_it4_res0.6_joined_annotated_sparse.rds")
obj_human <- readRDS("n0f/final_figures/human_lymph.rds")

library(nichenetr)

mouseToHumanNames = convert_mouse_to_human_symbols( rownames( obj_mouse ) )
ii = !(is.na(mouseToHumanNames))
counts.mat = GetAssayData( obj_mouse, slot = 'counts' )
counts.mat = counts.mat[ ii, ]
rownames(counts.mat) = mouseToHumanNames[ii]
counts.mat = as(counts.mat, Class = "dgCMatrix")
counts.mat = counts.mat[ !duplicated(rownames(counts.mat)), ]
options(Seurat.object.assay.version = "v5")
obj_mouse = CreateSeuratObject( counts.mat, meta.data = obj_mouse@meta.data )
obj_mouse[["RNA"]]$counts = as(obj_mouse[["RNA"]]$counts, Class = "dgCMatrix")
obj_mouse@meta.data <- obj_mouse@meta.data %>%
  select(any_of(c("orig.ident", "nCount_RNA","nFeature_RNA","study",
                  "organ","cancer","model","sorting","site","percent.mt",
                  "manual.annotation")))
write_matrix_dir(mat = obj_mouse[["RNA"]]$counts, dir = "species/mouse_toAnchor_lymph_counts", overwrite=T)
obj_mouse[["RNA"]]$counts <- open_matrix_dir(dir = "species/mouse_toAnchor_lymph_counts")
saveRDS(obj_mouse, "species/mouse_toAnchor_lymph_counts.rds",compress = F)

obj_human = CreateSeuratObject( obj_human[["RNA"]]$counts, meta.data = obj_human@meta.data )
obj_human[["RNA"]]$counts = as(obj_human[["RNA"]]$counts, Class = "dgCMatrix")
obj_human@meta.data <- obj_human@meta.data %>%
  select(any_of(c("orig.ident", "nCount_RNA","nFeature_RNA","study",
                  "organ","cancer","sorting","site","percent.mt",
                  "manual.annotation")))
write_matrix_dir(mat = obj_human[["RNA"]]$counts, dir = "species/human_toAnchor_lymph_counts", overwrite=T)
obj_human[["RNA"]]$counts <- open_matrix_dir(dir = "species/human_toAnchor_lymph_counts")
saveRDS(obj_human, "species/human_toAnchor_lymph_counts.rds",compress = F)

rm(list = ls())
obj_mouse = readRDS("species/mouse_toAnchor_lymph_counts.rds")
obj_human = readRDS("species/human_toAnchor_lymph_counts.rds")
data_mouse = GetAssayData(obj_mouse, slot = 'counts' )
data_human = GetAssayData(obj_human, slot = 'counts' )
commonGenes = intersect( rownames(data_mouse), rownames(data_human) )
data_mouse = data_mouse[ rownames(data_mouse) %in% commonGenes, ]
data_human = data_human[ rownames(data_human) %in% commonGenes, ]
obj_mouse = CreateSeuratObject( data_mouse, meta.data = obj_mouse@meta.data )
obj_human = CreateSeuratObject( data_human, meta.data = obj_human@meta.data )
saveRDS( obj_mouse, "species/mouse_toAnchor_lymph_counts_commonGenes.rds", compress = F )
saveRDS( obj_human, "species/human_toAnchor_lymph_counts_commonGenes.rds", compress = F )

obj_mouse$species = 'mouse'
obj_human$species = 'human'

obj_mouse$species.study <-  paste(obj_mouse$species, obj_mouse$study, sep = "_")
obj_human$species.study <-  paste(obj_human$species, obj_human$study, sep = "_")

obj = merge( obj_mouse, obj_human )
obj = JoinLayers( obj )
saveRDS( obj, "species/human_mouse_lymph_merged.rds", compress = F )

# run downstream analysis
obj <- NormalizeData(obj)
obj[["RNA"]] <- split(obj[["RNA"]], f = obj$species.study)
obj <- FindVariableFeatures(obj, verbose = T )
obj <- SketchData(object = obj, ncells = 15000, method = "LeverageScore", sketched.assay = "sketch")
DefaultAssay(obj) <- "sketch"
obj <- FindVariableFeatures(obj, verbose = T)
obj <- ScaleData(obj, verbose = F)
obj <- RunPCA(obj, npcs = 30, verbose = F)
# need to specify group.by = 'study'
obj <- IntegrateLayers(obj, group.by = "species.study", method = HarmonyIntegration, orig.reduction = "pca", new.reduction = "harmony", kmeans_init_nstart=20, kmeans_init_iter_max=5000, verbose = T)
obj <- RunUMAP(obj, reduction = "harmony", dims = 1:30, reduction.name = "umap.harmony", return.model = T, min.dist = 0.001)
saveRDS(obj, "species/human_mouse_lymph_var_sketch_pca_harmony_umap.rds", compress = F)

obj <- ProjectIntegration(object = obj, reduction = "harmony")
options(future.globals.maxSize = 8000 * 1024^2)
obj <- ProjectData(object = obj, sketched.reduction = "harmony.full", full.reduction = "harmony.full", umap.model = "umap.harmony", dims = 1:30)
saveRDS(obj, "species/human_mouse_lymph_var_sketch_pca_harmony_umap_project.rds",compress = F)

DefaultAssay(obj) <- "sketch"
obj <- JoinLayers(obj)
DefaultAssay(obj) <- "RNA"
obj <- JoinLayers(obj)
saveRDS(obj, "species/human_mouse_lymph_var_sketch_pca_harmony_umap_project_joined.rds", compress = F)

### Human-Mouse Integration Nerve ----
obj_mouse <- readRDS("../newMouse/sketch_v5_all/final_figures/it2_union_nerve_sketch1000fullharmonycluster_findvar_pca_harmonyumap_res0.03_annotated_sparse.rds")
obj_human <- readRDS("n0h/final_figures/human_nerve.rds")

library(nichenetr)

mouseToHumanNames = convert_mouse_to_human_symbols( rownames( obj_mouse ) )
ii = !(is.na(mouseToHumanNames))
counts.mat = GetAssayData( obj_mouse, slot = 'counts' )
counts.mat = counts.mat[ ii, ]
rownames(counts.mat) = mouseToHumanNames[ii]
counts.mat = as(counts.mat, Class = "dgCMatrix")
counts.mat = counts.mat[ !duplicated(rownames(counts.mat)), ]
options(Seurat.object.assay.version = "v5")
obj_mouse = CreateSeuratObject( counts.mat, meta.data = obj_mouse@meta.data )
obj_mouse[["RNA"]]$counts = as(obj_mouse[["RNA"]]$counts, Class = "dgCMatrix")
obj_mouse@meta.data <- obj_mouse@meta.data %>%
  select(any_of(c("orig.ident", "nCount_RNA","nFeature_RNA","study",
                  "organ","cancer","model","sorting","site","percent.mt",
                  "manual.annotation")))
write_matrix_dir(mat = obj_mouse[["RNA"]]$counts, dir = "species/mouse_toAnchor_nerve_counts", overwrite=T)
obj_mouse[["RNA"]]$counts <- open_matrix_dir(dir = "species/mouse_toAnchor_nerve_counts")
saveRDS(obj_mouse, "species/mouse_toAnchor_nerve_counts.rds",compress = F)

obj_human = CreateSeuratObject( obj_human[["RNA"]]$counts, meta.data = obj_human@meta.data )
obj_human[["RNA"]]$counts = as(obj_human[["RNA"]]$counts, Class = "dgCMatrix")
obj_human@meta.data <- obj_human@meta.data %>%
  select(any_of(c("orig.ident", "nCount_RNA","nFeature_RNA","study",
                  "organ","cancer","sorting","site","percent.mt",
                  "manual.annotation")))
write_matrix_dir(mat = obj_human[["RNA"]]$counts, dir = "species/human_toAnchor_nerve_counts", overwrite=T)
obj_human[["RNA"]]$counts <- open_matrix_dir(dir = "species/human_toAnchor_nerve_counts")
saveRDS(obj_human, "species/human_toAnchor_nerve_counts.rds",compress = F)

rm(list = ls())
obj_mouse = readRDS("species/mouse_toAnchor_nerve_counts.rds")
obj_human = readRDS("species/human_toAnchor_nerve_counts.rds")
data_mouse = GetAssayData(obj_mouse, slot = 'counts' )
data_human = GetAssayData(obj_human, slot = 'counts' )
commonGenes = intersect( rownames(data_mouse), rownames(data_human) )
data_mouse = data_mouse[ rownames(data_mouse) %in% commonGenes, ]
data_human = data_human[ rownames(data_human) %in% commonGenes, ]
obj_mouse = CreateSeuratObject( data_mouse, meta.data = obj_mouse@meta.data )
obj_human = CreateSeuratObject( data_human, meta.data = obj_human@meta.data )
saveRDS( obj_mouse, "species/mouse_toAnchor_nerve_counts_commonGenes.rds", compress = F )
saveRDS( obj_human, "species/human_toAnchor_nerve_counts_commonGenes.rds", compress = F )

obj_mouse$species = 'mouse'
obj_human$species = 'human'

obj_mouse$species.study <-  paste(obj_mouse$species, obj_mouse$study, sep = "_")
obj_human$species.study <-  paste(obj_human$species, obj_human$study, sep = "_")

obj = merge( obj_mouse, obj_human )
obj = JoinLayers( obj )
saveRDS( obj, "species/human_mouse_nerve_merged.rds", compress = F )

# run downstream analysis
obj <- NormalizeData(obj)
obj[["RNA"]] <- split(obj[["RNA"]], f = obj$species.study)
obj <- FindVariableFeatures(obj, verbose = T )
obj <- ScaleData(obj, verbose = T)
obj <- RunPCA(obj, npcs = 30, verbose = T)
obj <- IntegrateLayers(obj, method = HarmonyIntegration, orig.reduction = "pca", new.reduction = "harmony" )
obj <- RunUMAP(obj, reduction = "harmony", dims = 1:30, reduction.name = "umap.harmony", return.model = T, min.dist = 0.001)
obj <- JoinLayers(obj)
saveRDS( obj, 'species/human_mouse_nerve_harmony_umap_joined.rds', compress = F ) 
