### MOUSE META-ANALYSIS CODE
### Written by John Lu
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
  # download mouse gene names, symbols for GRCm39 using https://www.ensembl.org/biomart/martview/877296f74eb1c2a0b108975726c2386b
  synToName <- read.csv("GRCm39_geneNames_geneSynonyms.csv", header = T)
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
# download mouse gene names, symbols for GRCm39 using https://www.ensembl.org/biomart/martview/877296f74eb1c2a0b108975726c2386b
synToName <- read.csv("GRCm39_geneNames_geneSynonyms.csv", header = T)
counts.mat <- counts.mat[rownames(counts.mat) %in% synToName$Gene.name,]
counts.mat <- as(counts.mat, Class = "dgCMatrix")

### BPCell-ize object
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
options(Seurat.object.assay.version = "v5")

Sys.setenv( RETICULATE_PYTHON = "/home/johnlu/.conda/envs/scvi-env/bin/python" )
library(reticulate)

### sketch-based integration
obj <- readRDS("sketch_v5_all/unfiltered_counts_obj_findvar_sketch1000_pca.rds")
obj <- IntegrateLayers(
  object = obj, method = scVIIntegration,
  orig.reduction = "pca", new.reduction = "integrated.scvi", conda_env = '~/.conda/envs/scvi-env/', verbose = T
)
obj <- RunUMAP(obj, reduction = "integrated.scvi", dims = 1:30, reduction.name = "umap.scvi", return.model = T, min.dist = 0.001)
saveRDS(obj, "sketch_v5_all/unfiltered_counts_obj_findvar_sketch1000_scvi_umap.rds",compress = F)

### clustering
obj <- FindNeighbors(obj, reduction = "integrated.scvi", dims = 1:30)
for (i in seq(0.08, 0.09, 0.1)) {  #choose 0.08 since 0.07 does not distinguish mural vs fibro
  DefaultAssay(obj) = 'RNA'
  obj <- FindClusters(obj, resolution = i)
  DimPlot(obj, raster = TRUE)
  ggsave(paste0("sketch_v5_all/unfiltered_counts_obj_findvar_sketch1000_scvi_umap_res",i,".jpg"), width = 5, height = 5, units = "in", limitsize = FALSE)
  saveRDS(obj, paste0("sketch_v5_all/unfiltered_counts_obj_findvar_sketch1000_scvi_umap_res",i,".rds"), compress = FALSE)
}

### find markers
obj <- readRDS("sketch_v5_all/unfiltered_counts_obj_findvar_sketch1000_scvi_umap_res0.08.rds")
DefaultAssay(obj) <- "sketch"
obj <- JoinLayers(obj)
markers <- FindAllMarkers(obj, max.cells.per.ident = 1000, only.pos = T)
write.csv(markers, "sketch_v5_all/unfiltered_counts_obj_findvar_sketch1000_scvi_umap_res0.08_markers.csv")
markers$rank <- markers$pct.1 / (markers$pct.2 + 0.001) * markers$avg_log2FC
topmarkers <- markers %>%
  group_by(cluster) %>%
  slice_max(n = 100, order_by = rank)
write.csv(topmarkers, "sketch_v5_all/unfiltered_counts_obj_findvar_sketch1000_scvi_umap_res0.08_topmarkers.csv")

# res 0.08: cluster 0 = epi, 1 = fibro, 2 = myeloid, 3 = T cell, 4 = endo, 
# 5 = myeloid, 6 = B cells, 7 = mural, 8 = hepatocyte (epi), 9 = erythro (myeloid), 
# 10 = goblet cells (epi), 11 = lymphatic (endo), 12 = schwann (epi)
# 13 = ionocyte (epi), 14 = pDC (myeloid), 15 = mast (myeloid), 
# 16 = fibroblast (fibroblast), 17 = platelet (myeloid)

### project to full dataset
obj <- readRDS("sketch_v5_all/unfiltered_counts_obj_findvar_sketch1000_scvi_umap_res0.08.rds")
obj <- ProjectIntegration(object = obj, reduction = "integrated.scvi")
options(future.globals.maxSize = 8000 * 1024^2)
obj <- ProjectData(object = obj, sketched.reduction = "integrated.scvi.full", full.reduction = "integrated.scvi.full", umap.model = "umap.scvi", dims = 1:30, refdata = list(scvi.cluster.full = "seurat_clusters"))
saveRDS(obj, "sketch_v5_all/unfiltered_counts_obj_findvar_sketch1000_scvi_umap_res0.08_project.rds",compress = F)

### join layers
DefaultAssay(obj) <- "sketch"
obj <- JoinLayers(obj)
DefaultAssay(obj) <- "RNA"
obj <- JoinLayers(obj)
saveRDS(obj, "sketch_v5_all/unfiltered_counts_obj_findvar_sketch1000_scvi_umap_res0.08_project_joined.rds", compress = F)
saveRDS(obj$scvi.cluster.full, "sketch_v5_all/unfiltered_counts_obj_findvar_sketch1000_scvi_umap_res0.08_clusterList.rds")


### Integrate using Harmony----

obj <- readRDS("sketch_v5_all/unfiltered_counts_obj_findvar_sketch1000_pca.rds")
obj <- IntegrateLayers(obj, method = HarmonyIntegration, orig.reduction = "pca", new.reduction = "harmony", kmeans_init_nstart=20, kmeans_init_iter_max=5000, verbose = F)
obj <- RunUMAP(obj, reduction = "harmony", dims = 1:30, reduction.name = "umap.harmony", return.model = T, min.dist = 0.001)
saveRDS(obj, "sketch_v5_all/unfiltered_counts_obj_findvar_sketch1000_harmony_umap.rds",compress = F)

### clustering
obj <- FindNeighbors(obj, reduction = "harmony", dims = 1:30)
for (i in seq(0.01, 0.1, 0.01) ) { # choose 0.08
  obj <- FindClusters(obj, resolution = i)
  saveRDS(obj, paste0("sketch_v5_all/unfiltered_counts_obj_findvar_sketch1000_harmony_umap_res",i,".rds"), compress = FALSE)
}

### find markers
obj <- readRDS("sketch_v5_all/unfiltered_counts_obj_findvar_sketch1000_harmony_umap_res0.08.rds")
DefaultAssay(obj) <- "sketch"
obj <- JoinLayers(obj)
markers <- FindAllMarkers(obj, max.cells.per.ident = 1000, only.pos = T)
write.csv(markers, "sketch_v5_all/unfiltered_counts_obj_findvar_sketch1000_harmony_umap_res0.08_markers.csv")
markers$rank <- markers$pct.1 / (markers$pct.2 + 0.001) * markers$avg_log2FC
topmarkers <- markers %>%
  group_by(cluster) %>%
  slice_max(n = 100, order_by = rank)
write.csv(topmarkers, "sketch_v5_all/unfiltered_counts_obj_findvar_sketch1000_harmony_umap_res0.08_topmarkers.csv")

# cluster 0 = epi, 1 = fibro, 2 = myeloid?, 3 = T cells, 4 = endo,
# 5 = myeloid, 6 = B cells, 7 = pericyte, 8 = mesothelial, 9 = Schwann (epi)

### project to full dataset
obj <- readRDS("sketch_v5_all/unfiltered_counts_obj_findvar_sketch1000_harmony_umap_res0.08.rds")
obj <- ProjectIntegration(object = obj, reduction = "harmony")
options(future.globals.maxSize = 8000 * 1024^2)
obj <- ProjectData(object = obj, sketched.reduction = "harmony.full", full.reduction = "harmony.full", umap.model = "umap.harmony", dims = 1:30, refdata = list(harmony.cluster.full = "seurat_clusters"))
saveRDS(obj, "sketch_v5_all/unfiltered_counts_obj_findvar_sketch1000_harmony_umap_res0.08_project.rds",compress = F)

### join layers
DefaultAssay(obj) <- "sketch"
obj <- JoinLayers(obj)
DefaultAssay(obj) <- "RNA"
obj <- JoinLayers(obj)
saveRDS(obj, "sketch_v5_all/unfiltered_counts_obj_findvar_sketch1000_harmony_umap_res0.08_project_joined.rds", compress = F)
saveRDS(obj$harmony.cluster.full, "sketch_v5_all/unfiltered_counts_obj_findvar_sketch1000_harmony_umap_res0.08_clusterList.rds")

### Integrate using MNN----

obj <- readRDS("sketch_v5_all/unfiltered_counts_obj_findvar_sketch1000_pca.rds")
obj <- IntegrateLayers(
  object = obj, method = FastMNNIntegration,
  orig.reduction = "pca", new.reduction = "integrated.mnn", verbose = T
)
obj <- RunUMAP(obj, reduction = "integrated.mnn", dims = 1:30, reduction.name = "umap.mnn", return.model = T, min.dist = 0.001)
saveRDS(obj, "sketch_v5_all/unfiltered_counts_obj_findvar_sketch1000_mnn_umap.rds",compress = F)

### clustering
obj <- FindNeighbors(obj, reduction = "integrated.mnn", dims = 1:30)
DefaultAssay(obj) <- "RNA"
for (i in seq(0.05, 0.08, 0.01)) {  #choose 0.08
  obj <- FindClusters(obj, resolution = i)
  saveRDS(obj, paste0("sketch_v5_all/unfiltered_counts_obj_findvar_sketch1000_mnn_umap_res",i,".rds"), compress = FALSE)
}

### find markers
obj <- readRDS("sketch_v5_all/unfiltered_counts_obj_findvar_sketch1000_mnn_umap_res0.08.rds")
DefaultAssay(obj) <- "sketch"
obj <- JoinLayers(obj)
markers <- FindAllMarkers(obj, max.cells.per.ident = 1000, only.pos = T)
write.csv(markers, "sketch_v5_all/unfiltered_counts_obj_findvar_sketch1000_mnn_umap_res0.08_markers.csv")
markers$rank <- markers$pct.1 / (markers$pct.2 + 0.001) * markers$avg_log2FC
topmarkers <- markers %>%
  group_by(cluster) %>%
  slice_max(n = 100, order_by = rank)
write.csv(topmarkers, "sketch_v5_all/unfiltered_counts_obj_findvar_sketch1000_mnn_umap_res0.08_topmarkers.csv")

# cluster 0 = fibro, 1 = epi, 2 = myeloid, 3 = endo, 4 = T cell 
# 5 = myeloid, 6 = B cell, 7 = pericyte, 8 = neuron? (epi), 9 = acinar (epi)

### project
obj <- readRDS("sketch_v5_all/unfiltered_counts_obj_findvar_sketch1000_mnn_umap_res0.08.rds")
obj <- ProjectIntegration(object = obj, reduction = "integrated.mnn")
options(future.globals.maxSize = 8000 * 1024^2)
obj <- ProjectData(object = obj, sketched.reduction = "integrated.mnn.full", full.reduction = "integrated.mnn.full", umap.model = "umap.mnn", dims = 1:30, refdata = list(mnn.cluster.full = "seurat_clusters"))
saveRDS(obj, "sketch_v5_all/unfiltered_counts_obj_findvar_sketch1000_mnn_umap_res0.08_project.rds",compress = F)

### join layers
obj <- readRDS("sketch_v5_all/unfiltered_counts_obj_findvar_sketch1000_mnn_umap_res0.08.rds")
DefaultAssay(obj) <- "sketch"
obj <- JoinLayers(obj)
DefaultAssay(obj) <- "RNA"
obj <- JoinLayers(obj)
saveRDS(obj, "sketch_v5_all/unfiltered_counts_obj_findvar_sketch1000_mnn_umap_res0.08_project_joined.rds", compress = F)
saveRDS(obj$mnn.cluster.full, "sketch_v5_all/unfiltered_counts_obj_findvar_sketch1000_mnn_umap_res0.08_clusterList.rds")

### Fibroblasts iteration 1 -----

obj <- readRDS("sketch_v5_all/unfiltered_counts_obj_findvar_sketch1000_harmony_umap_res0.08_project.rds")
mnn <- readRDS("sketch_v5_all/unfiltered_counts_obj_findvar_sketch1000_mnn_umap_res0.08_clusterList.rds")
scvi <- readRDS("sketch_v5_all/unfiltered_counts_obj_findvar_sketch1000_scvi_umap_res0.08_clusterList.rds")
obj$mnn.cluster.full <- mnn
obj$scvi.cluster.full <- scvi

obj$fibro <- "false"
obj$fibro[obj$harmony.cluster.full %in% c(1,8)] <- "true"
obj$fibro[obj$mnn.cluster.full %in% c(0)] <- "true"
obj$fibro[obj$scvi.cluster.full %in% c(1,16)] <- "true"

DefaultAssay(obj) <- "RNA"
obj <- JoinLayers(obj)
obj <- subset(obj, subset = fibro == "true")

# re-generate seurat object and remove extraneous data
options(Seurat.object.assay.version = "v5")
obj <- CreateSeuratObject(counts = obj[["RNA"]]$counts, meta.data = obj@meta.data)
obj@meta.data <- obj@meta.data %>%
  select(any_of(c("orig.ident", "nCount_RNA","nFeature_RNA","study",
                  "organ","cancer","model","sorting","site","percent.mt","lin_neg")))

# re-perform BPcells
write_matrix_dir(mat = obj[["RNA"]]$counts, dir = "sketch_v5_all/union_fibro_sketch1000fullharmonycluster_counts")
counts.mat <- open_matrix_dir(dir = "sketch_v5_all/union_fibro_sketch1000fullharmonycluster_counts")
obj[["RNA"]]$counts <- counts.mat
saveRDS(obj, "sketch_v5_all/union_fibro_sketch1000fullharmonycluster_obj.rds",compress = F)

# run downstream analysis
obj <- NormalizeData(obj)
obj[["RNA"]] <- split(obj[["RNA"]], f = obj$study)
obj <- FindVariableFeatures(obj, verbose = FALSE)
obj <- ScaleData(obj, verbose = F)
obj <- RunPCA(obj, npcs = 30, verbose = F)
obj <- IntegrateLayers(obj, method = HarmonyIntegration, orig.reduction = "pca", new.reduction = "harmony", kmeans_init_nstart=20, kmeans_init_iter_max=5000, verbose = F)
obj <- RunUMAP(obj, reduction = "harmony", dims = 1:30, reduction.name = "umap.harmony", min.dist = 0.001, seed.use = 10000)
saveRDS(obj, "sketch_v5_all/union_fibro_sketch1000fullharmonycluster_findvar_pca_harmonyumap.rds",compress = F)

# clustering
obj <- FindNeighbors(obj, reduction = "harmony", dims = 1:30)
for (i in c(0.3, 0.1, 0.03, 0.01, 0.04, 0.05, 0.06)) {   # choose 0.04
  obj <- FindClusters(obj, resolution = i)
  saveRDS(obj, paste0("sketch_v5_all/union_fibro_sketch1000fullharmonycluster_findvar_pca_harmonyumap_res",i,".rds"), compress = FALSE)
}

# find markers
obj <- readRDS("sketch_v5_all/union_fibro_sketch1000fullharmonycluster_findvar_pca_harmonyumap_res0.04.rds")
DefaultAssay(obj) <- "RNA"
obj <- JoinLayers(obj)
obj[["RNA"]]$data <- as(object = obj[["RNA"]]$data, Class = "dgCMatrix")
markers <- FindAllMarkers(obj, max.cells.per.ident = 1000, only.pos = T)
write.csv(markers, "sketch_v5_all/union_fibro_sketch1000fullharmonycluster_findvar_pca_harmonyumap_res0.04_markers.csv")
markers$rank <- markers$pct.1 / (markers$pct.2 + 0.001) * markers$avg_log2FC
topmarkers <- markers %>%
  group_by(cluster) %>%
  slice_max(n = 100, order_by = rank)
write.csv(topmarkers, "sketch_v5_all/union_fibro_sketch1000fullharmonycluster_findvar_pca_harmonyumap_res0.04_topmarkers.csv")

# cluster 0 = fibro, 1 = tumor, 2 = immune, 3 = nerve, 4 = meso, 5 = endo, 6 = tumor

### Fibroblasts iteration 2 -----
obj <- readRDS("sketch_v5_all/union_fibro_sketch1000fullharmonycluster_findvar_pca_harmonyumap_res0.04.rds")

obj$fibromeso <- "false"
obj$fibromeso[obj$seurat_clusters %in% c(0, 4)] <- "true" # 0 fibro, 4 meso

DefaultAssay(obj) <- "RNA"
obj <- JoinLayers(obj)
obj <- subset(obj, subset = fibromeso == "true")
saveRDS(obj, "sketch_v5_all/union_fibro_sketch1000fullharmonycluster_findvar_pca_harmonyumap_res0.04_fibromesocutoff.rds", compress = F)

# re-generate seurat object and remove extraneous data
options(Seurat.object.assay.version = "v5")
obj <- CreateSeuratObject(counts = obj[["RNA"]]$counts, meta.data = obj@meta.data)
obj_backup <- obj
obj@meta.data <- obj@meta.data %>%
  select(any_of(c("orig.ident", "nCount_RNA","nFeature_RNA","study",
                  "organ","cancer","model","sorting","site","percent.mt","lin_neg")))

# re-perform BPcells
write_matrix_dir(mat = obj[["RNA"]]$counts, dir = "sketch_v5_all/it2_union_fibro_sketch1000fullharmonycluster_counts")
counts.mat <- open_matrix_dir(dir = "sketch_v5_all/it2_union_fibro_sketch1000fullharmonycluster_counts")
obj[["RNA"]]$counts <- counts.mat
saveRDS(obj, "sketch_v5_all/it2_union_fibro_sketch1000fullharmonycluster_obj.rds",compress = F)

# run downstream analysis
obj <- NormalizeData(obj)
obj[["RNA"]] <- split(obj[["RNA"]], f = obj$study)
obj <- FindVariableFeatures(obj, verbose = FALSE)
obj <- ScaleData(obj, verbose = F)
obj <- RunPCA(obj, npcs = 30, verbose = F)
obj <- IntegrateLayers(obj, method = HarmonyIntegration, orig.reduction = "pca", new.reduction = "harmony", kmeans_init_nstart=20, kmeans_init_iter_max=5000, verbose = F)
obj <- RunUMAP(obj, reduction = "harmony", dims = 1:30, reduction.name = "umap.harmony", min.dist = 0.001, seed.use = 10000)
saveRDS(obj, "sketch_v5_all/it2_union_fibro_sketch1000fullharmonycluster_findvar_pca_harmonyumap.rds",compress = F)

# clustering
obj <- FindNeighbors(obj, reduction = "harmony", dims = 1:30)
for (i in c(0.3, 0.2, 0.1, 0.09, 0.08, 0.07, 0.06, 0.05, 0.04)) {  # choose 0.05
  obj <- FindClusters(obj, resolution = i)
  DimPlot(obj, raster = TRUE)
  ggsave(paste0("sketch_v5_all/it2_union_fibro_sketch1000fullharmonycluster_findvar_pca_harmonyumap_res",i,".jpg"), width = 5, height = 5, units = "in", limitsize = FALSE)
  saveRDS(obj, paste0("sketch_v5_all/it2_union_fibro_sketch1000fullharmonycluster_findvar_pca_harmonyumap_res",i,".rds"), compress = FALSE)
}

# find clusters
obj <- readRDS("sketch_v5_all/it2_union_fibro_sketch1000fullharmonycluster_findvar_pca_harmonyumap_res0.05.rds")
obj <- JoinLayers(obj)
obj[["RNA"]]$data <- as(object = obj[["RNA"]]$data, Class = "dgCMatrix")
markers <- FindAllMarkers(obj, max.cells.per.ident = 1000, only.pos = T)
write.csv(markers, "sketch_v5_all/it2_union_fibro_sketch1000fullharmonycluster_findvar_pca_harmonyumap_res0.05_markers.csv")
markers$rank <- markers$pct.1 / (markers$pct.2 + 0.001) * markers$avg_log2FC
topmarkers <- markers %>%
  group_by(cluster) %>%
  slice_max(n = 100, order_by = rank)
write.csv(topmarkers, "sketch_v5_all/it2_union_fibro_sketch1000fullharmonycluster_findvar_pca_harmonyumap_res0.05_topmarkers.csv")

# cluster 0 = fibro, 1 = prolif, 2 = meso/apCAF, 3 = pericyte
# keep 0, 2


### Fibroblasts iteration 3 -----
obj <- readRDS("sketch_v5_all/it2_union_fibro_sketch1000fullharmonycluster_findvar_pca_harmonyumap_res0.05.rds")

obj$fibromeso <- "false"
obj$fibromeso[obj$seurat_clusters %in% c(0, 2)] <- "true" # 0 fibro, 1 meso, 2 mural, 3 = one study

DefaultAssay(obj) <- "RNA"
obj <- JoinLayers(obj)
obj <- subset(obj, subset = fibromeso == "true")
saveRDS(obj, "sketch_v5_all/it2_union_fibro_sketch1000fullharmonycluster_findvar_pca_harmonyumap_res0.05_fibromesocutoff.rds", compress = F)

# re-generate seurat object and remove extraneous data
options(Seurat.object.assay.version = "v5")
obj <- CreateSeuratObject(counts = obj[["RNA"]]$counts, meta.data = obj@meta.data)
obj_backup <- obj
obj@meta.data <- obj@meta.data %>%
  select(any_of(c("orig.ident", "nCount_RNA","nFeature_RNA","study",
                  "organ","cancer","model","sorting","site","percent.mt","lin_neg")))
# re-perform BPcells
write_matrix_dir(mat = obj[["RNA"]]$counts, dir = "sketch_v5_all/it3_union_fibro_sketch1000fullharmonycluster_counts")
counts.mat <- open_matrix_dir(dir = "sketch_v5_all/it3_union_fibro_sketch1000fullharmonycluster_counts")
obj[["RNA"]]$counts <- counts.mat
saveRDS(obj, "sketch_v5_all/it3_union_fibro_sketch1000fullharmonycluster_obj.rds",compress = F)

# run downstream analysis
obj <- NormalizeData(obj)
obj[["RNA"]] <- split(obj[["RNA"]], f = obj$study)
obj <- FindVariableFeatures(obj, verbose = FALSE)
obj <- ScaleData(obj, verbose = F)
obj <- RunPCA(obj, npcs = 30, verbose = F)
obj <- IntegrateLayers(obj, method = HarmonyIntegration, orig.reduction = "pca", new.reduction = "harmony", kmeans_init_nstart=20, kmeans_init_iter_max=5000, verbose = F)
obj <- RunUMAP(obj, reduction = "harmony", dims = 1:30, reduction.name = "umap.harmony", min.dist = 0.001, seed.use = 10000)
saveRDS(obj, "sketch_v5_all/it3_union_fibro_sketch1000fullharmonycluster_findvar_pca_harmonyumap.rds",compress = F)

# clustering
obj <- FindNeighbors(obj, reduction = "harmony", dims = 1:30)
for (i in seq(0.03, 0.25, 0.01)) { ## choose 0.1
  obj <- FindClusters(obj, resolution = i)
  saveRDS(obj, paste0("sketch_v5_all/it3_union_fibro_sketch1000fullharmonycluster_findvar_pca_harmonyumap_res",i,".rds"), compress = FALSE)
}

# find clusters
obj <- readRDS("sketch_v5_all/it3_union_fibro_sketch1000fullharmonycluster_findvar_pca_harmonyumap_res0.1.rds")
obj <- JoinLayers(obj)
obj[["RNA"]]$data <- as(object = obj[["RNA"]]$data, Class = "dgCMatrix")
markers <- FindAllMarkers(obj, max.cells.per.ident = 1000, only.pos = T)
write.csv(markers, "sketch_v5_all/it3_union_fibro_sketch1000fullharmonycluster_findvar_pca_harmonyumap_res0.1_markers.csv")
markers$rank <- markers$pct.1 / (markers$pct.2 + 0.001) * markers$avg_log2FC
topmarkers <- markers %>%
  group_by(cluster) %>%
  slice_max(n = 100, order_by = rank)
write.csv(topmarkers, "sketch_v5_all/it3_union_fibro_sketch1000fullharmonycluster_findvar_pca_harmonyumap_res0.1_topmarkers.csv")

### Fibroblasts iteration 4 -----
obj <- readRDS("sketch_v5_all/it3_union_fibro_sketch1000fullharmonycluster_findvar_pca_harmonyumap_res0.1.rds")

obj$fibromeso <- "false"
obj$fibromeso[obj$seurat_clusters %in% c(0,1,3)] <- "true" #2 is low-quallity

DefaultAssay(obj) <- "RNA"
obj <- JoinLayers(obj)
obj <- subset(obj, subset = fibromeso == "true")
saveRDS(obj, "sketch_v5_all/it3_union_fibro_sketch1000fullharmonycluster_findvar_pca_harmonyumap_res0.1_fibromesocutoff.rds", compress = F)

# re-generate seurat object and remove extraneous data
options(Seurat.object.assay.version = "v5")
obj <- CreateSeuratObject(counts = obj[["RNA"]]$counts, meta.data = obj@meta.data)
obj_backup <- obj
obj@meta.data <- obj@meta.data %>%
  select(any_of(c("orig.ident", "nCount_RNA","nFeature_RNA","study",
                  "organ","cancer","model","sorting","site","percent.mt","lin_neg")))
# re-perform BPcells
write_matrix_dir(mat = obj[["RNA"]]$counts, dir = "sketch_v5_all/it4_union_fibro_sketch1000fullharmonycluster_counts")
counts.mat <- open_matrix_dir(dir = "sketch_v5_all/it4_union_fibro_sketch1000fullharmonycluster_counts")
obj[["RNA"]]$counts <- counts.mat
saveRDS(obj, "sketch_v5_all/it4_union_fibro_sketch1000fullharmonycluster_obj.rds",compress = F)

# run downstream analysis
obj <- NormalizeData(obj)
obj[["RNA"]] <- split(obj[["RNA"]], f = obj$study)
obj <- FindVariableFeatures(obj, verbose = FALSE)
obj <- ScaleData(obj, verbose = F)
obj <- RunPCA(obj, npcs = 30, verbose = F)
obj <- IntegrateLayers(obj, method = HarmonyIntegration, orig.reduction = "pca", new.reduction = "harmony", kmeans_init_nstart=20, kmeans_init_iter_max=5000, verbose = F)
obj <- RunUMAP(obj, reduction = "harmony", dims = 1:30, reduction.name = "umap.harmony", min.dist = 0.001, seed.use = 10000)
saveRDS(obj, "sketch_v5_all/it4_union_fibro_sketch1000fullharmonycluster_findvar_pca_harmonyumap.rds",compress = F)

# clustering
obj <- FindNeighbors(obj, reduction = "harmony", dims = 1:30)
for (i in seq(0.03, 0.3, 0.01)) { # choose 0.26
  obj <- FindClusters(obj, resolution = i)
  saveRDS(obj, paste0("sketch_v5_all/it4_union_fibro_sketch1000fullharmonycluster_findvar_pca_harmonyumap_res",i,".rds"), compress = FALSE)
}

# find clusters
obj <- readRDS("sketch_v5_all/it4_union_fibro_sketch1000fullharmonycluster_findvar_pca_harmonyumap_res0.26.rds")
obj <- JoinLayers(obj)
saveRDS(obj, "sketch_v5_all/it4_union_fibro_sketch1000fullharmonycluster_findvar_pca_harmonyumap_res0.26_joined.rds", compress = F)
obj[["RNA"]]$data <- as(object = obj[["RNA"]]$data, Class = "dgCMatrix")
markers <- FindAllMarkers(obj, max.cells.per.ident = 1000, only.pos = T)
write.csv(markers, "sketch_v5_all/it4_union_fibro_sketch1000fullharmonycluster_findvar_pca_harmonyumap_res0.26_markers.csv")
markers$rank <- markers$pct.1 / (markers$pct.2 + 0.001) * markers$avg_log2FC
topmarkers <- markers %>%
  group_by(cluster) %>%
  slice_max(n = 100, order_by = rank)
write.csv(topmarkers, "sketch_v5_all/it4_union_fibro_sketch1000fullharmonycluster_findvar_pca_harmonyumap_res0.26_topmarkers.csv")

# save final object
obj <- readRDS("sketch_v5_all/it4_union_fibro_sketch1000fullharmonycluster_findvar_pca_harmonyumap_res0.26_joined.rds")
obj <- subset(obj, subset = seurat_clusters %in% c(0:7))
obj$manual.annotation <- "tbd"
obj$manual.annotation[obj$seurat_clusters == 0] <- "ssCAF Cxcl12+"
obj$manual.annotation[obj$seurat_clusters == 1] <- "mCAF Lrrc15+"
obj$manual.annotation[obj$seurat_clusters == 2] <- "ssCAF Pi16+"
obj$manual.annotation[obj$seurat_clusters == 3] <- "ssCAF Wif1+"
obj$manual.annotation[obj$seurat_clusters == 4] <- "iCAF Spp1+"
obj$manual.annotation[obj$seurat_clusters == 5] <- "apCAF Msln+"
obj$manual.annotation[obj$seurat_clusters == 6] <- "iCAF Isg15+"
obj$manual.annotation[obj$seurat_clusters == 7] <- "apCAF Cd74+"
saveRDS(obj, "sketch_v5_all/final_figures/mouse_union_fibro_it4_res0.26_joined_annotated.rds", compress = F)

# convert BP-cellized layers to sparse matrices
obj[["RNA"]]$counts <- as(obj[["RNA"]]$counts, Class = "dgCMatrix" )
obj[["RNA"]]$data <- as(obj[["RNA"]]$data, Class = "dgCMatrix" )
obj[["RNA"]]$scale.data <- as(obj[["RNA"]]$scale.data, Class = "dgCMatrix" )
saveRDS(obj, "sketch_v5_all/final_figures/mouse_union_fibro_it4_res0.26_joined_annotated_sparse.rds", compress = F)

### Endothelial iteration 1 -----
obj <- readRDS("sketch_v5_all/unfiltered_counts_obj_findvar_sketch1000_harmony_umap_res0.08_project.rds")
mnn <- readRDS("sketch_v5_all/unfiltered_counts_obj_findvar_sketch1000_mnn_umap_res0.08_clusterList.rds")
scvi <- readRDS("sketch_v5_all/unfiltered_counts_obj_findvar_sketch1000_scvi_umap_res0.08_clusterList.rds")
obj$mnn.cluster.full <- mnn
obj$scvi.cluster.full <- scvi

obj$endo <- "false"
obj$endo[obj$harmony.cluster.full %in% c(4)] <- "true"
obj$endo[obj$mnn.cluster.full %in% c(3)] <- "true"
obj$endo[obj$scvi.cluster.full %in% c(4, 11)] <- "true"

DefaultAssay(obj) <- "RNA"
obj <- JoinLayers(obj)
obj <- subset(obj, subset = endo == "true")
saveRDS(obj, "sketch_v5_all/unfiltered_counts_obj_findvar_sketch1000_harmony_umap_res0.08_project_endocutoff.rds", compress = F)

# re-generate seurat object and remove extraneous data
options(Seurat.object.assay.version = "v5")
obj <- CreateSeuratObject(counts = obj[["RNA"]]$counts, meta.data = obj@meta.data)
obj_backup <- obj
obj@meta.data <- obj@meta.data %>%
  select(any_of(c("orig.ident", "nCount_RNA","nFeature_RNA","study",
                  "organ","cancer","model","sorting","site","percent.mt","lin_neg")))
# re-perform BPcells
write_matrix_dir(mat = obj[["RNA"]]$counts, dir = "sketch_v5_all/union_endo_sketch1000fullharmonycluster_counts")
counts.mat <- open_matrix_dir(dir = "sketch_v5_all/union_endo_sketch1000fullharmonycluster_counts")
obj[["RNA"]]$counts <- counts.mat
saveRDS(obj, "sketch_v5_all/union_endo_sketch1000fullharmonycluster_obj.rds",compress = F)

# run downstream analysis
obj <- NormalizeData(obj)
obj[["RNA"]] <- split(obj[["RNA"]], f = obj$study)
obj <- FindVariableFeatures(obj, verbose = FALSE)
obj <- ScaleData(obj, verbose = F)
obj <- RunPCA(obj, npcs = 30, verbose = F)
obj <- IntegrateLayers(obj, method = HarmonyIntegration, orig.reduction = "pca", new.reduction = "harmony", kmeans_init_nstart=20, kmeans_init_iter_max=5000, verbose = F)
obj <- RunUMAP(obj, reduction = "harmony", dims = 1:30, reduction.name = "umap.harmony", min.dist = 0.001, seed.use = 10000)
saveRDS(obj, "sketch_v5_all/union_endo_sketch1000fullharmonycluster_findvar_pca_harmonyumap.rds",compress = F)

# clustering
obj <- FindNeighbors(obj, reduction = "harmony", dims = 1:30)
for (i in seq(0.03, 0.1, 0.01)) { # choose 0.05
  obj <- FindClusters(obj, resolution = i)
  saveRDS(obj, paste0("sketch_v5_all/union_endo_sketch1000fullharmonycluster_findvar_pca_harmonyumap_res",i,".rds"), compress = FALSE)
}

# find clusters
obj <- readRDS("sketch_v5_all/union_endo_sketch1000fullharmonycluster_findvar_pca_harmonyumap_res0.05.rds")
obj <- JoinLayers(obj)
obj[["RNA"]]$data <- as(object = obj[["RNA"]]$data, Class = "dgCMatrix")
markers <- FindAllMarkers(obj, max.cells.per.ident = 1000, only.pos = T)
write.csv(markers, "sketch_v5_all/union_endo_sketch1000fullharmonycluster_findvar_pca_harmonyumap_res0.05_markers.csv")
markers$rank <- markers$pct.1 / (markers$pct.2 + 0.001) * markers$avg_log2FC
topmarkers <- markers %>%
  group_by(cluster) %>%
  slice_max(n = 100, order_by = rank)
write.csv(topmarkers, "sketch_v5_all/union_endo_sketch1000fullharmonycluster_findvar_pca_harmonyumap_res0.05_topmarkers.csv")

# cluster 0 = endo, 1 = hepatocyte, 2 = endo, 3 = fibro, 4 = immune

### Endothelial iteration 2 -----
obj <- readRDS("sketch_v5_all/union_endo_sketch1000fullharmonycluster_findvar_pca_harmonyumap_res0.05.rds")

obj$endo <- "false"
obj$endo[obj$seurat_clusters %in% c(0, 2)] <- "true"
DimPlot(obj, group.by = "endo", raster = T, reduction = "umap.harmony")
ggsave("sketch_v5_all/union_endo_sketch1000fullharmonycluster_findvar_pca_harmonyumap_res0.05_endocutoff.jpg")

DefaultAssay(obj) <- "RNA"
obj <- JoinLayers(obj)
obj <- subset(obj, subset = endo == "true")
saveRDS(obj, "sketch_v5_all/union_endo_sketch1000fullharmonycluster_findvar_pca_harmonyumap_res0.05_endocutoff.rds", compress = F)

# re-generate seurat object and remove extraneous data
options(Seurat.object.assay.version = "v5")
obj <- CreateSeuratObject(counts = obj[["RNA"]]$counts, meta.data = obj@meta.data)
obj_backup <- obj
obj@meta.data <- obj@meta.data %>%
  select(any_of(c("orig.ident", "nCount_RNA","nFeature_RNA","study",
                  "organ","cancer","model","sorting","site","percent.mt","lin_neg")))
# re-perform BPcells
write_matrix_dir(mat = obj[["RNA"]]$counts, dir = "sketch_v5_all/it2_union_endo_sketch1000fullharmonycluster_counts")
counts.mat <- open_matrix_dir(dir = "sketch_v5_all/it2_union_endo_sketch1000fullharmonycluster_counts")
obj[["RNA"]]$counts <- counts.mat
saveRDS(obj, "sketch_v5_all/it2_union_endo_sketch1000fullharmonycluster_obj.rds",compress = F)

# run downstream analysis
obj <- NormalizeData(obj)
obj[["RNA"]] <- split(obj[["RNA"]], f = obj$study)
obj <- FindVariableFeatures(obj, verbose = FALSE)
obj <- ScaleData(obj, verbose = F)
obj <- RunPCA(obj, npcs = 30, verbose = F)
obj <- IntegrateLayers(obj, method = HarmonyIntegration, orig.reduction = "pca", new.reduction = "harmony", kmeans_init_nstart=20, kmeans_init_iter_max=5000, verbose = F)
obj <- RunUMAP(obj, reduction = "harmony", dims = 1:30, reduction.name = "umap.harmony", min.dist = 0.001, seed.use = 10000)
saveRDS(obj, "sketch_v5_all/it2_union_endo_sketch1000fullharmonycluster_findvar_pca_harmonyumap.rds",compress = F)

# clustering
obj <- FindNeighbors(obj, reduction = "harmony", dims = 1:30)
#for (i in c(0.3, 0.1, 0.03, 0.01)) {  
for (i in seq(0.03, 0.1, 0.01) ) {  # choose 0.04
  obj <- FindClusters(obj, resolution = i)
  saveRDS(obj, paste0("sketch_v5_all/it2_union_endo_sketch1000fullharmonycluster_findvar_pca_harmonyumap_res",i,".rds"), compress = FALSE)
}

# find clusters
obj <- readRDS("sketch_v5_all/it2_union_endo_sketch1000fullharmonycluster_findvar_pca_harmonyumap_res0.04.rds")
obj <- JoinLayers(obj)
obj[["RNA"]]$data <- as(object = obj[["RNA"]]$data, Class = "dgCMatrix")
markers <- FindAllMarkers(obj, max.cells.per.ident = 1000, only.pos = T)
write.csv(markers, "sketch_v5_all/it2_union_endo_sketch1000fullharmonycluster_findvar_pca_harmonyumap_res0.04_markers.csv")
markers$rank <- markers$pct.1 / (markers$pct.2 + 0.001) * markers$avg_log2FC
topmarkers <- markers %>%
  group_by(cluster) %>%
  slice_max(n = 100, order_by = rank)
write.csv(topmarkers, "sketch_v5_all/it2_union_endo_sketch1000fullharmonycluster_findvar_pca_harmonyumap_res0.04_topmarkers.csv")

# 0 = endo, 1 = lymphatic, 2 = tumor, 3 = prolif, 4 = liver sinusoidal cells

### Endothelial iteration 3 -----
obj <- readRDS("sketch_v5_all/it2_union_endo_sketch1000fullharmonycluster_findvar_pca_harmonyumap_res0.04.rds")

obj$endo <- "false"
obj$endo[obj$seurat_clusters %in% c(0, 1)] <- "true" # 0, 1 endo

DefaultAssay(obj) <- "RNA"
obj <- JoinLayers(obj)
obj <- subset(obj, subset = endo == "true")
saveRDS(obj, "sketch_v5_all/it2_union_endo_sketch1000fullharmonycluster_findvar_pca_harmonyumap_res0.04_endocutoff.rds", compress = F)

# re-generate seurat object and remove extraneous data
options(Seurat.object.assay.version = "v5")
obj <- CreateSeuratObject(counts = obj[["RNA"]]$counts, meta.data = obj@meta.data)
obj_backup <- obj
obj@meta.data <- obj@meta.data %>%
  select(any_of(c("orig.ident", "nCount_RNA","nFeature_RNA","study",
                  "organ","cancer","model","sorting","site","percent.mt","lin_neg")))
# re-perform BPcells
write_matrix_dir(mat = obj[["RNA"]]$counts, dir = "sketch_v5_all/it3_union_endo_sketch1000fullharmonycluster_counts")
counts.mat <- open_matrix_dir(dir = "sketch_v5_all/it3_union_endo_sketch1000fullharmonycluster_counts")
obj[["RNA"]]$counts <- counts.mat
saveRDS(obj, "sketch_v5_all/it3_union_endo_sketch1000fullharmonycluster_obj.rds",compress = F)

# run downstream analysis
obj <- NormalizeData(obj)
obj[["RNA"]] <- split(obj[["RNA"]], f = obj$study)
obj <- FindVariableFeatures(obj, verbose = FALSE)
obj <- ScaleData(obj, verbose = F)
obj <- RunPCA(obj, npcs = 30, verbose = F)
obj <- IntegrateLayers(obj, method = HarmonyIntegration, orig.reduction = "pca", new.reduction = "harmony", kmeans_init_nstart=20, kmeans_init_iter_max=5000, verbose = F)
obj <- RunUMAP(obj, reduction = "harmony", dims = 1:30, reduction.name = "umap.harmony", min.dist = 0.001, seed.use = 10000)
saveRDS(obj, "sketch_v5_all/it3_union_endo_sketch1000fullharmonycluster_findvar_pca_harmonyumap.rds",compress = F)

# clustering
obj <- FindNeighbors(obj, reduction = "harmony", dims = 1:30)
for (i in seq(0.03, 0.35, 0.01)) {  # choose 0.22 to resolve tip cells
  obj <- FindClusters(obj, resolution = i)
  saveRDS(obj, paste0("sketch_v5_all/it3_union_endo_sketch1000fullharmonycluster_findvar_pca_harmonyumap_res",i,".rds"), compress = FALSE)
}

# find clusters
obj <- readRDS("sketch_v5_all/it3_union_endo_sketch1000fullharmonycluster_findvar_pca_harmonyumap_res0.22.rds")
obj <- JoinLayers(obj)
obj[["RNA"]]$data <- as(object = obj[["RNA"]]$data, Class = "dgCMatrix")
markers <- FindAllMarkers(obj, max.cells.per.ident = 1000, only.pos = T)
write.csv(markers, "sketch_v5_all/it3_union_endo_sketch1000fullharmonycluster_findvar_pca_harmonyumap_res0.22_markers.csv")
markers$rank <- markers$pct.1 / (markers$pct.2 + 0.001) * markers$avg_log2FC
topmarkers <- markers %>%
  group_by(cluster) %>%
  slice_max(n = 100, order_by = rank)
write.csv(topmarkers, "sketch_v5_all/it3_union_endo_sketch1000fullharmonycluster_findvar_pca_harmonyumap_res0.22_topmarkers.csv")


obj <- readRDS("sketch_v5_all/it3_union_endo_sketch1000fullharmonycluster_findvar_pca_harmonyumap_res0.22.rds")
obj <- JoinLayers(obj)
obj <- subset(obj, subset = seurat_clusters %in% c(0:6))

obj$manual.annotation <- "tbd"
obj$manual.annotation[obj$seurat_clusters == 0] <- "Capillary Car4+"
obj$manual.annotation[obj$seurat_clusters == 1] <- "Capillary Car4+"
obj$manual.annotation[obj$seurat_clusters == 2] <- "Lymphatic Prox1+"
obj$manual.annotation[obj$seurat_clusters == 3] <- "Vein Ackr1+"
obj$manual.annotation[obj$seurat_clusters == 4] <- "Tip Esm1+"
obj$manual.annotation[obj$seurat_clusters == 5] <- "Artery Gja5+"
obj$manual.annotation[obj$seurat_clusters == 6] <- "Capillary Car4+"

obj[["RNA"]]$counts <- as(obj[["RNA"]]$counts, Class = "dgCMatrix" )
obj[["RNA"]]$data <- as(obj[["RNA"]]$data, Class = "dgCMatrix" )
obj[["RNA"]]$scale.data <- as(obj[["RNA"]]$scale.data, Class = "dgCMatrix" )

saveRDS(obj, "sketch_v5_all/final_figures/mouse_union_endo_it3_res0.22_joined_annotated_sparse.rds", compress = F)


### Mural iteration 1 -----
obj <- readRDS("sketch_v5_all/unfiltered_counts_obj_findvar_sketch1000_harmony_umap_res0.08_project.rds")
mnn <- readRDS("sketch_v5_all/unfiltered_counts_obj_findvar_sketch1000_mnn_umap_res0.08_clusterList.rds")
scvi <- readRDS("sketch_v5_all/unfiltered_counts_obj_findvar_sketch1000_scvi_umap_res0.08_clusterList.rds")
obj$mnn.cluster.full <- mnn
obj$scvi.cluster.full <- scvi

obj$mural <- "false"
obj$mural[obj$harmony.cluster.full %in% c(7)] <- "true"
obj$mural[obj$mnn.cluster.full %in% c(7)] <- "true"
obj$mural[obj$scvi.cluster.full %in% c(7)] <- "true"

DefaultAssay(obj) <- "RNA"
obj <- JoinLayers(obj)
obj <- subset(obj, subset = mural == "true")
saveRDS(obj, "sketch_v5_all/unfiltered_counts_obj_findvar_sketch1000_harmony_umap_res0.08_project_muralcutoff.rds", compress = F)

# re-generate seurat object and remove extraneous data
options(Seurat.object.assay.version = "v5")
obj <- CreateSeuratObject(counts = obj[["RNA"]]$counts, meta.data = obj@meta.data)
obj_backup <- obj
obj@meta.data <- obj@meta.data %>%
  select(any_of(c("orig.ident", "nCount_RNA","nFeature_RNA","study",
                  "organ","cancer","model","sorting","site","percent.mt","lin_neg")))
# re-perform BPcells
write_matrix_dir(mat = obj[["RNA"]]$counts, dir = "sketch_v5_all/union_mural_sketch1000fullharmonycluster_counts")
counts.mat <- open_matrix_dir(dir = "sketch_v5_all/union_mural_sketch1000fullharmonycluster_counts")
obj[["RNA"]]$counts <- counts.mat
saveRDS(obj, "sketch_v5_all/union_mural_sketch1000fullharmonycluster_obj.rds",compress = F)

# run downstream analysis
obj <- NormalizeData(obj)
obj[["RNA"]] <- split(obj[["RNA"]], f = obj$study)
obj <- FindVariableFeatures(obj, verbose = FALSE)
obj <- ScaleData(obj, verbose = F)
obj <- RunPCA(obj, npcs = 30, verbose = F)
obj <- IntegrateLayers(obj, method = HarmonyIntegration, orig.reduction = "pca", new.reduction = "harmony", kmeans_init_nstart=20, kmeans_init_iter_max=5000, verbose = F)
obj <- RunUMAP(obj, reduction = "harmony", dims = 1:30, reduction.name = "umap.harmony", min.dist = 0.001, seed.use = 10000)
saveRDS(obj, "sketch_v5_all/union_mural_sketch1000fullharmonycluster_findvar_pca_harmonyumap.rds",compress = F)

# clustering
obj <- FindNeighbors(obj, reduction = "harmony", dims = 1:30)
#for (i in c(0.3, 0.1, 0.03, 0.01)) { 
for (i in seq(0.04, 0.09, 0.01) ) { # choose 0.05
  obj <- FindClusters(obj, resolution = i)
  saveRDS(obj, paste0("sketch_v5_all/union_mural_sketch1000fullharmonycluster_findvar_pca_harmonyumap_res",i,".rds"), compress = FALSE)
}

# find markers
obj <- readRDS("sketch_v5_all/union_mural_sketch1000fullharmonycluster_findvar_pca_harmonyumap_res0.05.rds")
obj <- JoinLayers(obj)
obj[["RNA"]]$data <- as(object = obj[["RNA"]]$data, Class = "dgCMatrix")
markers <- FindAllMarkers(obj, max.cells.per.ident = 1000, only.pos = T)
write.csv(markers, "sketch_v5_all/union_mural_sketch1000fullharmonycluster_findvar_pca_harmonyumap_res0.05_markers.csv")
markers$rank <- markers$pct.1 / (markers$pct.2 + 0.001) * markers$avg_log2FC
topmarkers <- markers %>%
  group_by(cluster) %>%
  slice_max(n = 100, order_by = rank)
write.csv(topmarkers, "sketch_v5_all/union_mural_sketch1000fullharmonycluster_findvar_pca_harmonyumap_res0.05_topmarkers.csv")

# cluster 0 = pericyte , 1 = hepatocyte, 2 = fibroblast, 3 = smooth muscle, 4 = endo, 5 = prolif, 6 = stellate cell (fibro)

### Mural iteration 2 -----
obj <- readRDS("sketch_v5_all/union_mural_sketch1000fullharmonycluster_findvar_pca_harmonyumap_res0.05.rds")

obj$mural <- "false"
obj$mural[obj$seurat_clusters %in% c(0,3)] <- "true"

DefaultAssay(obj) <- "RNA"
obj <- JoinLayers(obj)
obj <- subset(obj, subset = mural == "true")
saveRDS(obj, "sketch_v5_all/union_mural_sketch1000fullharmonycluster_findvar_pca_harmonyumap_res0.05_muralcutoff.rds", compress = F)

# re-generate seurat object and remove extraneous data
options(Seurat.object.assay.version = "v5")
obj <- CreateSeuratObject(counts = obj[["RNA"]]$counts, meta.data = obj@meta.data)
obj_backup <- obj
obj@meta.data <- obj@meta.data %>%
  select(any_of(c("orig.ident", "nCount_RNA","nFeature_RNA","study",
                  "organ","cancer","model","sorting","site","percent.mt","lin_neg")))
# re-perform BPcells
write_matrix_dir(mat = obj[["RNA"]]$counts, dir = "sketch_v5_all/it2_union_mural_sketch1000fullharmonycluster_counts")
counts.mat <- open_matrix_dir(dir = "sketch_v5_all/it2_union_mural_sketch1000fullharmonycluster_counts")
obj[["RNA"]]$counts <- counts.mat
saveRDS(obj, "sketch_v5_all/it2_union_mural_sketch1000fullharmonycluster_obj.rds",compress = F)

# run downstream analysis
obj <- NormalizeData(obj)
obj[["RNA"]] <- split(obj[["RNA"]], f = obj$study)
obj <- FindVariableFeatures(obj, verbose = FALSE)
obj <- ScaleData(obj, verbose = F)
obj <- RunPCA(obj, npcs = 30, verbose = F)
obj <- IntegrateLayers(obj, method = HarmonyIntegration, orig.reduction = "pca", new.reduction = "harmony", kmeans_init_nstart=20, kmeans_init_iter_max=5000, verbose = F)
obj <- RunUMAP(obj, reduction = "harmony", dims = 1:30, reduction.name = "umap.harmony", min.dist = 0.001, seed.use = 10000)
saveRDS(obj, "sketch_v5_all/it2_union_mural_sketch1000fullharmonycluster_findvar_pca_harmonyumap.rds",compress = F)

# clustering
obj <- FindNeighbors(obj, reduction = "harmony", dims = 1:30)
#for (i in c(0.3, 0.1, 0.03, 0.01)) {  # choose 0.05
for (i in c(0.04, 0.05, 0.06, 0.07, 0.08, 0.09)) {
  obj <- FindClusters(obj, resolution = i)
  saveRDS(obj, paste0("sketch_v5_all/it2_union_mural_sketch1000fullharmonycluster_findvar_pca_harmonyumap_res",i,".rds"), compress = FALSE)
}

# find clusters
obj <- readRDS("sketch_v5_all/it2_union_mural_sketch1000fullharmonycluster_findvar_pca_harmonyumap_res0.05.rds")
obj <- JoinLayers(obj)
obj[["RNA"]]$data <- as(object = obj[["RNA"]]$data, Class = "dgCMatrix")
markers <- FindAllMarkers(obj, max.cells.per.ident = 1000, only.pos = T)
write.csv(markers, "sketch_v5_all/it2_union_mural_sketch1000fullharmonycluster_findvar_pca_harmonyumap_res0.05_markers.csv")
markers$rank <- markers$pct.1 / (markers$pct.2 + 0.001) * markers$avg_log2FC
topmarkers <- markers %>%
  group_by(cluster) %>%
  slice_max(n = 100, order_by = rank)
write.csv(topmarkers, "sketch_v5_all/it2_union_mural_sketch1000fullharmonycluster_findvar_pca_harmonyumap_res0.05_topmarkers.csv")


# cluster 0 = pericyte, 1 = vascular smooth muscle, 2 = hepatic stellate


### Mural iteration 3 -----
obj <- readRDS("sketch_v5_all/it2_union_mural_sketch1000fullharmonycluster_findvar_pca_harmonyumap_res0.05.rds")

obj$mural <- "false"
obj$mural[obj$seurat_clusters %in% c(0, 1)] <- "true" # 0, 1 mural

DefaultAssay(obj) <- "RNA"
obj <- JoinLayers(obj)
obj <- subset(obj, subset = mural == "true")
saveRDS(obj, "sketch_v5_all/it2_union_mural_sketch1000fullharmonycluster_findvar_pca_harmonyumap_res0.05_muralcutoff.rds", compress = F)

# re-generate seurat object and remove extraneous data
options(Seurat.object.assay.version = "v5")
obj <- CreateSeuratObject(counts = obj[["RNA"]]$counts, meta.data = obj@meta.data)
obj_backup <- obj
obj@meta.data <- obj@meta.data %>%
  select(any_of(c("orig.ident", "nCount_RNA","nFeature_RNA","study",
                  "organ","cancer","model","sorting","site","percent.mt","lin_neg")))
# re-perform BPcells
write_matrix_dir(mat = obj[["RNA"]]$counts, dir = "sketch_v5_all/it3_union_mural_sketch1000fullharmonycluster_counts")
counts.mat <- open_matrix_dir(dir = "sketch_v5_all/it3_union_mural_sketch1000fullharmonycluster_counts")
obj[["RNA"]]$counts <- counts.mat
saveRDS(obj, "sketch_v5_all/it3_union_mural_sketch1000fullharmonycluster_obj.rds",compress = F)

# run downstream analysis
obj <- NormalizeData(obj)
obj[["RNA"]] <- split(obj[["RNA"]], f = obj$study)
obj <- FindVariableFeatures(obj, verbose = FALSE)
obj <- ScaleData(obj, verbose = F)
obj <- RunPCA(obj, npcs = 30, verbose = F)
obj <- IntegrateLayers(obj, method = HarmonyIntegration, orig.reduction = "pca", new.reduction = "harmony", kmeans_init_nstart=20, kmeans_init_iter_max=5000, verbose = F)
obj <- RunUMAP(obj, reduction = "harmony", dims = 1:30, reduction.name = "umap.harmony", min.dist = 0.001, seed.use = 10000)
saveRDS(obj, "sketch_v5_all/it3_union_mural_sketch1000fullharmonycluster_findvar_pca_harmonyumap.rds",compress = F)

# clustering
obj <- FindNeighbors(obj, reduction = "harmony", dims = 1:30)
#for (i in c(0.3, 0.1, 0.03, 0.01)) {  # choose 0.05
#for (i in c(0.04, 0.05, 0.06, 0.07, 0.08, 0.09)) {
for (i in seq(0.12, 0.28, 0.02)) {
  obj <- FindClusters(obj, resolution = i)
  DimPlot(obj, raster = TRUE)
  ggsave(paste0("sketch_v5_all/it3_union_mural_sketch1000fullharmonycluster_findvar_pca_harmonyumap_res",i,".jpg"), width = 5, height = 5, units = "in", limitsize = FALSE)
  saveRDS(obj, paste0("sketch_v5_all/it3_union_mural_sketch1000fullharmonycluster_findvar_pca_harmonyumap_res",i,".rds"), compress = FALSE)
}

# find clusters
obj <- readRDS("sketch_v5_all/it3_union_mural_sketch1000fullharmonycluster_findvar_pca_harmonyumap_res0.05.rds")
obj <- JoinLayers(obj)
obj[["RNA"]]$data <- as(object = obj[["RNA"]]$data, Class = "dgCMatrix")
markers <- FindAllMarkers(obj, max.cells.per.ident = 1000, only.pos = T)
write.csv(markers, "sketch_v5_all/it3_union_mural_sketch1000fullharmonycluster_findvar_pca_harmonyumap_res0.05_markers.csv")
markers$rank <- markers$pct.1 / (markers$pct.2 + 0.001) * markers$avg_log2FC
topmarkers <- markers %>%
  group_by(cluster) %>%
  slice_max(n = 100, order_by = rank)
write.csv(topmarkers, "sketch_v5_all/it3_union_mural_sketch1000fullharmonycluster_findvar_pca_harmonyumap_res0.05_topmarkers.csv")

# save final object
obj <- readRDS("sketch_v5_all/it3_union_mural_sketch1000fullharmonycluster_findvar_pca_harmonyumap_res0.05.rds")
obj <- JoinLayers(obj)

obj$manual.annotation <- "tbd"
obj$manual.annotation[obj$seurat_clusters == 0] <- "Pericyte Ccl2+"
obj$manual.annotation[obj$seurat_clusters == 1] <- "Smooth muscle"
obj$manual.annotation[obj$seurat_clusters == 2] <- "Pericyte Cd248+"

saveRDS(obj, "sketch_v5_all/final_figures/mouse_union_mural_it3_res0.05_joined_annotated.rds", compress = F)

obj[["RNA"]]$counts <- as(obj[["RNA"]]$counts, Class = "dgCMatrix" )
obj[["RNA"]]$data <- as(obj[["RNA"]]$data, Class = "dgCMatrix" )
obj[["RNA"]]$scale.data <- as(obj[["RNA"]]$scale.data, Class = "dgCMatrix" )

saveRDS(obj, "sketch_v5_all/final_figures/mouse_union_mural_it3_res0.05_joined_annotated_sparse.rds", compress = F)


### Lymphoid iteration 1 -----
obj <- readRDS("sketch_v5_all/unfiltered_counts_obj_findvar_sketch1000_harmony_umap_res0.08_project.rds")
mnn <- readRDS("sketch_v5_all/unfiltered_counts_obj_findvar_sketch1000_mnn_umap_res0.08_clusterList.rds")
scvi <- readRDS("sketch_v5_all/unfiltered_counts_obj_findvar_sketch1000_scvi_umap_res0.08_clusterList.rds")
obj$mnn.cluster.full <- mnn
obj$scvi.cluster.full <- scvi

obj$lymph <- "false"
obj$lymph[obj$harmony.cluster.full %in% c(3,6)] <- "true"
obj$lymph[obj$mnn.cluster.full %in% c(4,6)] <- "true"
obj$lymph[obj$scvi.cluster.full %in% c(3,6)] <- "true"

DefaultAssay(obj) <- "RNA"
obj <- JoinLayers(obj)
obj <- subset(obj, subset = lymph == "true")
saveRDS(obj, "sketch_v5_all/unfiltered_counts_obj_findvar_sketch1000_harmony_umap_res0.08_project_lymphcutoff.rds", compress = F)

# re-generate seurat object and remove extraneous data
options(Seurat.object.assay.version = "v5")
obj <- CreateSeuratObject(counts = obj[["RNA"]]$counts, meta.data = obj@meta.data)
obj_backup <- obj
obj@meta.data <- obj@meta.data %>%
  select(any_of(c("orig.ident", "nCount_RNA","nFeature_RNA","study",
                  "organ","cancer","model","sorting","site","percent.mt","lin_neg")))
# re-perform BPcells
write_matrix_dir(mat = obj[["RNA"]]$counts, dir = "sketch_v5_all/union_lymph_sketch1000fullharmonycluster_counts")
counts.mat <- open_matrix_dir(dir = "sketch_v5_all/union_lymph_sketch1000fullharmonycluster_counts")
obj[["RNA"]]$counts <- counts.mat
saveRDS(obj, "sketch_v5_all/union_lymph_sketch1000fullharmonycluster_obj.rds",compress = F)

# run downstream analysis
obj <- NormalizeData(obj)
obj[["RNA"]] <- split(obj[["RNA"]], f = obj$study)
obj <- FindVariableFeatures(obj, verbose = FALSE)
obj <- ScaleData(obj, verbose = F)
obj <- RunPCA(obj, npcs = 30, verbose = F)
obj <- IntegrateLayers(obj, method = HarmonyIntegration, orig.reduction = "pca", new.reduction = "harmony", kmeans_init_nstart=20, kmeans_init_iter_max=5000, verbose = F)
obj <- RunUMAP(obj, reduction = "harmony", dims = 1:30, reduction.name = "umap.harmony", min.dist = 0.001, seed.use = 10000)
saveRDS(obj, "sketch_v5_all/union_lymph_sketch1000fullharmonycluster_findvar_pca_harmonyumap.rds",compress = F)

# clustering
obj <- FindNeighbors(obj, reduction = "harmony", dims = 1:30)
#for (i in c(0.3, 0.1, 0.03, 0.01)) { 
#for (i in c(0.09, 0.08, 0.07, 0.06, 0.05, 0.04)) { # choose 0.06
for (i in c(0.12, 0.14, 0.16, 0.18, 0.2)) {
  obj <- FindClusters(obj, resolution = i)
  DimPlot(obj, raster = TRUE)
  ggsave(paste0("sketch_v5_all/union_lymph_sketch1000fullharmonycluster_findvar_pca_harmonyumap_res",i,".jpg"), width = 5, height = 5, units = "in", limitsize = FALSE)
  saveRDS(obj, paste0("sketch_v5_all/union_lymph_sketch1000fullharmonycluster_findvar_pca_harmonyumap_res",i,".rds"), compress = FALSE)
}

# find clusters
obj <- readRDS("sketch_v5_all/union_lymph_sketch1000fullharmonycluster_findvar_pca_harmonyumap_res0.06.rds")
obj <- JoinLayers(obj)
obj[["RNA"]]$data <- as(object = obj[["RNA"]]$data, Class = "dgCMatrix")
markers <- FindAllMarkers(obj, max.cells.per.ident = 1000, only.pos = T)
write.csv(markers, "sketch_v5_all/union_lymph_sketch1000fullharmonycluster_findvar_pca_harmonyumap_res0.06_markers.csv")
markers$rank <- markers$pct.1 / (markers$pct.2 + 0.001) * markers$avg_log2FC
topmarkers <- markers %>%
  group_by(cluster) %>%
  slice_max(n = 100, order_by = rank)
write.csv(topmarkers, "sketch_v5_all/union_lymph_sketch1000fullharmonycluster_findvar_pca_harmonyumap_res0.06_topmarkers.csv")

# res 0.06: 0 = NK/T, 1 = B, 2 = Treg/T, 3 = myeloid, 4 = epi, 5 = prolif, 6 = fibro, 7 = plasma, 8 = endo, 9 = rbc, 10 = single study

### Lymphoid iteration 2 -----
obj <- readRDS("sketch_v5_all/union_lymph_sketch1000fullharmonycluster_findvar_pca_harmonyumap_res0.06.rds")

obj$lymph <- "false"
obj$lymph[obj$seurat_clusters %in% c(0,1,2,7)] <- "true"
DimPlot(obj, group.by = "lymph", raster = T, reduction = "umap.harmony")
ggsave("sketch_v5_all/union_lymph_sketch1000fullharmonycluster_findvar_pca_harmonyumap_res0.06_lymphcutoff.jpg")

DefaultAssay(obj) <- "RNA"
obj <- JoinLayers(obj)
obj <- subset(obj, subset = lymph == "true")
saveRDS(obj, "sketch_v5_all/union_lymph_sketch1000fullharmonycluster_findvar_pca_harmonyumap_res0.06_lymphcutoff.rds", compress = F)

# re-generate seurat object and remove extraneous data
options(Seurat.object.assay.version = "v5")
obj <- CreateSeuratObject(counts = obj[["RNA"]]$counts, meta.data = obj@meta.data)
obj_backup <- obj
obj@meta.data <- obj@meta.data %>%
  select(any_of(c("orig.ident", "nCount_RNA","nFeature_RNA","study",
                  "organ","cancer","model","sorting","site","percent.mt","lin_neg")))
# re-perform BPcells
write_matrix_dir(mat = obj[["RNA"]]$counts, dir = "sketch_v5_all/it2_union_lymph_sketch1000fullharmonycluster_counts")
counts.mat <- open_matrix_dir(dir = "sketch_v5_all/it2_union_lymph_sketch1000fullharmonycluster_counts")
obj[["RNA"]]$counts <- counts.mat
saveRDS(obj, "sketch_v5_all/it2_union_lymph_sketch1000fullharmonycluster_obj.rds",compress = F)

# run downstream analysis
obj <- NormalizeData(obj)
obj[["RNA"]] <- split(obj[["RNA"]], f = obj$study)
obj <- FindVariableFeatures(obj, verbose = FALSE)
obj <- ScaleData(obj, verbose = F)
obj <- RunPCA(obj, npcs = 30, verbose = F)
obj <- IntegrateLayers(obj, method = HarmonyIntegration, orig.reduction = "pca", new.reduction = "harmony", kmeans_init_nstart=20, kmeans_init_iter_max=5000, verbose = F)
obj <- RunUMAP(obj, reduction = "harmony", dims = 1:30, reduction.name = "umap.harmony", min.dist = 0.001, seed.use = 10000)
saveRDS(obj, "sketch_v5_all/it2_union_lymph_sketch1000fullharmonycluster_findvar_pca_harmonyumap.rds",compress = F)

# clustering
obj <- FindNeighbors(obj, reduction = "harmony", dims = 1:30)
#for (i in c(0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1, 0.12, 0.14, 0.16, 0.18, 0.2, 0.22, 0.24, 0.26, 0.28, 0.3)) {  
for (i in c(0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0)) {  # choose 0.8
  obj <- FindClusters(obj, resolution = i)
  DimPlot(obj, raster = TRUE, label = T)
  ggsave(paste0("sketch_v5_all/it2_union_lymph_sketch1000fullharmonycluster_findvar_pca_harmonyumap_res",i,".jpg"), width = 5, height = 5, units = "in", limitsize = FALSE)
  saveRDS(obj, paste0("sketch_v5_all/it2_union_lymph_sketch1000fullharmonycluster_findvar_pca_harmonyumap_res",i,".rds"), compress = FALSE)
}

# find clusters
obj <- readRDS("sketch_v5_all/it2_union_lymph_sketch1000fullharmonycluster_findvar_pca_harmonyumap_res0.8.rds")
obj <- JoinLayers(obj)
obj[["RNA"]]$data <- as(object = obj[["RNA"]]$data, Class = "dgCMatrix")
markers <- FindAllMarkers(obj, max.cells.per.ident = 1000, only.pos = T)
write.csv(markers, "sketch_v5_all/it2_union_lymph_sketch1000fullharmonycluster_findvar_pca_harmonyumap_res0.8_markers.csv")
markers$rank <- markers$pct.1 / (markers$pct.2 + 0.001) * markers$avg_log2FC
topmarkers <- markers %>%
  group_by(cluster) %>%
  slice_max(n = 100, order_by = rank)
write.csv(topmarkers, "sketch_v5_all/it2_union_lymph_sketch1000fullharmonycluster_findvar_pca_harmonyumap_res0.8_topmarkers.csv")

### Lymphoid iteration 3 -----
obj <- readRDS("sketch_v5_all/it2_union_lymph_sketch1000fullharmonycluster_findvar_pca_harmonyumap_res0.8.rds")

obj$lymph <- "false"
obj$lymph[obj$seurat_clusters %in% c(0:13, 15:18, 20)] <- "true"

DefaultAssay(obj) <- "RNA"
obj <- JoinLayers(obj)
obj <- subset(obj, subset = lymph == "true")
saveRDS(obj, "sketch_v5_all/it2_union_lymph_sketch1000fullharmonycluster_findvar_pca_harmonyumap_res0.8_lymphcutoff.rds", compress = F)

# re-generate seurat object and remove extraneous data
options(Seurat.object.assay.version = "v5")
obj <- CreateSeuratObject(counts = obj[["RNA"]]$counts, meta.data = obj@meta.data)
obj_backup <- obj
obj@meta.data <- obj@meta.data %>%
  select(any_of(c("orig.ident", "nCount_RNA","nFeature_RNA","study",
                  "organ","cancer","model","sorting","site","percent.mt","lin_neg")))
# re-perform BPcells
write_matrix_dir(mat = obj[["RNA"]]$counts, dir = "sketch_v5_all/it3_union_lymph_sketch1000fullharmonycluster_counts")
counts.mat <- open_matrix_dir(dir = "sketch_v5_all/it3_union_lymph_sketch1000fullharmonycluster_counts")
obj[["RNA"]]$counts <- counts.mat
saveRDS(obj, "sketch_v5_all/it3_union_lymph_sketch1000fullharmonycluster_obj.rds",compress = F)

# run downstream analysis
obj <- NormalizeData(obj)
obj[["RNA"]] <- split(obj[["RNA"]], f = obj$study)
obj <- FindVariableFeatures(obj, verbose = FALSE)
obj <- ScaleData(obj, verbose = F)
obj <- RunPCA(obj, npcs = 30, verbose = F)
obj <- IntegrateLayers(obj, method = HarmonyIntegration, orig.reduction = "pca", new.reduction = "harmony", kmeans_init_nstart=20, kmeans_init_iter_max=5000, verbose = F)
obj <- RunUMAP(obj, reduction = "harmony", dims = 1:30, reduction.name = "umap.harmony", min.dist = 0.001, seed.use = 10000)
saveRDS(obj, "sketch_v5_all/it3_union_lymph_sketch1000fullharmonycluster_findvar_pca_harmonyumap.rds",compress = F)

# clustering
obj <- FindNeighbors(obj, reduction = "harmony", dims = 1:30)
#for (i in c(0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1, 0.12, 0.14, 0.16, 0.18, 0.2, 0.22, 0.24, 0.26, 0.28, 0.3)) {  
for (i in c(0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0)) {  # choose 0.8
  obj <- FindClusters(obj, resolution = i)
  DimPlot(obj, raster = TRUE, label = T)
  ggsave(paste0("sketch_v5_all/it3_union_lymph_sketch1000fullharmonycluster_findvar_pca_harmonyumap_res",i,".jpg"), width = 5, height = 5, units = "in", limitsize = FALSE)
  saveRDS(obj, paste0("sketch_v5_all/it3_union_lymph_sketch1000fullharmonycluster_findvar_pca_harmonyumap_res",i,".rds"), compress = FALSE)
}

# find markers
obj <- readRDS("sketch_v5_all/it3_union_lymph_sketch1000fullharmonycluster_findvar_pca_harmonyumap_res0.8.rds")
obj <- JoinLayers(obj)
obj[["RNA"]]$data <- as(object = obj[["RNA"]]$data, Class = "dgCMatrix")
markers <- FindAllMarkers(obj, max.cells.per.ident = 1000, only.pos = T)
write.csv(markers, "sketch_v5_all/it3_union_lymph_sketch1000fullharmonycluster_findvar_pca_harmonyumap_res0.8_markers.csv")
markers$rank <- markers$pct.1 / (markers$pct.2 + 0.001) * markers$avg_log2FC
topmarkers <- markers %>%
  group_by(cluster) %>%
  slice_max(n = 100, order_by = rank)
write.csv(topmarkers, "sketch_v5_all/it3_union_lymph_sketch1000fullharmonycluster_findvar_pca_harmonyumap_res0.8_topmarkers.csv")

### Lymphoid iteration 4 -----
obj <- readRDS("sketch_v5_all/it3_union_lymph_sketch1000fullharmonycluster_findvar_pca_harmonyumap_res0.8.rds")

obj$lymph <- "false"
obj$lymph[obj$seurat_clusters %in% c(0,3:12, 14:16)] <- "true"
DimPlot(obj, group.by = "lymph", raster = T, reduction = "umap.harmony")
ggsave("sketch_v5_all/it3_union_lymph_sketch1000fullharmonycluster_findvar_pca_harmonyumap_res0.8_lymphcutoff.jpg")

obj@reductions <- list()
DefaultAssay(obj) <- "RNA"
obj[["RNA"]]$data <- NULL
obj[["sketch"]] <- NULL
obj <- JoinLayers(obj)
obj <- subset(obj, subset = lymph == "true")
saveRDS(obj, "sketch_v5_all/it3_union_lymph_sketch1000fullharmonycluster_findvar_pca_harmonyumap_res0.8_lymphcutoff.rds", compress = F)

# re-generate seurat object and remove extraneous data
options(Seurat.object.assay.version = "v5")
obj <- CreateSeuratObject(counts = obj[["RNA"]]$counts, meta.data = obj@meta.data)
obj_backup <- obj
obj@meta.data <- obj@meta.data %>%
  select(any_of(c("orig.ident", "nCount_RNA","nFeature_RNA","study",
                  "organ","cancer","model","sorting","site","percent.mt","lin_neg")))
# re-perform BPcells
write_matrix_dir(mat = obj[["RNA"]]$counts, dir = "sketch_v5_all/it4_union_lymph_sketch1000fullharmonycluster_counts")
counts.mat <- open_matrix_dir(dir = "sketch_v5_all/it4_union_lymph_sketch1000fullharmonycluster_counts")
obj[["RNA"]]$counts <- counts.mat
saveRDS(obj, "sketch_v5_all/it4_union_lymph_sketch1000fullharmonycluster_obj.rds",compress = F)

# run downstream analysis
obj <- NormalizeData(obj)
obj[["RNA"]] <- split(obj[["RNA"]], f = obj$study)
obj <- FindVariableFeatures(obj, verbose = FALSE)
obj <- ScaleData(obj, verbose = F)
obj <- RunPCA(obj, npcs = 30, verbose = F)
obj <- IntegrateLayers(obj, method = HarmonyIntegration, orig.reduction = "pca", new.reduction = "harmony", kmeans_init_nstart=20, kmeans_init_iter_max=5000, verbose = F)
obj <- RunUMAP(obj, reduction = "harmony", dims = 1:30, reduction.name = "umap.harmony", min.dist = 0.001, seed.use = 10000)
saveRDS(obj, "sketch_v5_all/it4_union_lymph_sketch1000fullharmonycluster_findvar_pca_harmonyumap.rds",compress = F)

# clustering
obj <- FindNeighbors(obj, reduction = "harmony", dims = 1:30)
#for (i in c(0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1, 0.12, 0.14, 0.16, 0.18, 0.2, 0.22, 0.24, 0.26, 0.28, 0.3)) {  
for (i in rev(c(0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0))) {  
  obj <- FindClusters(obj, resolution = i)
  saveRDS(obj, paste0("sketch_v5_all/it4_union_lymph_sketch1000fullharmonycluster_findvar_pca_harmonyumap_res",i,".rds"), compress = FALSE)
}

# find clusters
obj <- readRDS("sketch_v5_all/it4_union_lymph_sketch1000fullharmonycluster_findvar_pca_harmonyumap_res1.rds")
obj <- JoinLayers(obj)
obj[["RNA"]]$data <- as(object = obj[["RNA"]]$data, Class = "dgCMatrix")
markers <- FindAllMarkers(obj, max.cells.per.ident = 1000, only.pos = T)
write.csv(markers, "sketch_v5_all/it4_union_lymph_sketch1000fullharmonycluster_findvar_pca_harmonyumap_res1_markers.csv")
markers$rank <- markers$pct.1 / (markers$pct.2 + 0.001) * markers$avg_log2FC
topmarkers <- markers %>%
  group_by(cluster) %>%
  slice_max(n = 100, order_by = rank)
write.csv(topmarkers, "sketch_v5_all/it4_union_lymph_sketch1000fullharmonycluster_findvar_pca_harmonyumap_res1_topmarkers.csv")

obj <- readRDS("sketch_v5_all/it4_union_lymph_sketch1000fullharmonycluster_findvar_pca_harmonyumap_res0.6.rds")
obj <- JoinLayers(obj)

# manual annotations for entire dataset
obj$manual.annotation <- "tbd"
obj$manual.annotation[obj$seurat_clusters == 0] <- "B Naive Cd19+"
obj$manual.annotation[obj$seurat_clusters == 1] <- "NK Klrk1+"
obj$manual.annotation[obj$seurat_clusters == 2] <- "Treg Foxp3+"
obj$manual.annotation[obj$seurat_clusters == 3] <- "T Cd4+ Naive"
obj$manual.annotation[obj$seurat_clusters == 4] <- "T Cd8+ Naive"
obj$manual.annotation[obj$seurat_clusters == 5] <- "Th17 Il17a+"
obj$manual.annotation[obj$seurat_clusters == 6] <- "T Cd8+ Naive"
obj$manual.annotation[obj$seurat_clusters == 7] <- "Plasma Jchain+"
obj$manual.annotation[obj$seurat_clusters == 8] <- "Tex Pdcd1+"
obj$manual.annotation[obj$seurat_clusters == 9] <- "B Pre Vpreb3+"
obj$manual.annotation[obj$seurat_clusters == 10] <- "B Naive Cd19+"
obj$manual.annotation[obj$seurat_clusters == 11] <- "T Cd8+/Isg15+"
obj$manual.annotation[obj$seurat_clusters == 12] <- "B Mem Zbtb32+"
obj$manual.annotation[obj$seurat_clusters == 13] <- "Plasma Jchain+"
obj$manual.annotation[obj$seurat_clusters == 14] <- "Th2 Gata3+"
obj$manual.annotation[obj$seurat_clusters == 15] <- "B Naive Cd19+"
obj$manual.annotation[obj$seurat_clusters == 16] <- "B Isg15+"
obj$manual.annotation[obj$seurat_clusters == 17] <- "T Cd8+ Naive"

saveRDS(obj, "sketch_v5_all/final_figures/mouse_union_lymph_it4_res0.6_joined_annotated.rds", compress = F)

obj[["RNA"]]$counts <- as(obj[["RNA"]]$counts, Class = "dgCMatrix" )
obj[["RNA"]]$data <- as(obj[["RNA"]]$data, Class = "dgCMatrix" )
obj[["RNA"]]$scale.data <- as(obj[["RNA"]]$scale.data, Class = "dgCMatrix" )

saveRDS(obj, "sketch_v5_all/final_figures/mouse_union_lymph_it4_res0.6_joined_annotated_sparse.rds", compress = F)


### Myeloid iteration 1 -----
obj <- readRDS("sketch_v5_all/unfiltered_counts_obj_findvar_sketch1000_harmony_umap_res0.08_project.rds")
mnn <- readRDS("sketch_v5_all/unfiltered_counts_obj_findvar_sketch1000_mnn_umap_res0.08_clusterList.rds")
scvi <- readRDS("sketch_v5_all/unfiltered_counts_obj_findvar_sketch1000_scvi_umap_res0.08_clusterList.rds")
obj$mnn.cluster.full <- mnn
obj$scvi.cluster.full <- scvi

obj$myel <- "false"
obj$myel[obj$harmony.cluster.full %in% c(2,5)] <- "true"
obj$myel[obj$mnn.cluster.full %in% c(2,5)] <- "true"
obj$myel[obj$scvi.cluster.full %in% c(2,5,14,17)] <- "true"
DimPlot(obj, group.by = "myel", raster = T, reduction = "full.umap.harmony", alpha = 0.05)
ggsave("sketch_v5_all/unfiltered_counts_obj_findvar_sketch1000_harmony_umap_res0.08_project_myelcutoff.jpg")

DefaultAssay(obj) <- "RNA"
obj <- JoinLayers(obj)
obj <- subset(obj, subset = myel == "true")
saveRDS(obj, "sketch_v5_all/unfiltered_counts_obj_findvar_sketch1000_harmony_umap_res0.08_project_myelcutoff.rds", compress = F)

# re-generate seurat object and remove extraneous data
options(Seurat.object.assay.version = "v5")
obj <- CreateSeuratObject(counts = obj[["RNA"]]$counts, meta.data = obj@meta.data)
obj_backup <- obj
obj@meta.data <- obj@meta.data %>%
  select(any_of(c("orig.ident", "nCount_RNA","nFeature_RNA","study",
                  "organ","cancer","model","sorting","site","percent.mt","lin_neg")))
# re-perform BPcells
write_matrix_dir(mat = obj[["RNA"]]$counts, dir = "sketch_v5_all/union_myel_sketch1000fullharmonycluster_counts")
counts.mat <- open_matrix_dir(dir = "sketch_v5_all/union_myel_sketch1000fullharmonycluster_counts")
obj[["RNA"]]$counts <- counts.mat
saveRDS(obj, "sketch_v5_all/union_myel_sketch1000fullharmonycluster_obj.rds",compress = F)

# run downstream analysis
obj <- NormalizeData(obj)
obj[["RNA"]] <- split(obj[["RNA"]], f = obj$study)
obj <- FindVariableFeatures(obj, verbose = FALSE)
obj <- ScaleData(obj, verbose = F)
obj <- RunPCA(obj, npcs = 30, verbose = F)
obj <- IntegrateLayers(obj, method = HarmonyIntegration, orig.reduction = "pca", new.reduction = "harmony", kmeans_init_nstart=20, kmeans_init_iter_max=5000, verbose = F)
obj <- RunUMAP(obj, reduction = "harmony", dims = 1:30, reduction.name = "umap.harmony", min.dist = 0.001, seed.use = 10000)
saveRDS(obj, "sketch_v5_all/union_myel_sketch1000fullharmonycluster_findvar_pca_harmonyumap.rds",compress = F)

# clustering
obj <- FindNeighbors(obj, reduction = "harmony", dims = 1:30)
#for (i in c(0.3, 0.1, 0.03, 0.01)) { 
#for (i in c(0.09, 0.08, 0.07, 0.06, 0.05, 0.04)) { # choose 0.08
for (i in c(0.1, 0.12, 0.14, 0.16, 0.18, 0.2)) {
  obj <- FindClusters(obj, resolution = i)
  DimPlot(obj, raster = TRUE)
  ggsave(paste0("sketch_v5_all/union_myel_sketch1000fullharmonycluster_findvar_pca_harmonyumap_res",i,".jpg"), width = 5, height = 5, units = "in", limitsize = FALSE)
  saveRDS(obj, paste0("sketch_v5_all/union_myel_sketch1000fullharmonycluster_findvar_pca_harmonyumap_res",i,".rds"), compress = FALSE)
}

# find clusters
obj <- readRDS("sketch_v5_all/union_myel_sketch1000fullharmonycluster_findvar_pca_harmonyumap_res0.08.rds")
obj <- JoinLayers(obj)
obj[["RNA"]]$data <- as(object = obj[["RNA"]]$data, Class = "dgCMatrix")
markers <- FindAllMarkers(obj, max.cells.per.ident = 1000, only.pos = T)
write.csv(markers, "sketch_v5_all/union_myel_sketch1000fullharmonycluster_findvar_pca_harmonyumap_res0.08_markers.csv")
markers$rank <- markers$pct.1 / (markers$pct.2 + 0.001) * markers$avg_log2FC
topmarkers <- markers %>%
  group_by(cluster) %>%
  slice_max(n = 100, order_by = rank)
write.csv(topmarkers, "sketch_v5_all/union_myel_sketch1000fullharmonycluster_findvar_pca_harmonyumap_res0.08_topmarkers.csv")

### Myeloid iteration 2 -----
obj <- readRDS("sketch_v5_all/union_myel_sketch1000fullharmonycluster_findvar_pca_harmonyumap_res0.08.rds")

obj$myel <- "false"
obj$myel[obj$seurat_clusters %in% c(0,1,3)] <- "true"

DefaultAssay(obj) <- "RNA"
obj <- JoinLayers(obj)
obj <- subset(obj, subset = myel == "true")
saveRDS(obj, "sketch_v5_all/union_myel_sketch1000fullharmonycluster_findvar_pca_harmonyumap_res0.08_myelcutoff.rds", compress = F)

# re-generate seurat object and remove extraneous data
options(Seurat.object.assay.version = "v5")
obj <- CreateSeuratObject(counts = obj[["RNA"]]$counts, meta.data = obj@meta.data)
obj_backup <- obj
obj@meta.data <- obj@meta.data %>%
  select(any_of(c("orig.ident", "nCount_RNA","nFeature_RNA","study",
                  "organ","cancer","model","sorting","site","percent.mt","lin_neg")))
# re-perform BPcells
write_matrix_dir(mat = obj[["RNA"]]$counts, dir = "sketch_v5_all/it2_union_myel_sketch1000fullharmonycluster_counts")
counts.mat <- open_matrix_dir(dir = "sketch_v5_all/it2_union_myel_sketch1000fullharmonycluster_counts")
obj[["RNA"]]$counts <- counts.mat
saveRDS(obj, "sketch_v5_all/it2_union_myel_sketch1000fullharmonycluster_obj.rds",compress = F)

# run downstream analysis
obj <- NormalizeData(obj)
obj[["RNA"]] <- split(obj[["RNA"]], f = obj$study)
obj <- FindVariableFeatures(obj, verbose = FALSE)
obj <- ScaleData(obj, verbose = F)
obj <- RunPCA(obj, npcs = 30, verbose = F)
obj <- IntegrateLayers(obj, method = HarmonyIntegration, orig.reduction = "pca", new.reduction = "harmony", kmeans_init_nstart=20, kmeans_init_iter_max=5000, verbose = F)
obj <- RunUMAP(obj, reduction = "harmony", dims = 1:30, reduction.name = "umap.harmony", min.dist = 0.001, seed.use = 10000)
saveRDS(obj, "sketch_v5_all/it2_union_myel_sketch1000fullharmonycluster_findvar_pca_harmonyumap.rds",compress = F)

# clustering
obj <- FindNeighbors(obj, reduction = "harmony", dims = 1:30)
for (i in c(seq(0.03, 0.1, 0.01), seq(0.12, 0.3, 0.02))) {  
  obj <- FindClusters(obj, resolution = i)
  saveRDS(obj, paste0("sketch_v5_all/it2_union_myel_sketch1000fullharmonycluster_findvar_pca_harmonyumap_res",i,".rds"), compress = FALSE)
}

# find clusters
obj <- readRDS("sketch_v5_all/it2_union_myel_sketch1000fullharmonycluster_findvar_pca_harmonyumap_res0.2.rds")
obj <- JoinLayers(obj)
obj[["RNA"]]$data <- as(object = obj[["RNA"]]$data, Class = "dgCMatrix")
markers <- FindAllMarkers(obj, max.cells.per.ident = 1000, only.pos = T)
write.csv(markers, "sketch_v5_all/it2_union_myel_sketch1000fullharmonycluster_findvar_pca_harmonyumap_res0.2_markers.csv")
markers$rank <- markers$pct.1 / (markers$pct.2 + 0.001) * markers$avg_log2FC
topmarkers <- markers %>%
  group_by(cluster) %>%
  slice_max(n = 100, order_by = rank)
write.csv(topmarkers, "sketch_v5_all/it2_union_myel_sketch1000fullharmonycluster_findvar_pca_harmonyumap_res0.2_topmarkers.csv")

# 2 is low qual --> remove

### Myeloid iteration 3 -----
obj <- readRDS("sketch_v5_all/it2_union_myel_sketch1000fullharmonycluster_findvar_pca_harmonyumap_res0.2.rds")

obj$myel <- "false"
obj$myel[obj$seurat_clusters != "2"] <- "true"

DefaultAssay(obj) <- "RNA"
obj <- JoinLayers(obj)
obj <- subset(obj, subset = myel == "true")
saveRDS(obj, "sketch_v5_all/it2_union_myel_sketch1000fullharmonycluster_findvar_pca_harmonyumap_res0.2_myelcutoff.rds", compress = F)

# re-generate seurat object and remove extraneous data
options(Seurat.object.assay.version = "v5")
obj <- CreateSeuratObject(counts = obj[["RNA"]]$counts, meta.data = obj@meta.data)
obj_backup <- obj
obj@meta.data <- obj@meta.data %>%
  select(any_of(c("orig.ident", "nCount_RNA","nFeature_RNA","study",
                  "organ","cancer","model","sorting","site","percent.mt","lin_neg")))
# re-perform BPcells
write_matrix_dir(mat = obj[["RNA"]]$counts, dir = "sketch_v5_all/it3_union_myel_sketch1000fullharmonycluster_counts")
counts.mat <- open_matrix_dir(dir = "sketch_v5_all/it3_union_myel_sketch1000fullharmonycluster_counts")
obj[["RNA"]]$counts <- counts.mat
saveRDS(obj, "sketch_v5_all/it3_union_myel_sketch1000fullharmonycluster_obj.rds",compress = F)

# run downstream analysis
obj <- NormalizeData(obj)
obj[["RNA"]] <- split(obj[["RNA"]], f = obj$study)
obj <- FindVariableFeatures(obj, verbose = FALSE)
obj <- ScaleData(obj, verbose = F)
obj <- RunPCA(obj, npcs = 30, verbose = F)
obj <- IntegrateLayers(obj, method = HarmonyIntegration, orig.reduction = "pca", new.reduction = "harmony", kmeans_init_nstart=20, kmeans_init_iter_max=5000, verbose = F)
obj <- RunUMAP(obj, reduction = "harmony", dims = 1:30, reduction.name = "umap.harmony", min.dist = 0.001, seed.use = 10000)
saveRDS(obj, "sketch_v5_all/it3_union_myel_sketch1000fullharmonycluster_findvar_pca_harmonyumap.rds",compress = F)

# clustering
obj <- FindNeighbors(obj, reduction = "harmony", dims = 1:30)
for (i in c(1.0, 0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1, 0.03)) {  # choose 0.4
  obj <- FindClusters(obj, resolution = i)
  DimPlot(obj, raster = TRUE, label = T)
  ggsave(paste0("sketch_v5_all/it3_union_myel_sketch1000fullharmonycluster_findvar_pca_harmonyumap_res",i,".jpg"), width = 5, height = 5, units = "in", limitsize = FALSE)
  saveRDS(obj, paste0("sketch_v5_all/it3_union_myel_sketch1000fullharmonycluster_findvar_pca_harmonyumap_res",i,".rds"), compress = FALSE)
}

# find clusters
obj <- readRDS("sketch_v5_all/it3_union_myel_sketch1000fullharmonycluster_findvar_pca_harmonyumap_res0.4.rds")
obj <- JoinLayers(obj)
obj[["RNA"]]$data <- as(object = obj[["RNA"]]$data, Class = "dgCMatrix")
markers <- FindAllMarkers(obj, max.cells.per.ident = 1000, only.pos = T)
write.csv(markers, "sketch_v5_all/it3_union_myel_sketch1000fullharmonycluster_findvar_pca_harmonyumap_res0.4_markers.csv")
markers$rank <- markers$pct.1 / (markers$pct.2 + 0.001) * markers$avg_log2FC
topmarkers <- markers %>%
  group_by(cluster) %>%
  slice_max(n = 100, order_by = rank)
write.csv(topmarkers, "sketch_v5_all/it3_union_myel_sketch1000fullharmonycluster_findvar_pca_harmonyumap_res0.4_topmarkers.csv")


### Myeloid iteration 4 -----
obj <- readRDS("sketch_v5_all/it3_union_myel_sketch1000fullharmonycluster_findvar_pca_harmonyumap_res0.4.rds")

obj$myel <- "false"
obj$myel[obj$seurat_clusters %in% c(0:6, 8:15)] <- "true"

DefaultAssay(obj) <- "RNA"
obj <- JoinLayers(obj)
obj <- subset(obj, subset = myel == "true")
saveRDS(obj, "sketch_v5_all/it3_union_myel_sketch1000fullharmonycluster_findvar_pca_harmonyumap_res0.4_myelcutoff.rds", compress = F)

# re-generate seurat object and remove extraneous data
options(Seurat.object.assay.version = "v5")
obj <- CreateSeuratObject(counts = obj[["RNA"]]$counts, meta.data = obj@meta.data)
obj_backup <- obj
obj@meta.data <- obj@meta.data %>%
  select(any_of(c("orig.ident", "nCount_RNA","nFeature_RNA","study",
                  "organ","cancer","model","sorting","site","percent.mt","lin_neg")))
# re-perform BPcells
write_matrix_dir(mat = obj[["RNA"]]$counts, dir = "sketch_v5_all/it4_union_myel_sketch1000fullharmonycluster_counts")
counts.mat <- open_matrix_dir(dir = "sketch_v5_all/it4_union_myel_sketch1000fullharmonycluster_counts")
obj[["RNA"]]$counts <- counts.mat
saveRDS(obj, "sketch_v5_all/it4_union_myel_sketch1000fullharmonycluster_obj.rds",compress = F)

# run downstream analysis
obj <- NormalizeData(obj)
obj[["RNA"]] <- split(obj[["RNA"]], f = obj$study)
obj <- FindVariableFeatures(obj, verbose = FALSE)
obj <- ScaleData(obj, verbose = F)
obj <- RunPCA(obj, npcs = 30, verbose = F)
obj <- IntegrateLayers(obj, method = HarmonyIntegration, orig.reduction = "pca", new.reduction = "harmony", kmeans_init_nstart=20, kmeans_init_iter_max=5000, verbose = F)
obj <- RunUMAP(obj, reduction = "harmony", dims = 1:30, reduction.name = "umap.harmony", min.dist = 0.001, seed.use = 10000)
saveRDS(obj, "sketch_v5_all/it4_union_myel_sketch1000fullharmonycluster_findvar_pca_harmonyumap.rds",compress = F)

# clustering
obj <- FindNeighbors(obj, reduction = "harmony", dims = 1:30)
for (i in c(1.0, 0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1, 0.03)) {  # 0.6
  obj <- FindClusters(obj, resolution = i)
  DimPlot(obj, raster = TRUE, label = T)
  ggsave(paste0("sketch_v5_all/it4_union_myel_sketch1000fullharmonycluster_findvar_pca_harmonyumap_res",i,".jpg"), width = 5, height = 5, units = "in", limitsize = FALSE)
  saveRDS(obj, paste0("sketch_v5_all/it4_union_myel_sketch1000fullharmonycluster_findvar_pca_harmonyumap_res",i,".rds"), compress = FALSE)
}

# find markers
obj <- readRDS("sketch_v5_all/it4_union_myel_sketch1000fullharmonycluster_findvar_pca_harmonyumap_res0.6.rds")
obj <- JoinLayers(obj)
obj[["RNA"]]$data <- as(object = obj[["RNA"]]$data, Class = "dgCMatrix")
markers <- FindAllMarkers(obj, max.cells.per.ident = 1000, only.pos = T)
write.csv(markers, "sketch_v5_all/it4_union_myel_sketch1000fullharmonycluster_findvar_pca_harmonyumap_res0.6_markers.csv")
markers$rank <- markers$pct.1 / (markers$pct.2 + 0.001) * markers$avg_log2FC
topmarkers <- markers %>%
  group_by(cluster) %>%
  slice_max(n = 100, order_by = rank)
write.csv(topmarkers, "sketch_v5_all/it4_union_myel_sketch1000fullharmonycluster_findvar_pca_harmonyumap_res0.6_topmarkers.csv")


obj <- readRDS("sketch_v5_all/it4_union_myel_sketch1000fullharmonycluster_findvar_pca_harmonyumap_res0.6.rds")
obj <- JoinLayers(obj)

# manual annotations for entire dataset
obj$manual.annotation <- "tbd"
obj$manual.annotation[obj$seurat_clusters == 0] <- "Macrophage C1qc+"
obj$manual.annotation[obj$seurat_clusters == 1] <- "Neutrophil Csf3r+"
obj$manual.annotation[obj$seurat_clusters == 2] <- "Neutrophil Csf3r+"
obj$manual.annotation[obj$seurat_clusters == 3] <- "Macrophage Spp1+"
obj$manual.annotation[obj$seurat_clusters == 4] <- "Macrophage Lyve1+"
obj$manual.annotation[obj$seurat_clusters == 5] <- "Neutrophil Csf3r+"
obj$manual.annotation[obj$seurat_clusters == 6] <- "Monocyte Ly6c2+"
obj$manual.annotation[obj$seurat_clusters == 7] <- "DC cDC2 Cd209a+"
obj$manual.annotation[obj$seurat_clusters == 8] <- "Monocyte Isg15+"
obj$manual.annotation[obj$seurat_clusters == 9] <- "DC mregDC Ccr7+"
obj$manual.annotation[obj$seurat_clusters == 10] <- "DC cDC1 Clec9a+"
obj$manual.annotation[obj$seurat_clusters == 11] <- "Monocyte Ly6c2+"
obj$manual.annotation[obj$seurat_clusters == 12] <- "Macrophage Marco+"
obj$manual.annotation[obj$seurat_clusters == 13] <- "Neutrophil Csf3r+"
obj$manual.annotation[obj$seurat_clusters == 14] <- "Neutrophil Csf3r+"
obj$manual.annotation[obj$seurat_clusters == 15] <- "DC pDC Tcf4+"
obj$manual.annotation[obj$seurat_clusters == 16] <- "Osteoclast Ctsk+"
obj$manual.annotation[obj$seurat_clusters == 17] <- "Kupffer Clec4f+"
obj$manual.annotation[obj$seurat_clusters == 18] <- "Macrophage Marco+"
obj$manual.annotation[obj$seurat_clusters == 19] <- "DC cDC2 Cd209a+"


saveRDS(obj, "sketch_v5_all/final_figures/mouse_union_myel_it4_res0.6_joined_annotated.rds", compress = F)

obj[["RNA"]]$counts <- as(obj[["RNA"]]$counts, Class = "dgCMatrix" )
obj[["RNA"]]$data <- as(obj[["RNA"]]$data, Class = "dgCMatrix" )
obj[["RNA"]]$scale.data <- as(obj[["RNA"]]$scale.data, Class = "dgCMatrix" )

saveRDS(obj, "sketch_v5_all/final_figures/mouse_union_myel_it4_res0.6_joined_annotated_sparse.rds", compress = F)

### Nerve iteration 1 -----
obj <- readRDS("sketch_v5_all/unfiltered_counts_obj_findvar_sketch1000_harmony_umap_res0.08_project.rds")
mnn <- readRDS("sketch_v5_all/unfiltered_counts_obj_findvar_sketch1000_mnn_umap_res0.08_clusterList.rds")
scvi <- readRDS("sketch_v5_all/unfiltered_counts_obj_findvar_sketch1000_scvi_umap_res0.08_clusterList.rds")
obj$mnn.cluster.full <- mnn
obj$scvi.cluster.full <- scvi

obj$nerve <- "false"
obj$nerve[obj$harmony.cluster.full %in% c(9)] <- "true"
obj$nerve[obj$mnn.cluster.full %in% c()] <- "true"
obj$nerve[obj$scvi.cluster.full %in% c(12)] <- "true"

obj$nerve[obj$cancer %in% c("melanoma", "MPNST", "neuroblastoma", "SCLC")] <- "false"

DefaultAssay(obj) <- "RNA"
obj <- JoinLayers(obj)
obj <- subset(obj, subset = nerve == "true")
saveRDS(obj, "sketch_v5_all/unfiltered_counts_obj_findvar_sketch1000_harmony_umap_res0.08_project_nervecutoff.rds", compress = F)

# re-generate seurat object and remove extraneous data
options(Seurat.object.assay.version = "v5")
obj <- CreateSeuratObject(counts = obj[["RNA"]]$counts, meta.data = obj@meta.data)
obj_backup <- obj
obj@meta.data <- obj@meta.data %>%
  select(any_of(c("orig.ident", "nCount_RNA","nFeature_RNA","study",
                  "organ","cancer","model","sorting","site","percent.mt","lin_neg")))
# re-perform BPcells
write_matrix_dir(mat = obj[["RNA"]]$counts, dir = "sketch_v5_all/union_nerve_sketch1000fullharmonycluster_counts")
counts.mat <- open_matrix_dir(dir = "sketch_v5_all/union_nerve_sketch1000fullharmonycluster_counts")
obj[["RNA"]]$counts <- counts.mat
saveRDS(obj, "sketch_v5_all/union_nerve_sketch1000fullharmonycluster_obj.rds",compress = F)

# run downstream analysis
obj <- NormalizeData(obj)
obj[["RNA"]] <- split(obj[["RNA"]], f = obj$study)
obj <- FindVariableFeatures(obj, verbose = T)
obj <- ScaleData(obj, verbose = F)
obj <- RunPCA(obj, npcs = 30, verbose = F)
obj <- IntegrateLayers(obj, method = HarmonyIntegration, orig.reduction = "pca", new.reduction = "harmony", kmeans_init_nstart=20, kmeans_init_iter_max=5000, verbose = F)
obj <- RunUMAP(obj, reduction = "harmony", dims = 1:30, reduction.name = "umap.harmony", min.dist = 0.001, seed.use = 10000)
saveRDS(obj, "sketch_v5_all/union_nerve_sketch1000fullharmonycluster_findvar_pca_harmonyumap.rds",compress = F)

# clustering
obj <- FindNeighbors(obj, reduction = "harmony", dims = 1:30)
#for (i in c(0.3, 0.1, 0.03, 0.01)) { 
for (i in c(0.09, 0.08, 0.07, 0.06, 0.05, 0.04)) { # choose 0.05
  obj <- FindClusters(obj, resolution = i)
  saveRDS(obj, paste0("sketch_v5_all/union_nerve_sketch1000fullharmonycluster_findvar_pca_harmonyumap_res",i,".rds"), compress = FALSE)
}

# find clusters
obj <- readRDS("sketch_v5_all/union_nerve_sketch1000fullharmonycluster_findvar_pca_harmonyumap_res0.05.rds")
obj <- JoinLayers(obj)
obj[["RNA"]]$data <- as(object = obj[["RNA"]]$data, Class = "dgCMatrix")
markers <- FindAllMarkers(obj, max.cells.per.ident = 1000, only.pos = T)
write.csv(markers, "sketch_v5_all/union_nerve_sketch1000fullharmonycluster_findvar_pca_harmonyumap_res0.05_markers.csv")
markers$rank <- markers$pct.1 / (markers$pct.2 + 0.001) * markers$avg_log2FC
topmarkers <- markers %>%
  group_by(cluster) %>%
  slice_max(n = 100, order_by = rank)
write.csv(topmarkers, "sketch_v5_all/union_nerve_sketch1000fullharmonycluster_findvar_pca_harmonyumap_res0.05_topmarkers.csv")

### Nerve iteration 2 -----
obj <- readRDS("sketch_v5_all/union_nerve_sketch1000fullharmonycluster_findvar_pca_harmonyumap_res0.05.rds")

obj$nerve <- "false"
obj$nerve[obj$seurat_clusters %in% c(1,3)] <- "true"

DefaultAssay(obj) <- "RNA"
obj <- JoinLayers(obj)
obj <- subset(obj, subset = nerve == "true")
saveRDS(obj, "sketch_v5_all/union_nerve_sketch1000fullharmonycluster_findvar_pca_harmonyumap_res0.05_nervecutoff.rds", compress = F)

# re-generate seurat object and remove extraneous data
options(Seurat.object.assay.version = "v5")
obj <- CreateSeuratObject(counts = obj[["RNA"]]$counts, meta.data = obj@meta.data)
obj_backup <- obj
obj@meta.data <- obj@meta.data %>%
  select(any_of(c("orig.ident", "nCount_RNA","nFeature_RNA","study",
                  "organ","cancer","model","sorting","site","percent.mt","lin_neg")))
# re-perform BPcells
write_matrix_dir(mat = obj[["RNA"]]$counts, dir = "sketch_v5_all/it2_union_nerve_sketch1000fullharmonycluster_counts")
counts.mat <- open_matrix_dir(dir = "sketch_v5_all/it2_union_nerve_sketch1000fullharmonycluster_counts")
obj[["RNA"]]$counts <- counts.mat
saveRDS(obj, "sketch_v5_all/it2_union_nerve_sketch1000fullharmonycluster_obj.rds",compress = F)

# run downstream analysis
obj <- NormalizeData(obj)
obj[["RNA"]] <- split(obj[["RNA"]], f = obj$study)
obj <- FindVariableFeatures(obj, verbose = T)
obj <- ScaleData(obj, verbose = F)
obj <- RunPCA(obj, npcs = 30, verbose = F)
obj <- IntegrateLayers(obj, method = HarmonyIntegration, orig.reduction = "pca", new.reduction = "harmony", kmeans_init_nstart=20, kmeans_init_iter_max=5000, verbose = F)
obj <- RunUMAP(obj, reduction = "harmony", dims = 1:30, reduction.name = "umap.harmony", min.dist = 0.001, seed.use = 10000)
saveRDS(obj, "sketch_v5_all/it2_union_nerve_sketch1000fullharmonycluster_findvar_pca_harmonyumap.rds",compress = F)

# clustering
obj <- FindNeighbors(obj, reduction = "harmony", dims = 1:30)
for (i in c(0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1, 0.12, 0.14, 0.16, 0.18, 0.2, 0.3)) {  # choose 0.03
  obj <- FindClusters(obj, resolution = i)
  DimPlot(obj, raster = TRUE, label = T)
  ggsave(paste0("sketch_v5_all/it2_union_nerve_sketch1000fullharmonycluster_findvar_pca_harmonyumap_res",i,".jpg"), width = 5, height = 5, units = "in", limitsize = FALSE)
  saveRDS(obj, paste0("sketch_v5_all/it2_union_nerve_sketch1000fullharmonycluster_findvar_pca_harmonyumap_res",i,".rds"), compress = FALSE)
}

# find markers
obj <- readRDS("sketch_v5_all/it2_union_nerve_sketch1000fullharmonycluster_findvar_pca_harmonyumap_res0.03.rds")
obj <- JoinLayers(obj)
obj[["RNA"]]$data <- as(object = obj[["RNA"]]$data, Class = "dgCMatrix")
markers <- FindAllMarkers(obj, max.cells.per.ident = 1000, only.pos = T)
write.csv(markers, "sketch_v5_all/it2_union_nerve_sketch1000fullharmonycluster_findvar_pca_harmonyumap_res0.03_markers.csv")
markers$rank <- markers$pct.1 / (markers$pct.2 + 0.001) * markers$avg_log2FC
topmarkers <- markers %>%
  group_by(cluster) %>%
  slice_max(n = 100, order_by = rank)
write.csv(topmarkers, "sketch_v5_all/it2_union_nerve_sketch1000fullharmonycluster_findvar_pca_harmonyumap_res0.03_topmarkers.csv")

# res 0.03: 0 = non-myel Schwann, 1 = myel Schwann

obj <- readRDS("sketch_v5_all/it2_union_nerve_sketch1000fullharmonycluster_findvar_pca_harmonyumap_res0.03.rds")
obj <- JoinLayers(obj)

# manual annotations for entire dataset
obj$manual.annotation <- "tbd"
obj$manual.annotation[obj$seurat_clusters == 0] <- "Schwann non-myel Ngfr+"
obj$manual.annotation[obj$seurat_clusters == 1] <- "Schwann myel Mpz+"

saveRDS(obj, "sketch_v5_all/final_figures/it2_union_nerve_sketch1000fullharmonycluster_findvar_pca_harmonyumap_res0.03_annotated.rds", compress = F)

obj[["RNA"]]$counts <- as(obj[["RNA"]]$counts, Class = "dgCMatrix" )
obj[["RNA"]]$data <- as(obj[["RNA"]]$data, Class = "dgCMatrix" )
obj[["RNA"]]$scale.data <- as(obj[["RNA"]]$scale.data, Class = "dgCMatrix" )

saveRDS(obj, "sketch_v5_all/final_figures/it2_union_nerve_sketch1000fullharmonycluster_findvar_pca_harmonyumap_res0.03_annotated_sparse.rds", compress = F)


