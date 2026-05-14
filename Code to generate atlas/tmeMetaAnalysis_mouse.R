library(Seurat); library(ggplot2); library(stringr); library(grid); library(gridExtra); library( dplyr );library(BPCells);library(Matrix)
source('tmeFunctions.R')
options(Seurat.object.assay.version = "v3")



# read in files
obj_list <- list.files(paste0("readyForSeurat/"),"*object.rds", recursive = TRUE)
obj_list <- substring(obj_list,0,nchar(obj_list)-11)
write.csv(obj_list, "obj_list_01_preQC.csv")
pb <- txtProgressBar(min = 1, max = length(obj_list), style = 3)
count <- 0
for (obj in obj_list) {
  count <- count + 1
  setTxtProgressBar(pb, count)
  eval(call("<-", as.name(obj), readRDS(paste0("readyForSeurat/",obj,".object.rds"))))
}

# confirm col = cells; row = features
out <- data.frame(matrix(ncol=4,nrow=0))
for ( i in 1:length(obj_list) ) {; print(i)
  obj = obj_list[i]
  x <- get(obj)
x[['RNA']] = as(x[['RNA']], Class = 'Assay')
  out[nrow(out)+1,] <- c(obj,x@assays[["RNA"]]@counts@Dimnames[[1]][1],x@assays[["RNA"]]@counts@Dimnames[[1]][10000],x@assays[["RNA"]]@counts@Dimnames[[2]][1]) 
}
write.csv(out, "head.csv")

# confirm data has no non-integer values
out <- data.frame(matrix(ncol=2,nrow=0))
for ( i in 1:length(obj_list) ) {; print(i)
  obj = obj_list[i]
  x <- get(obj)
  out[nrow(out)+1,] <- c(obj,any(GetAssayData(object = x, slot = "counts")%%1!=0)) # tests if integer
}
write.csv(out, "max_counts.csv")

# get assay name
sink(file = "assays.txt")
for (obj in obj_list) {
  x <- get(obj)
  print(obj)
  show(x)
}
sink(file = NULL)

# QC 
for ( i in 1:length(obj_list) ) {; print(i)
  objName = obj_list[i]
  setTxtProgressBar(pb, i)
  obj <- get(objName)
  DefaultAssay(obj) = 'RNA'
  if( !( 'nFeature_RNA' %in% colnames( obj@meta.data) ) ) {
    counts = obj[['RNA']]$counts
    obj$nCount_RNA = colSums( counts )
    obj$nFeature_RNA = colSums( counts > 0 )
  }
  if( max(obj$nFeature_RNA) <= 200 ) {; next; }
  obj <- subset(obj, subset = nFeature_RNA > 200)
  if( dim(obj)[2] <= 1 ) {; next; }
  obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^mt-")
  VlnPlot(obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
  ggsave(paste0("01_qc/", gsub("/","_",objName), "_qc-01-pre-filt.jpg"))
  obj <- subset(obj, subset = nFeature_RNA > 200 & percent.mt < 20)
  if( dim(obj)[2] <= 1 ) {; next; }
  VlnPlot(obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
  ggsave(paste0("01_qc/", gsub("/","_",objName), "_qc-02-post-filt.jpg"))
  eval(call("<-", as.name(objName), obj))
  dir.create(paste0("processed/",strsplit(objName,"/")[[1]][1]),suppressWarnings(F))
  saveRDS(obj, paste0("processed/",objName, "_01_qc.rds"), compress = F)
}

# output number of cells in file
out <- data.frame(matrix(ncol=2,nrow=0))
for (obj in obj_list) {
  x <- get(obj)
  out[nrow(out)+1,] <- c(obj,dim(x)[2])
}
write.csv(out, "cell_number_01_qc.csv")



obj_list <- as.vector(read.csv("obj_list_01_preQC.csv", colClasses=c("NULL",NA))$x)
pb <- txtProgressBar(min = 1, max = length(obj_list), style = 3)
count <- 0

## format individual post-QC objects 
for (i in 1:length(obj_list)) {; print(i)
  count <- count + 1
  obj <- obj_list[i]
  setTxtProgressBar(pb, count)
  x <- readRDS(paste0("processed/",obj,"_01_qc.rds"))
  if( dim(x)[2] <= 2 ) {; next; };
  x <- DietSeurat(x)
  x@meta.data <- x@meta.data %>% select(any_of(c("orig.ident", "nCount_RNA","nFeature_RNA","study", 
        "POD","size","age", "sequencingmethod", "treated", "model","region", "sorting","site","percent.mt","type","strain")))
  x <- RenameCells(x, new.names = paste(x$orig.ident, Cells(x), sep = "_"))
  # deal with duplicate gene names and non-standard gene names
  counts.mat <- GetAssayData( x, slot = 'counts' )
  geneNames <- rownames(counts.mat)
  synToName <- read.csv("GRCm39_geneNames_geneSynonyms.csv", header = T)
  for (i in 1:length(geneNames)) {
    name <- geneNames[i]
    if (name %in% synToName$Gene.name) {
      next
    } else if (name %in% synToName$Gene.Synonym){
      geneNames[i] <- synToName$Gene.name[which(synToName$Gene.Synonym == name)]; }; }
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
  saveRDS(x, paste0("processed/",obj, "_01b_correctGeneName.rds"), compress = F)
}


#### Export sample files to python format

library (zellkonverter); library(SingleCellExperiment)
obj_list <- list.files(paste0("processed/"),"*_01b_correctGeneName.rds", recursive = T); 
for( i in 1:length(obj_list) ) {; print(i)
  obj = readRDS( paste0( "processed/", obj_list[i] ) )
  sce = as.SingleCellExperiment( obj )
  writeH5AD (sce, paste0("sketch_v5_all/H5AD/forScrublet_",unique(obj$study),'_',unique(obj$orig.ident),'.h5ad' ) )
}


#### Run scrublet in python

conda activate scib/anndata2ri
python
import scanpy as sc
import scib
import anndata2ri
import anndata
import rpy2
from rpy2.robjects import r
from rpy2.robjects.conversion import localconverter
anndata2ri.activate()
import anndata
import os
import scrublet as scr
os.chdir('sketch_v5_all/H5AD')
dbList = []
fileList = os.listdir()
badCount = 0

for i in range( 1, len(fileList) ):
  file = fileList[i]
  if file.endswith(".h5ad"):
      print(i)
      print(i)
      obj = anndata.read_h5ad(file)
      if( obj.n_obs <= 30 ):
        badCount = badCount+1
        print('bad')
        print(badCount)
      else:
        print(file)
        scrub = scr.Scrublet(obj.X)
        doublet_scores, predicted_doublets = scrub.scrub_doublets()
        obj.obs['doublet_scores'] = doublet_scores
        obj.obs['predicted_doublets'] = predicted_doublets
        dbList.append( obj.obs )

import pickle
with open('scrublets.pkl', 'wb' ) as file:
  pickle.dump(dbList, file)


#### Filter for doublets using scrublet results

require("reticulate")
source_python("pickle_reader.py")
pickle_data <- read_pickle_file("scrublets.pkl")
obj_list <- list.files(paste0("processed/"),"*_01b_correctGeneName.rds", recursive = TRUE); obj_list <- substring(obj_list,0,nchar(obj_list)-24)

temp = c()
for( i in 1:length( pickle_data) ) {
   temp = c(temp, paste0( as.character(unique( pickle_data[[i]]$study )), '_', as.character(unique( pickle_data[[i]]$orig.ident )) ) )
}

for( i in 1:length(obj_list) ) {; print(i)
for( i in start:end ) {; print(i)
  objName <- obj_list[i]
  obj <- readRDS(paste0("processed/",objName,"_01b_correctGeneName.rds"));
  denom = dim(obj)[2]

  j = which( paste0( unique(obj$study), '_', unique(obj$orig.ident) )  == temp  )
  obj$scrublet_doublet_scores = pickle_data[[j]]$doublet_scores
  obj$scrublet_predicted_doublets = pickle_data[[j]]$predicted_doublets
  obj = subset( obj, scrublet_predicted_doublets==F )
  num   = dim(obj)[2]
  print( c( 100*(denom-num)/denom, ncol(obj) ) )
  saveRDS(obj, paste0("processed/",objName, "_01c_singlets.rds"), compress = F)
}


#### Standardize rows into common gene space

library(dplyr)
obj_list <- list.files(paste0("processed/"), "*_01c_singlets.rds", recursive = T); obj_list <- substring(obj_list,0,nchar(obj_list)-17)
synToName <- read.csv("GRCm39_geneNames_geneSynonyms.csv", header = T)
keepGenes <- sort(unique(synToName$Gene.name))
keepGenes <- keepGenes[-1] # remove empty "" gene ; length = 40901
for (i in 1:length(obj_list)) {; print(i)
  obj <- obj_list[i]
  x <- readRDS(paste0("processed/",obj,"_01c_singlets.rds"))
  if( ncol(x) <= 1 ) {; next; };
  counts.mat <- GetAssayData( x, slot = 'counts' )
  counts.mat <- convert_matrix_type(counts.mat, type = "uint32_t" )
  counts.mat <- counts.mat[rownames(counts.mat) %in% synToName$Gene.name,]
  counts.mat <- as(counts.mat, Class = "dgCMatrix")
  mat1 <- counts.mat
  ii1 = match( keepGenes, rownames( mat1 ) )
  ii1[ is.na(ii1) ] = dim(mat1)[1]+1
  mat1 = rbind( mat1, rep(0,ncol(mat1)) )
  new1 = mat1[ ii1, ]
  rownames(new1) <- keepGenes
  x <- CreateSeuratObject(counts = new1, meta.data = x@meta.data)
  dir.create(paste0("processed/",strsplit(obj,"/")[[1]][1]),suppressWarnings(FALSE))
  saveRDS(x, paste0("processed/",obj, "_01d_stddim.rds"), compress = FALSE)
}


#### Merge samples for each study

study_list <- list.files(paste0("processed/"), recursive = F)
start <- 1; end <- length(study_list)
for (i in start:end) {; print(i)
  study <- study_list[i]
  study_list_samples <- list.files(paste0("processed/",study), "*_01d_stddim.rds")
  if( length(study_list_samples) == 0 ) {; next; };
  study_list_obj <- list()
  for (sample in study_list_samples) {
    x <- readRDS(paste("processed",study,sample, sep = "/"))
    print(paste(Layers(x), collapse = " "))
    study_list_obj <- append(study_list_obj, x)
  }
  if (length(study_list_obj) > 1 ) {
    obj <- merge(study_list_obj[[1]], study_list_obj[2:length(study_list_obj)])
    if( ("data" %in% Layers(obj)) & length(Layers(obj))==2  ) {
      print('problem')
    } else {
      obj = JoinLayers( obj )
    }
  } else {
    obj <- study_list_obj[[1]]
  }
  if (dim(obj)[2] > 1) {
    saveRDS(obj, paste0("processed/",study,"/",study,"_01e_study.rds"), compress = F)
  }
}


#### Find genes with at least minimal expression across datasets

obj_list <- list.files(paste0("processed/"), "*_01e_study.rds", recursive = T)
cell_count <- 0

obj <- readRDS( paste0("processed/",obj_list[1]) )
gene_counts <- as.data.frame( matrix(0, ncol = ncol(obj), nrow = 0) )
colnames(gene_counts) <- ncol(obj)

for( i in 1:length(obj_list) ) {; print(i)
  objName <- obj_list[i]
  if( !file.exists( paste0("processed/",objName) ) ) {; next; };
  obj <- readRDS( paste0("processed/",objName) )
  counts.mat <- obj[["RNA"]]$counts
  cell_count <- cell_count + dim(counts.mat)[2]
  thing = rowSums( counts.mat > 0)
  gene_counts <- rbind(gene_counts, thing )
}
nCount_Feature <- colSums( gene_counts, na.rm=T )
keepGenesIdx <- c(nCount_Feature >=0.001*cell_count)
keepGenesNames <- colnames(gene_counts)[keepGenesIdx]
write.csv(keepGenesIdx, "keepGenesIdx.csv", quote = F)
write.csv(keepGenesNames, "keepGenesNames.csv", quote = F, row.names = F)


#### Place all files in BP cell space

obj_list <- list.files(paste0("processed/"), "*_01e_study.rds", recursive = T)
start <- 1; end <- length(obj_list)
for (i in start:end) {; print(i)
  obj <- readRDS(paste0("processed/",obj_list[i]))
  counts.mat <- obj[["RNA"]]$counts
  counts.mat <- counts.mat[keepGenesIdx, ]
  # re-perform BP cells 
  write_matrix_dir(mat = counts.mat, dir = paste0("BP/",substring(obj_list[i], 1, nchar(obj_list[i])-14)))
  counts.mat <- open_matrix_dir(dir = paste0("BP/",substring(obj_list[i], 1, nchar(obj_list[i])-14)))
  obj <- CreateSeuratObject(counts = counts.mat, meta.data = obj@meta.data)
  dir.create(paste0("BP_object/",strsplit(obj_list[i],"/")[[1]][1]),suppressWarnings(FALSE))
  saveRDS(obj, paste0("BP_object/",substring(obj_list[i], 1, nchar(obj_list[i])-14),"_01f_prunedGenes.rds"), compress = F)
}


#### Merge as a new layer in seurat v5

options(Seurat.object.assay.version = "v5")
obj_list <- list.files(paste0("BP_object/"), "*_01f_prunedGenes.rds", recursive = T)
split_list <- strsplit(obj_list, "/")
first_elements <- sapply(split_list, function(x) x[1])
caf_list <- list()
cell_count = 0
for (i in 1:length(obj_list)) {; print(i)
  obj <- readRDS( paste0( "BP_object/",obj_list[i] ) ); print(first_elements[i])
  obj$study = first_elements[i]
  caf_list[[i]] <- obj
  cell_count <- cell_count + dim(obj)[2]
}

cell_count
merged <- merge(caf_list[[1]], caf_list[2:length(caf_list)])
saveRDS(merged, "sketch_v5_all/mouseMerged.rds", compress = F)


#### Attach meta-data
obj = merged
clinical <- read.csv("mouse.meta.data.csv", header = T)
clinical <- clinical[,-1]
mapping <- match(obj$orig.ident, clinical$orig.ident)
for (col in colnames(clinical)) {
  print(col)
  obj[[col]] <- clinical[[col]][mapping]
}


#### Integrate with Harmony

obj <- NormalizeData(obj)
saveRDS(obj, "sketch_v5_all/unfiltered_counts_obj_norm.rds", compress = F)
#obj[["RNA"]] <- split(obj[["RNA"]], f = obj$study)
obj <- FindVariableFeatures(obj, verbose = T)
saveRDS(obj, "sketch_v5_all/unfiltered_counts_obj_findvar.rds", compress = F)
obj <- SketchData(obj, ncells = 1000, method = "LeverageScore", sketched.assay = "sketch")
saveRDS(obj, "sketch_v5_all/unfiltered_counts_obj_findvar_sketch1000.rds", compress = F)
DefaultAssay(obj) <- "sketch"
obj <- FindVariableFeatures(obj, verbose = T)
saveRDS(obj, "sketch_v5_all/unfiltered_counts_obj_findvar_sketch1000_findvar.rds", compress = F)
obj <- ScaleData(obj, verbose = F)
saveRDS(obj, "sketch_v5_all/unfiltered_counts_obj_findvar_sketch1000_scale.rds", compress = F)
obj <- RunPCA(obj, npcs = 30, verbose = T)
saveRDS(obj, "sketch_v5_all/unfiltered_counts_obj_findvar_sketch1000_pca.rds", compress = F)

obj <- IntegrateLayers(obj, method = HarmonyIntegration, orig.reduction = "pca", new.reduction = "harmony", kmeans_init_nstart=20, kmeans_init_iter_max=5000, verbose = T)
saveRDS(obj, "sketch_v5_all/unfiltered_counts_obj_findvar_sketch1000_harmony.rds", compress = F)
obj <- RunUMAP(obj, reduction = "harmony", dims = 1:30, reduction.name = "umap.harmony", return.model = T, min.dist = 0.001)
saveRDS(obj, "sketch_v5_all/unfiltered_counts_obj_findvar_sketch1000_harmony_umap.rds",compress = F)
# clustering
obj <- FindNeighbors(obj, reduction = "harmony", dims = 1:30)
for (i in c(0.06, 0.07, 0.08, 0.09, 0.1, 0.12, 0.14, 0.16, 0.018, 0.2)) {
  obj <- FindClusters(obj, resolution = i, graph.name = 'sketch_snn')
  DimPlot(obj, raster = T, label=T)
  ggsave(paste0("sketch_v5_all/unfiltered_counts_obj_findvar_sketch1000YES_harmony_umap_res",i,".jpg"), width = 5, height = 5, units = "in", limitsize = FALSE)
  saveRDS(obj, paste0("sketch_v5_all/unfiltered_counts_obj_findvar_sketch1000YES_harmony_umap_res",i,".rds"), compress = FALSE)
}

# Project to full dataset
obj$seurat_clusters = obj$sketch_snn_res.0.2
obj <- ProjectIntegration(object = obj, reduction = "harmony")
options(future.globals.maxSize = 8000 * 1024^2)
obj <- ProjectData(object = obj, sketched.reduction = "harmony.full", full.reduction = "harmony.full", umap.model = "umap.harmony", dims = 1:30, refdata = list(harmony.cluster.full = "seurat_clusters"))
saveRDS(obj, "sketch_v5_all/unfiltered_counts_obj_findvar_sketch1000_harmony_umap_project.rds",compress = F)



#### Integrate with MNN

obj <- readRDS("sketch_v5_all/unfiltered_counts_obj_findvar_sketch1000_pca.rds")
library(SeuratWrappers)
obj <- IntegrateLayers(
  object = obj, method = FastMNNIntegration,
  orig.reduction = "pca", new.reduction = "integrated.mnn", verbose = T
)
saveRDS(obj, "sketch_v5_all/unfiltered_counts_obj_findvar_sketch1000_mnn.rds", compress = F)
obj <- RunUMAP(obj, reduction = "integrated.mnn", dims = 1:30, reduction.name = "umap.mnn", return.model = T, min.dist = 0.001)
saveRDS(obj, "sketch_v5_all/unfiltered_counts_obj_findvar_sketch1000_mnn_umap.rds",compress = F)

obj <- FindNeighbors(obj, reduction = "integrated.mnn", dims = 1:30)
for (i in c(0.05,0.06,0.07,0.09,1.0,1.5,2.0) ) {
  obj <- FindClusters(obj, resolution = i, graph.name = 'RNA_snn')
  DimPlot(obj, raster = TRUE, label=T, reduction =  'umap.mnn')
  ggsave(paste0("sketch_v5_all/unfiltered_counts_obj_findvar_sketch1000_mnn_umap_res",i,".jpg"), width = 5, height = 5, units = "in", limitsize = FALSE)
  saveRDS(obj, paste0("sketch_v5_all/unfiltered_counts_obj_findvar_sketch1000_mnn_umap_res",i,".rds"), compress = FALSE)
}

# Project to full dataset
obj$seurat_clusters <- obj$RNA_snn_res.0.15
obj <- ProjectIntegration(object = obj, reduction = "integrated.mnn")
options(future.globals.maxSize = 8000 * 1024^2)
obj <- ProjectData(object = obj, sketched.reduction = "integrated.mnn.full", full.reduction = "integrated.mnn.full", umap.model = "umap.mnn", dims = 1:30, refdata = list(mnn.cluster.full = "seurat_clusters"))
saveRDS(obj, "sketch_v5_all/unfiltered_counts_obj_findvar_sketch1000_mnn_umap_project.rds",compress = F)



#### Integrate with scVI

Sys.setenv( RETICULATE_PYTHON = "scvi-env/bin/python" )
library(reticulate)
obj <- readRDS("sketch_v5_all/unfiltered_counts_obj_findvar_sketch1000_pca.rds")
obj <- IntegrateLayers(
  object = obj, method = scVIIntegration,
  orig.reduction = "pca", new.reduction = "integrated.scvi", conda_env = '~/.conda/envs/scvi-env/', verbose = T
)
saveRDS(obj, "sketch_v5_all/unfiltered_counts_obj_findvar_sketch1000_scvi.rds", compress = F)
obj <- RunUMAP(obj, reduction = "integrated.scvi", dims = 1:30, reduction.name = "umap.scvi", return.model = T, min.dist = 0.001)
ggsave("sketch_v5_all/unfiltered_obj_sketch1000_scvi_featureplot.jpg", width = 8, height = 10, units = "in", limitsize = FALSE)
saveRDS(obj, "sketch_v5_all/unfiltered_counts_obj_findvar_sketch1000_scvi_umap.rds",compress = F)

obj <- FindNeighbors(obj, reduction = "integrated.scvi", dims = 1:30)
for (i in c(0.05,0.06,0.07,0.09,1.0,1.5,2.0) ) {
  obj <- FindClusters(obj, resolution = i, graph.name = 'RNA_snn')
  DimPlot(obj, raster = TRUE, reduction =  'umap.scvi', label=T)
  ggsave(paste0("sketch_v5_all/unfiltered_counts_obj_findvar_sketch1000_scvi_umap_res",i,".jpg"), width = 5, height = 5, units = "in", limitsize = FALSE)
  saveRDS(obj, paste0("sketch_v5_all/unfiltered_counts_obj_findvar_sketch1000_scvi_umap_res",i,".rds"), compress = FALSE)
}

# Project to full dataset
obj$seurat_clusters <- obj$RNA_snn_res.0.4
obj <- ProjectIntegration(object = obj, reduction = "integrated.scvi")
options(future.globals.maxSize = 8000 * 1024^2)
obj <- ProjectData(object = obj, sketched.reduction = "integrated.scvi.full", full.reduction = "integrated.scvi.full", umap.model = "umap.scvi", dims = 1:30, refdata = list(scvi.cluster.full = "seurat_clusters"))
saveRDS(obj, "sketch_v5_all/unfiltered_counts_obj_findvar_sketch1000_scvi_umap_res0.4_project.rds",compress = F)







# Manual annotations for entire dataset (Harmony)
temp = as.character( obj$harmony.cluster.full )
temp[ temp == 0 ] = 'Epithelial'
temp[ temp == 1 ] = 'Fibroblast'
temp[ temp == 2 ] = 'Myeloid'
temp[ temp == 3 ] = 'Lymphoid'
temp[ temp == 4 ] = 'Endothelial'
temp[ temp == 5 ] = 'Cycling'
temp[ temp == 6 ] = 'Myeloid'
temp[ temp == 7 ] = 'Lymphoid'
temp[ temp == 8 ] = 'Fibroblast'
temp[ temp == 9 ] = 'Mural'
temp[ temp ==10 ] = 'Fibroblast'
temp[ temp ==11 ] = 'Epithelial'
temp[ temp ==12 ] = 'Epithelial' 
temp[ temp ==13 ] = 'Nerve'
temp[ temp ==14 ] = 'Epithelial'
temp[ temp ==15 ] = 'Lymphoid'
temp[ temp ==16 ] = 'Lymphoid'
obj$cellType = temp

# Manual annotations for entire dataset (MNN)
temp = as.character( obj$mnn.cluster.full )
temp[ temp == 0 ] = 'Epithelial'
temp[ temp == 1 ] = 'Fibroblast'
temp[ temp == 2 ] = 'Myeloid'
temp[ temp == 3 ] = 'Lymphoid'
temp[ temp == 4 ] = 'Endothelial'
temp[ temp == 5 ] = 'Myeloid'
temp[ temp == 6 ] = 'Lymphoid'
temp[ temp == 7 ] = 'Epithelial'
temp[ temp == 8 ] = 'Epithelial'
temp[ temp == 9 ] = 'Mural'
temp[ temp ==10 ] = 'Fibroblast'
temp[ temp ==11 ] = 'Epithelial'
temp[ temp ==12 ] = 'Epithelial'
temp[ temp ==13 ] = 'Epithelial' 
temp[ temp ==14 ] = 'Epithelial'
temp[ obj$RNA_snn_res.4 == 42 ] = 'Nerve'
obj$cellType = temp

# Manual annotations for entire dataset (scVI)
temp = as.character( obj$scvi.cluster.full )
temp[ temp == 0 ] = 'Fibroblast'
temp[ temp == 1 ] = 'Myeloid'
temp[ temp == 2 ] = 'Lymphoid'
temp[ temp == 3 ] = 'Fibroblast'
temp[ temp == 4 ] = 'Endothelial'
temp[ temp == 5 ] = 'Myeloid'
temp[ temp == 6 ] = 'Epithelial'
temp[ temp == 7 ] = 'Epithelial'
temp[ temp == 8 ] = 'Epithelial'
temp[ temp == 9 ] = 'Epithelial'
temp[ temp ==10 ] = 'Epithelial'
temp[ temp ==11 ] = 'Lymphoid'
temp[ temp ==12 ] = 'Myeloid'
temp[ temp ==13 ] = 'Mural'
temp[ temp ==14 ] = 'Nerve'
temp[ temp ==15 ] = 'Fibroblast'
temp[ temp ==16 ] = 'Epithelial'
temp[ temp ==17 ] = 'Myeloid'
temp[ temp ==18 ] = 'Epithelial'
temp[ temp ==19 ] = 'Lymphoid'
temp[ temp ==20 ] = 'Epithelial'
temp[ temp ==21 ] = 'Endothelial'
temp[ temp ==22 ] = 'Epithelial'
temp[ temp ==23 ] = 'Nerve'
temp[ temp ==24 ] = 'Lymphoid'
temp[ temp ==25 ] = 'Epithelial'
temp[ temp ==26 ] = 'Mural' 
temp[ temp ==27 ] = 'Myeloid'
temp[ temp ==28 ] = 'Epithelial'
temp[ temp ==29 ] = 'Lymphoid'
temp[ temp ==30 ] = 'Myeloid'
temp[ temp ==31 ] = 'Fibroblast'
temp[ temp ==32 ] = 'Lymphoid'
temp[ temp ==33 ] = 'Myeloid'
temp[ temp ==34 ] = 'Epithelial'
temp[ temp ==35 ] = 'Epithelial' 
obj$cellType = temp





obj = combineIntegration('Fibroblast')
obj = analyzeObj( obj, 'caf_it1')
obj = combineIntegration('Endothelial')
obj = analyzeObj( obj, 'cec_it1')
obj = combineIntegration('Mural')
obj = analyzeObj( obj, 'cmc_it1')
obj = combineIntegration('Nerve')
obj = analyzeObj( obj, 'can_it1')
obj = combineIntegration('Lymphoid')
obj = analyzeObj( obj, 'cal_it1')
obj = combineIntegration('Myeloid')
obj = analyzeObj( obj, 'cam_it1')



###################### Fibroblast ###############################

# Fibroblast_it1
temp = as.character( obj$RNA_snn_res.0.09 )
temp[ temp == 0 ] = 'Fibroblast'
temp[ temp == 1 ] = 'Cycling'
temp[ temp == 2 ] = 'Epithelial'
temp[ temp == 3 ] = 'Myeloid'
temp[ temp == 4 ] = 'Fibroblast'
temp[ temp == 5 ] = 'Endothelial'
temp[ temp == 6 ] = 'Nerve'
temp[ temp == 7 ] = 'Fibroblast'
temp[ temp == 8 ] = 'Lymphoid'
temp[ temp == 9 ] = 'Low Quality'
temp[ temp ==10 ] = 'Fibroblast'
temp[ temp ==11 ] = 'Epithelial'
temp[ temp ==12 ] = 'Epithelial'
temp[ temp ==13 ] = 'Myeloid'
temp[ temp ==14 ] = 'Epithelial'
temp[ temp ==15 ] = 'Low Quality'
temp[ temp ==16 ] = 'Fibroblast'
temp[ temp ==17 ] = 'Low Quality'
obj$cellType = temp
saveRDS( obj, 'sketch_v5_all/caf_it1.rds', compress = F )
obj = subset( obj, cellType == 'Fibroblast' )
obj = analyzeObj( obj, 'caf_it2' )

# Fibroblast_it2
temp = as.character( obj$RNA_snn_res.0.06 )
temp[ temp == 0 ] = 'Fibroblast'
temp[ temp == 1 ] = 'Fibroblast'
temp[ temp == 2 ] = 'Fibroblast'
temp[ temp == 3 ] = 'Cycling'
temp[ temp == 4 ] = 'Lymphoid'
temp[ temp == 5 ] = 'Lymphoid'
obj$cellType = temp
saveRDS( obj, 'sketch_v5_all/caf_it2.rds', compress = F )
obj = subset( obj, cellType == 'Fibroblast' )
obj = analyzeObj( obj, 'caf_it3' )

# Fibroblast_it3
temp = as.character( obj$RNA_snn_res.0.18 )
temp[ temp == 0 ] = 'Fibroblast'
temp[ temp == 1 ] = 'Fibroblast'
temp[ temp == 2 ] = 'Fibroblast'
temp[ temp == 3 ] = 'Fibroblast'
temp[ temp == 4 ] = 'Fibroblast'
temp[ temp == 5 ] = 'Fibroblast'
temp[ temp == 6 ] = 'Epithelial'
temp[ temp == 7 ] = 'Low Quality'
temp[ temp == 8 ] = 'Epithelial'
temp[ temp == 9 ] = 'Fibroblast'
temp[ temp ==10 ] = 'Fibroblast'
temp[ temp ==11 ] = 'Fibroblast'
temp[ temp ==12 ] = 'Fibroblast'
temp[ temp ==13 ] = 'Low Quality'
obj$cellType = temp
saveRDS( obj, 'sketch_v5_all/caf_it3.rds', compress = F )
obj = subset( obj, cellType == 'Fibroblast' )
obj = analyzeObj( obj, 'caf_it4' )

# Fibroblast_it4
temp = as.character( obj$RNA_snn_res.0.2 )
temp[ temp == 0 ] = 'mCAF Lrrc15+'
temp[ temp == 1 ] = 'ssCAF Cxcl12+' 
temp[ temp == 2 ] = 'Low Quality'
temp[ temp == 3 ] = 'ssCAF Pi16+'
temp[ temp == 4 ] = 'apCAF Cd74+'
temp[ temp == 5 ] = 'iCAF Isg15+'
temp[ temp == 6 ] = 'ssCAF Pi16+'
temp[ temp == 7 ] = 'Pericyte'
temp[ temp == 8 ] = 'Low Quality'
temp[ temp == 9 ] = 'apCAF Cd74+'
temp[ temp ==10 ] = 'ssCAF Pi16+'
temp[ temp ==11 ] = 'mCAF Lrrc15+'
temp[ temp ==12 ] = 'ssCAF Cxcl12+'
temp[ temp ==13 ] = 'Epithelial'
obj$cellType = temp
saveRDS( obj, 'sketch_v5_all/caf_it4.rds', compress = F )
obj = subset( obj, cellType != 'Low Quality' & cellType != 'Pericyte' & cellType != 'Epithelial' )
obj = analyzeObj( obj, 'caf_it5' )

# Fibroblast_it5
temp = as.character( obj$RNA_snn_res.0.16 )
temp[ temp == 0 ] = 'ssCAF Cxcl12+'
temp[ temp == 1 ] = 'mCAF Lrrc15+'
temp[ temp == 2 ] = 'ssCAF Pi16+'
temp[ temp == 3 ] = 'iCAF Spp1+'
temp[ temp == 4 ] = 'apCAF Cd74+'
temp[ temp == 5 ] = 'iCAF Isg15+'
temp[ temp == 6 ] = 'ssCAF Cxcl12+'
temp[ temp == 7 ] = 'apCAF Cd74+'
temp[ temp == 8 ] = 'apCAF Cd74+'
temp[ temp == 9 ] = 'ssCAF Pi16+'
temp[ temp ==10 ] = 'mCAF Lrrc15+'
temp[ temp ==11 ] = 'apCAF Cd74+'
obj$cellType = temp
Idents(obj) = obj$cellType
saveRDS( obj, 'sketch_v5_all/caf_final.rds', compress = F )
makeFinalFigures( obj_caf, 'mouse_Fibroblast',c('Cd74','Isg15','Spp1','Lrrc15','Cxcl12','Pi16') )



###################### Endothelial ###############################

# Endothelial_it1
temp = as.character( obj$RNA_snn_res.0.1 )
temp[ temp == 0 ] = 'Endothelial'
temp[ temp == 1 ] = 'Low Quality'
temp[ temp == 2 ] = 'Fibroblast'
temp[ temp == 3 ] = 'Endothelial'
temp[ temp == 4 ] = 'Myeloid'
temp[ temp == 5 ] = 'Lymphoid'
temp[ temp == 6 ] = 'Epithelial'
temp[ temp == 7 ] = 'Endothelial'
temp[ temp == 8 ] = 'Lymphoid'
temp[ temp == 9 ] = 'Cycling'
temp[ temp ==10 ] = 'Myeloid'
obj$cellType = temp
saveRDS( obj, 'sketch_v5_all/cec_it1.rds', compress = F )
obj = subset( obj, cellType == 'Endothelial' )
obj = analyzeObj( obj, 'cec_it2' )

# Endothelial_it2
temp = as.character( obj$RNA_snn_res.0.1 )
temp[ temp == 0 ] = 'Endothelial'
temp[ temp == 1 ] = 'Endothelial'
temp[ temp == 2 ] = 'Endothelial'
temp[ temp == 3 ] = 'Endothelial'
temp[ temp == 4 ] = 'Lymphoid'
temp[ temp == 5 ] = 'Low Quality'
temp[ temp == 6 ] = 'Lymphoid'
obj$cellType = temp
saveRDS( obj, 'sketch_v5_all/cec_it2.rds', compress = F )
obj = subset( obj, cellType == 'Endothelial' )
obj = analyzeObj( obj, 'cec_it3' )

# Endothelial_it3
temp = as.character( obj$RNA_snn_res.0.14 )
temp[ temp == 0 ] = 'Capillary Car4+'
temp[ temp == 1 ] = 'Vein Ackr1+'
temp[ temp == 2 ] = 'Tip Esm1+'
temp[ temp == 3 ] = 'Lymphatic Prox1+'
temp[ temp == 4 ] = 'Artery Gja5+'
temp[ temp == 5 ] = 'Capillary Car4+'
obj$cellType = temp
Idents(obj) = obj$cellType
saveRDS( obj, 'sketch_v5_all/cec_final.rds', compress = F )
makeFinalFigures( obj_cec, 'mouse_Endothelial',c('Gja5','Car4','Prox1','Esm1','Ackr1') )



######################## Mural ###################################

# Mural_it1
temp = as.character( obj$RNA_snn_res.0.06 )
temp[ temp == 0 ] = 'Mural'
temp[ temp == 1 ] = 'Fibroblast'
temp[ temp == 2 ] = 'Low Quality'
temp[ temp == 3 ] = 'Fibroblast'
temp[ temp == 4 ] = 'Myeloid'
temp[ temp == 5 ] = 'Cycling'
temp[ temp == 6 ] = 'Endothelial'
temp[ temp == 7 ] = 'Fibroblast'
temp[ temp == 8 ] = 'Low Quality'
obj$cellType = temp
saveRDS( obj, 'sketch_v5_all/cmc_it1.rds', compress = F )
obj = subset( obj, cellType == 'Mural' )
obj = analyzeObj( obj, 'cmc_it2' )

# Mural_it2
temp = as.character( obj$RNA_snn_res.0.16 )
temp[ temp == 0 ] = 'Smooth muscle'
temp[ temp == 1 ] = 'Pericyte Cd248+'
temp[ temp == 2 ] = 'Pericyte Ccl2+'
temp[ temp == 3 ] = 'Pericyte Cd248+'
temp[ temp == 4 ] = 'Smooth muscle'
obj$cellType = temp
saveRDS( obj, 'sketch_v5_all/cmc_final.rds', compress = F )
makeFinalFigures( obj_cmc, 'mouse_Mural',c('Ccl2','Cd248','Myh11'))



######################## Nerve ###################################

# Nerve_it1
temp = as.character( obj$RNA_snn_res.0.16 )
temp[ temp == 0 ] = 'Low Quality'
temp[ temp == 1 ] = 'Nerve'
temp[ temp == 2 ] = 'Cycling'
temp[ temp == 3 ] = 'Cycling'
temp[ temp == 4 ] = 'Cycling'
temp[ temp == 5 ] = 'Epithelial'
temp[ temp == 6 ] = 'Epithelial'
temp[ temp == 7 ] = 'Myeloid'
temp[ temp == 8 ] = 'Myeloid'
temp[ temp == 9 ] = 'Epithelial'
temp[ temp ==10 ] = 'Nerve'
temp[ temp ==11 ] = 'Endothelial'
temp[ temp ==12 ] = 'Epithelial'
obj$cellType = temp
saveRDS( obj, 'sketch_v5_all/can_it1.rds', compress = F )
obj = subset( obj, cellType == 'Nerve' )
obj = analyzeObj( obj, 'can_it2' )

# Nerve_it2
temp = as.character( obj$RNA_snn_res.0.06 )
temp[ temp == 0 ] = 'Low Quality'
temp[ temp == 1 ] = 'Nerve'
temp[ temp == 2 ] = 'Nerve'
temp[ temp == 3 ] = 'Nerve'
temp[ temp == 4 ] = 'Cycling'
temp[ temp == 5 ] = 'Mural'
temp[ temp == 6 ] = 'Epithelial'
obj$cellType = temp
saveRDS( obj, 'sketch_v5_all/can_it2.rds', compress = F )
obj = subset( obj, cellType == 'Nerve' )
obj = analyzeObj( obj, 'can_it3' )

# Nerve_it3
temp = as.character( obj$RNA_snn_res.0.06 )
temp[ temp == 0 ] = 'Low Quality'
temp[ temp == 1 ] = 'Nerve'
temp[ temp == 2 ] = 'Nerve'
temp[ temp == 3 ] = 'Low Quality'
obj$cellType = temp
saveRDS( obj, 'sketch_v5_all/can_it3.rds', compress = F )
obj = subset( obj, cellType == 'Nerve' )
obj = analyzeObj( obj, 'can_it4' )

# Nerve_it4
temp = as.character( obj$RNA_snn_res.0.1 )
temp[ temp == 0 ] = 'Schwann non-myel Ngfr+'
temp[ temp == 1 ] = 'Schwann myel Mpz+'
temp[ temp == 2 ] = 'Schwann non-myel Ngfr+'
temp[ temp == 3 ] = 'Low Quality'
obj$cellType = temp
obj = subset(obj, cellType != 'Low Quality')
Idents(obj) = obj$cellType
saveRDS( obj, 'sketch_v5_all/can_final.rds', compress = F )
makeFinalFigures( obj_can, 'mouse_Nerve',c('Mpz','Ngfr') )



####################### Lymphoid #################################

# Lymphoid_it1
temp = as.character( obj$RNA_snn_res.0.06 )
temp[ temp == 0 ] = 'Lymphoid'
temp[ temp == 1 ] = 'Lymphoid'
temp[ temp == 2 ] = 'Fibroblast'
temp[ temp == 3 ] = 'Myeloid'
temp[ temp == 4 ] = 'Cycling'
temp[ temp == 5 ] = 'Lymphoid'
temp[ temp == 6 ] = 'Endothelial'
temp[ temp == 7 ] = 'Erythroid'
temp[ temp == 8 ] = 'Low Quality'
temp[ temp == 9 ] = 'Lymphoid'
temp[ temp ==10 ] = 'Low Quality'
obj$cellType = temp
saveRDS( obj, 'sketch_v5_all/cal_it1.rds', compress = F )
obj = subset( obj, cellType == 'Lymphoid' )
obj = analyzeObj( obj, 'cal_it2' )

# Lymphoid_it2
temp = as.character( obj$RNA_snn_res.0.6 )
temp[ temp == 0 ] = 'B Naive Ighm+'
temp[ temp == 1 ] = 'T Cd4+ Naive'
temp[ temp == 2 ] = 'T Cd8+ Naive'
temp[ temp == 3 ] = 'NK Klrk1+'
temp[ temp == 4 ] = 'B Naive Ighm+'
temp[ temp == 5 ] = 'Treg Foxp3+'
temp[ temp == 6 ] = 'Tex Pdcd1+'
temp[ temp == 7 ] = 'Treg Foxp3+'
temp[ temp == 8 ] = 'Th17 Il17a+'
temp[ temp == 9 ] = 'Plasma Jchain+'
temp[ temp ==10 ] = 'B Mem Bank1+'
temp[ temp ==11 ] = 'T Cd8+/Isg15+'
temp[ temp ==12 ] = 'B Naive Ighm+'
temp[ temp ==13 ] = 'B Naive Ighm+'
temp[ temp ==14 ] = 'B Mem Bank1+'
temp[ temp ==15 ] = 'Th2 Gata3+'
temp[ temp ==16 ] = 'Plasma Jchain+'
temp[ temp ==17 ] = 'Treg Foxp3+'
temp[ temp ==18 ] = 'T Cd8+ Naive'
temp[ temp ==19 ] = 'NK Klrk1+'
temp[ temp ==20 ] = 'Plasma Jchain+'
temp[ temp ==21 ] = 'T Cd8+ Naive'
temp[ temp ==22 ] = 'T Cd4+ Naive'
obj$cellType = temp
Idents(obj) = obj$cellType
saveRDS( obj, 'sketch_v5_all/cal_final.rds', compress = F )
makeFinalFigures( obj_cal, 'mouse_Lymphod',c('Bank1','Ighm','Klrk1','Jchain','Cd4','Cd8a','Isg15','Pdcd1','Il17a','Gata3','Foxp3') )



####################### Myeloid #################################

# Myeloid_it1
temp = as.character( obj$RNA_snn_res.0.14 )
temp[ temp == 0 ] = 'Myeloid'
temp[ temp == 1 ] = 'Myeloid'
temp[ temp == 2 ] = 'Low Quality'
temp[ temp == 3 ] = 'Epithelial'
temp[ temp == 4 ] = 'Myeloid'
temp[ temp == 5 ] = 'Endothelial'
temp[ temp == 6 ] = 'Cycling'
temp[ temp == 7 ] = 'Lymphoid'
temp[ temp == 8 ] = 'Myeloid'
temp[ temp == 9 ] = 'Lymphoid'
temp[ temp ==10 ] = 'Endothelial'
temp[ temp ==11 ] = 'Myeloid'
temp[ temp ==12 ] = 'Myeloid'
temp[ temp ==13 ] = 'Myeloid'
temp[ temp ==14 ] = 'Low Quality'
temp[ temp ==15 ] = 'Myeloid'
temp[ temp ==16 ] = 'Myeloid'
temp[ temp ==17 ] = 'Low Quality'
temp[ temp ==18 ] = 'Low Quality'
temp[ temp ==19 ] = 'Myeloid'
temp[ temp ==20 ] = 'Low Quality'
obj$cellType = temp
saveRDS( obj, 'sketch_v5_all/cam_it1.rds', compress = F )
obj = subset( obj, cellType == 'Myeloid' )
obj = analyzeObj( obj, 'cam_it2' )

# Myeloid_it2
temp = as.character( obj$RNA_snn_res.0.8 )
temp[ temp == 0 ] = 'Neutrophil Csf3r+'
temp[ temp == 1 ] = 'Macrophage C1qc+'
temp[ temp == 2 ] = 'Neutrophil Csf3r+'
temp[ temp == 3 ] = 'Macrophage Lyve1+'
temp[ temp == 4 ] = 'Monocyte Ly6c2+'
temp[ temp == 5 ] = 'Monocyte Isg15+'
temp[ temp == 6 ] = 'Macrophage Spp1+'
temp[ temp == 7 ] = 'Neutrophil Csf3r+'
temp[ temp == 8 ] = 'Macrophage Spp1+'
temp[ temp == 9 ] = 'DC cDC2 Cd209a+'
temp[ temp ==10 ] = 'Neutrophil Csf3r+'
temp[ temp ==11 ] = 'DC mregDC Ccr7+'
temp[ temp ==12 ] = 'Monocyte Ly6c2+'
temp[ temp ==13 ] = 'Monocyte Cd16+' # Ear2
temp[ temp ==14 ] = 'DC cDC1 Clec9a+'
temp[ temp ==15 ] = 'Neutrophil Csf3r+'
temp[ temp ==16 ] = 'DC pDC Tcf4+'
temp[ temp ==17 ] = 'DC cDC2 Cd209a+'
temp[ temp ==18 ] = 'Mast Cpa3+'
temp[ temp ==19 ] = 'Mast Cpa3+'
temp[ temp ==20 ] = 'Cycling'
temp[ temp ==21 ] = 'Monocyte Ly6c2+'
temp[ temp ==22 ] = 'Macrophage C1qc+'
temp[ temp ==23 ] = 'Macrophage C1qc+'
temp[ temp ==24 ] = 'Neutrophil Csf3r+'
temp[ temp ==25 ] = 'Osteoclast Ctsk+'
temp[ temp ==26 ] = 'Macrophage Spp1+'
temp[ temp ==27 ] = 'Macrophage Marco+'
temp[ temp ==28 ] = 'Neutrophil Csf3r+'
temp[ temp ==29 ] = 'Macrophage Marco+'
temp[ temp ==30 ] = 'Neutrophil Csf3r+'
temp[ temp ==31 ] = 'Osteoclast Ctsk+'
temp[ temp ==32 ] = 'Osteoclast Ctsk+'
temp[ temp ==33 ] = 'DC cDC2 Cd209a+'
obj$cellType = temp
Idents(obj) = obj$cellType
saveRDS( obj, 'sketch_v5_all/cam_it2.rds', compress = F )
obj = subset( obj, cellType != 'Cycling' )
saveRDS( obj, 'sketch_v5_all/cam_final.rds', compress = F )
makeFinalFigures( obj_cam, 'mouse_Myeloid',c('Clec9a','Cd209a', 'Ccr7', 'Tcf4', 'C1qc','Lyve1','Marco','Spp1','Cpa3','Ear2','Isg15','Ly6c2','Csf3r','Ctsk'))




############ CAF ############
obj = readRDS('sketch_v5_all/caf_final_joined.rds')
genes <- c("Pdgfra", "Lum", "Dcn",
           "Cd74","H2-Ab1", "Msln",
           "Isg15", "Cxcl10", "Ifit3",
           "Spp1", "Eno1", "Pgk1",
           "Lrrc15","Sdc1","Cthrc1",
           "Cxcl12", "Clec3b", "Igf1",
           "Pi16", "Dpp4", "Cd34")
neg <- c("Pecam1", "Rgs5", "Epcam",
           "Plp1","Ptprc")
makeSuppPanel( obj, 'caf', genes, neg, dotColor='black' )

############ CEC ############
obj = readRDS('sketch_v5_all/cec_final_joined.rds')
genes <- c("Pecam1", "Tie1", "Cdh5",
           "Gja5","Hey1", "Sema3g",
           "Car4", "Cd36","Pltp",
           "Prox1","Fxyd6","Pard6g",
           "Esm1","Inhbb","Nid2",
           "Ackr1", "Selp", "Nrp2",
           "Plvap", "Igfbp3", "Hspg2",
           "Col4a1", "Kdr", "Postn")
neg <- c("Pdgfra", "Rgs5","Epcam",
           "Plp1","Ptprc")
makeSuppPanel( obj, 'cec', genes, neg, dotColor = 'black' )

############ CMC ############
obj = readRDS('sketch_v5_all/cmc_final_joined.rds')
genes <- c("Rgs5", "Pdgfrb", "Kcnj8",
           "Cspg4","Procr", "Mcam",
           "Ccl2", "Ccl19", "Ccl11",
           "Cd248","Thy1","Col1a1",
           "Myh11", "Sorbs2", "Sncg")
neg <- c("Epcam", "Pdgfra", "Pecam1",
           "Plp1","Ptprc")
makeSuppPanel( obj, 'cmc', genes, neg, dotColor = 'black', rasterFeaturePlot = F )

############ CAN ############
obj = readRDS('sketch_v5_all/can_final_joined.rds')
genes <- c(#"PLP1", "GPM6B", "S100B",
           "Mpz","Mbp", "Ncmap",
           "Ngfr","Scn7a","Ncam1")
neg <- c("Epcam", "Pdgfra", "Pecam1",
           "Rgs5","Ptprc")
makeSuppPanel( obj, 'can', genes, neg, dotColor = 'black', rasterFeaturePlot = F )

############ CAL ############
obj = readRDS('sketch_v5_all/cal_final_joined.rds')
genes <- c("Ms4a1", "Cd79a","Ly6d",    # All B cell pops
           "Bank1","Ighm", "Cd19",  # Specific B cell pops
           "Klrd1","Fcer1g","Nkg7",    # NK
           "Jchain", "Mzb1", "Derl3",  # Plasma
           "Cd4", "Gata3", "Foxp3",    # Cd4 subpops
           "Cd8b1", "Isg15", "Pdcd1")  # CD8 subpops
neg <- c("Epcam", "Pdgfra", "Pecam1",
           "Rgs5","Plp1")
makeSuppPanel( obj, 'cal', genes, neg, dotColor = 'black' )

############ CAM ############
obj = readRDS('sketch_v5_all/cam_final.rds')
genes <- c("Clec9a", "Cd209a", "Ccr7", # DC subpops
           "Tcf4","Itgax", "Flt3",     # DC subpops (1), all (2-3)
           "C1qc", "Lyve1", "Marco",   # Mac subpop
           "Spp1", "Adgre1", "Csf1r",  # Mac subpop (10, all (2-3)
           "Ear2","Isg15","Ly6c2",     # Mono subpop
           "Cpa3", "Csf3r", "Ctsk")    # Mast, Neut, Osteo
neg <- c("Epcam", "Pdgfra", "Pecam1",
           "Rgs5","Plp1")
makeSuppPanel( obj, 'cam', genes, neg, dotColor = 'black' )



