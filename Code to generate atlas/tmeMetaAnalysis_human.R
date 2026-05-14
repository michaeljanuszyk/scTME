
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
  obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^MT-")
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

### format individual post-QC objects (streamline metadata, rename gene names); seurat objects are v5
for (i in 1:length(obj_list)) {; print(i)
  count <- count + 1
  obj <- obj_list[i]
  setTxtProgressBar(pb, count)
  x <- readRDS(paste0("processed/",obj,"_01_qc.rds"))
  if( dim(x)[2] <= 2 ) {; next; };
  x <- DietSeurat(x)
  x@meta.data <- x@meta.data %>%
    select(any_of(c("orig.ident", "nCount_RNA","nFeature_RNA","study",
                    "organ","cancer","model","sorting","site","percent.mt")))
  x <- RenameCells(x, new.names = paste(x$orig.ident, Cells(x), sep = "_"))
  # deal with duplicate gene names and non-standard gene names
  counts.mat <- GetAssayData( x, slot = 'counts' )
  geneNames <- rownames(counts.mat)
  synToName <- read.csv("GRCh38.p14_geneNames_geneSynonyms.csv", header = T)
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
synToName <- read.csv("GRCh38.p14_geneNames_geneSynonyms.csv", header = T)
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
  obj <- NormalizeData(obj)
  obj$study = first_elements[i]
  print(unique(obj$study))
  print(table(obj$orig.ident))
  caf_list[[i]] <- obj
  cell_count <- cell_count + dim(obj)[2]
}

print(cell_count)
merged <- merge(caf_list[[1]], caf_list[2:length(caf_list)])
saveRDS(merged, "sketch_v5_all/humanMerged.rds", compress = F)

#### Attach meta-data
obj = merged
clinical <- read.csv("human.meta.data.csv", header = T)
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
obj <- ScaleData(obj, verbose = F)
obj <- RunPCA(obj, npcs = 30, verbose = T)
saveRDS(obj, "sketch_v5_all/unfiltered_counts_obj_findvar_sketch1000_pca.rds", compress = F)

obj <- IntegrateLayers(obj, method = HarmonyIntegration, orig.reduction = "pca", new.reduction = "harmony", kmeans_init_nstart=20, kmeans_init_iter_max=5000, verbose = T)
saveRDS(obj, "sketch_v5_all/unfiltered_counts_obj_findvar_sketch1000_harmony.rds", compress = F)
obj <- RunUMAP(obj, reduction = "harmony", dims = 1:30, reduction.name = "umap.harmony", return.model = T, min.dist = 0.001)
saveRDS(obj, "sketch_v5_all/unfiltered_counts_obj_findvar_sketch1000_harmony_umap.rds",compress = F)
obj <- FindNeighbors(obj, reduction = "harmony", dims = 1:30)
for (i in c(0.06, 0.07, 0.08, 0.09, 0.1, 0.12, 0.14, 0.16, 0.018, 0.2)) {
  obj <- FindClusters(obj, resolution = i, graph.name = 'sketch_snn')
  DimPlot(obj, raster = T, label=T)
  ggsave(paste0("sketch_v5_all/unfiltered_counts_obj_findvar_sketch1000_harmony_umap_res",i,".jpg"), width = 5, height = 5, units = "in", limitsize = FALSE)
  saveRDS(obj, paste0("sketch_v5_all/unfiltered_counts_obj_findvar_sketch1000_harmony_umap_res",i,".rds"), compress = FALSE)
}

# Project to full dataset
obj$seurat_clusters = obj$sketch_snn_res.0.06
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
for (i in c(0.05,0.06,0.07,0.09,0.1,0.15,0.2) ) {
  obj <- FindClusters(obj, resolution = i, graph.name = 'RNA_snn')
  DimPlot(obj, raster = TRUE, label=T, reduction =  'umap.mnn')
  ggsave(paste0("sketch_v5_all/unfiltered_counts_obj_findvar_sketch1000_mnn_umap_res",i,".jpg"), width = 5, height = 5, units = "in", limitsize = FALSE)
  saveRDS(obj, paste0("sketch_v5_all/unfiltered_counts_obj_findvar_sketch1000_mnn_umap_res",i,".rds"), compress = FALSE)
}
obj <- readRDS("sketch_v5_all/unfiltered_counts_obj_findvar_sketch1000_mnn_umap_res0.5.rds")
# Project to full dataset
obj$seurat_clusters <- obj$RNA_snn_res.0.3; Idents(obj) <- obj$seurat_clusters 
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
for (i in c(0.05,0.06,0.07,0.09,0.1,0.15,0.2) ) {
  obj <- FindClusters(obj, resolution = i, graph.name = 'RNA_snn')
  DimPlot(obj, raster = T, reduction =  'umap.scvi', label=T, group.by = paste0('RNA_snn_res.',i))
  ggsave(paste0("sketch_v5_all/unfiltered_counts_obj_findvar_sketch1000_scvi_umap_res",i,".jpg"), width = 5, height = 5, units = "in", limitsize = FALSE)
  saveRDS(obj, paste0("sketch_v5_all/unfiltered_counts_obj_findvar_sketch1000_scvi_umap_res",i,".rds"), compress = FALSE)
}
obj$seurat_clusters <- obj$RNA_snn_res.0.2; Idents(obj) <- obj$seurat_clusters
# Project to full dataset
obj <- ProjectIntegration(object = obj, reduction = "integrated.scvi")
options(future.globals.maxSize = 8000 * 1024^2)
obj <- ProjectData(object = obj, sketched.reduction = "integrated.scvi.full", full.reduction = "integrated.scvi.full", umap.model = "umap.scvi", dims = 1:30, refdata = list(scvi.cluster.full = "seurat_clusters"))
saveRDS(obj, "sketch_v5_all/unfiltered_counts_obj_findvar_sketch1000_scvi_umap_project.rds",compress = F)


########################################################################
##################### Coarse Cell Type Annotations ##################### 
########################################################################

# All cells Harmony
temp = as.character( obj$harmony.cluster.full )
temp[ temp == 0 ] = 'Epithelial'
temp[ temp == 1 ] = 'Lymphoid'
temp[ temp == 2 ] = 'Myeloid'
temp[ temp == 3 ] = 'Fibroblast'
temp[ temp == 4 ] = 'Endothelial'
temp[ temp == 5 ] = 'Cycling'
temp[ temp == 6 ] = 'Mural'
temp[ temp == 7 ] = 'Lymphoid'
temp[ temp == 8 ] = 'Lymphoid'
temp[ temp == 9 ] = 'Nerve'
temp[ temp ==10 ] = 'Myeloid'
obj$cellType = temp; DefaultAssay(obj) = obj$cellType
saveRDS( obj, 'sketch_v5_all/allCells_harmony.rds', compress = F )

# All cells MNN
temp = as.character( obj$mnn.cluster.full )
temp[ temp == 0 ] = 'Epithelial'
temp[ temp == 1 ] = 'Lymphoid'
temp[ temp == 2 ] = 'Myeloid'
temp[ temp == 3 ] = 'Fibroblast'
temp[ temp == 4 ] = 'Epithelial'
temp[ temp == 5 ] = 'Endothelial'
temp[ temp == 6 ] = 'Epithelial'
temp[ temp == 7 ] = 'Epithelial'
temp[ temp == 8 ] = 'Mural'
temp[ temp == 9 ] = 'Nerve'
temp[ temp ==10 ] = 'Lymphoid'
temp[ temp ==11 ] = 'Epithelial'
temp[ temp ==12 ] = 'Lymphoid'
temp[ temp ==13 ] = 'Myeloid'
temp[ temp ==14 ] = 'Cycling'
temp[ temp ==15 ] = 'Myeloid'
obj$cellType = temp
Idents(obj) = obj$cellType
saveRDS( obj, 'sketch_v5_all/allCells_mnn.rds', compress = F )



# All cells scVI
temp = as.character( obj$scvi.cluster.full )
temp[ temp == 0 ] = 'Myeloid'
temp[ temp == 1 ] = 'Fibroblast'
temp[ temp == 2 ] = 'Lymphoid'
temp[ temp == 3 ] = 'Epithelial'
temp[ temp == 4 ] = 'Endothelial'
temp[ temp == 5 ] = 'Epithelial'
temp[ temp == 6 ] = 'Epithelial'
temp[ temp == 7 ] = 'Mural'
temp[ temp == 8 ] = 'Lymphoid'
temp[ temp == 9 ] = 'Epithelial'
temp[ temp ==10 ] = 'Lymphoid'
temp[ temp ==11 ] = 'Epithelial'
temp[ temp ==12 ] = 'Epithelial'
temp[ temp ==13 ] = 'Nerve'
temp[ temp ==14 ] = 'Lymphoid'
temp[ temp ==15 ] = 'Epithelial'
temp[ temp ==16 ] = 'Lymphoid'
temp[ temp ==17 ] = 'Myeloid'
temp[ temp ==18 ] = 'Epithelial'
temp[ temp ==19 ] = 'Epithelial'
temp[ temp ==20 ] = 'Fibroblast'
temp[ temp ==21 ] = 'Endothelial'
temp[ temp ==22 ] = 'Epithelial'
temp[ temp ==23 ] = 'Epithelial'
temp[ temp ==24 ] = 'Epithelial'
temp[ temp ==25 ] = 'Lymphoid'
temp[ temp ==26 ] = 'Epithelial'
temp[ temp ==27 ] = 'Epithelial'
temp[ temp ==28 ] = 'Epithelial'
temp[ temp ==29 ] = 'Epithelial'
temp[ temp ==30 ] = 'Lymphoid'
temp[ temp ==31 ] = 'Myeloid'
obj$cellType = temp
Idents(obj) = obj$cellType
saveRDS( obj, 'sketch_v5_all/allCells_scvi.rds', compress = F )



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



########################################################################
##################### Fine Cell Type Annotations ###################### 
########################################################################

# Fibroblast_it1
temp = as.character( obj$harmony.cluster.full )
temp[ temp == 0 ] = 'Epithelial'
temp[ temp == 1 ] = 'Fibroblast'
temp[ temp == 2 ] = 'Fibroblast'
temp[ temp == 3 ] = 'Fibroblast'
temp[ temp == 4 ] = 'Fibroblast'
temp[ temp == 5 ] = 'Lymphoid'
temp[ temp == 6 ] = 'Myeloid'
temp[ temp == 7 ] = 'Endothelial'
temp[ temp == 8 ] = 'Lymphoid'
temp[ temp == 9 ] = 'Nerve'
temp[ temp ==10 ] = 'Epithelial'
temp[ temp ==11 ] = 'Low Quality'
temp[ temp ==12 ] = 'Nerve'
temp[ temp ==13 ] = 'Low Quality'
temp[ temp ==14 ] = 'Low Quality'
obj$cellType = temp
Idents(obj) = obj$cellType
saveRDS( obj, 'sketch_v5_all/caf_it1.rds', compress = F )
obj = buildFromSamples(  target = 'Fibroblast', iteration = 2, useSketch = T, addIdents = F )
obj = iterativeFilter(     obj = obj, name = 'caf_it2', skipToScale = T )
obj = projectSketchedData( obj = obj, name = 'caf_it2', sketchedLabels = obj$cellType )

# Fibroblast_it2 
temp = as.character( obj$harmony.cluster.full )
temp[ temp == 0 ] = 'Fibroblast'
temp[ temp == 1 ] = 'Low Quality'
temp[ temp == 2 ] = 'Fibroblast'
temp[ temp == 3 ] = 'Mural'
temp[ temp == 4 ] = 'Fibroblast'
temp[ temp == 5 ] = 'Low Quality' # Doublet?
temp[ temp == 6 ] = 'Fibroblast'
temp[ temp == 7 ] = 'Fibroblast'
temp[ temp == 8 ] = 'Cycling'
temp[ temp == 9 ] = 'Low Quality'
temp[ temp ==10 ] = 'Low Quality'
temp[ temp ==11 ] = 'Low Quality'
temp[ temp ==12 ] = 'Low Quality'
obj$cellType = temp
Idents(obj) = obj$cellType
saveRDS( obj, 'sketch_v5_all/caf_it2.rds', compress = F )
obj = buildFromSamples(  target = 'Fibroblast', iteration = 3, useSketch = T, addIdents = F )
obj = iterativeFilter(     obj = obj, name = 'caf_it3', skipToScale = T )
obj = projectSketchedData( obj = obj, name = 'caf_it3', sketchedLabels = obj$cellType )

# Fibroblast_it3
temp = as.character( obj$harmony.cluster.full )
temp[ temp == 0 ] = 'Fibroblast'
temp[ temp == 1 ] = 'Fibroblast'
temp[ temp == 2 ] = 'Low Quality'
temp[ temp == 3 ] = 'Fibroblast'
temp[ temp == 4 ] = 'Fibroblast'
temp[ temp == 5 ] = 'Fibroblast'
temp[ temp == 6 ] = 'Fibroblast'
temp[ temp == 7 ] = 'Fibroblast'
temp[ temp == 8 ] = 'Fibroblast'
temp[ temp == 9 ] = 'Fibroblast'
temp[ temp ==10 ] = 'Fibroblast'
temp[ temp ==11 ] = 'Fibroblast'
temp[ temp ==12 ] = 'Fibroblast'
temp[ temp ==13 ] = 'Low Quality'
temp[ temp ==14 ] = 'Fibroblast'
temp[ temp ==15 ] = 'Low Quality'
temp[ temp ==16 ] = 'Low Quality'
obj$cellType = temp
Idents(obj) = obj$cellType
saveRDS( obj, 'sketch_v5_all/caf_it3.rds', compress = F )
obj = subset( obj, cellType == 'Fibroblast' )
obj = analyzeObj( obj, 'caf_it4' )

# Fibroblast_it4
temp = as.character( obj$harmony.cluster.full )
temp[ temp == 0 ] = 'Fibroblast' # 'mCAF LRRC15+' # CTHRC1
temp[ temp == 1 ] = 'Fibroblast' # 'ssCAF PI16+'
temp[ temp == 2 ] = 'Low Quality' # LQ
temp[ temp == 3 ] = 'Low Quality' # Lyedig cells?
temp[ temp == 4 ] = 'Fibroblast'
temp[ temp == 5 ] = 'Fibroblast'
temp[ temp == 6 ] = 'Fibroblast'
temp[ temp == 7 ] = 'Fibroblast'
temp[ temp == 8 ] = 'Fibroblast'
temp[ temp == 9 ] = 'Fibroblast'
temp[ temp ==10 ] = 'Fibroblast'
temp[ temp ==11 ] = 'Low Quality' # CDH19+ LAMA2+ Fibroblast?
temp[ temp ==12 ] = 'Low Quality'
temp[ temp ==13 ] = 'Low Quality' # Pericyte
temp[ temp ==14 ] = 'Low Quality' # Doublet?
obj$cellType = temp
Idents(obj) = obj$cellType
saveRDS( obj, 'sketch_v5_all/caf_it4.rds', compress = F )
obj = subset( obj, cellType == 'Fibroblast' )
obj = analyzeObj( obj, 'caf_it4' )

# Fibroblast_it5
temp = as.character( obj$harmony.cluster.full )
temp[ temp == 0 ] = 'Fibroblast'
temp[ temp == 1 ] = 'Fibroblast'
temp[ temp == 2 ] = 'Fibroblast'
temp[ temp == 3 ] = 'Fibroblast'
temp[ temp == 4 ] = 'Fibroblast'
temp[ temp == 5 ] = 'Tumor'
temp[ temp == 6 ] = 'Fibroblast'
temp[ temp == 7 ] = 'Fibroblast'
temp[ temp == 8 ] = 'Fibroblast'
temp[ temp == 9 ] = 'Fibroblast'
temp[ temp ==10 ] = 'Fibroblast'
temp[ temp ==11 ] = 'Low Quality'
obj$cellType = temp
Idents(obj) = obj$cellType
saveRDS( obj, 'sketch_v5_all/caf_it5.rds', compress = F )
obj = subset( obj, cellType == 'Fibroblast' )
obj = analyzeObj( obj, 'caf_it5' )

# Fibroblast_it6
temp = as.character( obj$harmony.cluster.full )
temp[ temp == 0 ] = 'mCAF LRRC15+'
temp[ temp == 1 ] = 'ssCAF PI16+'
temp[ temp == 2 ] = 'iCAF CXCL8+'
temp[ temp == 3 ] = 'mCAF LRRC15+'
temp[ temp == 4 ] = 'mCAF COL4A1+'
temp[ temp == 5 ] = 'ssCAF PI16+'
temp[ temp == 6 ] = 'ssCAF PI16+'
temp[ temp == 7 ] = 'Low Quality'
temp[ temp == 8 ] = 'iCAF CXCL8+' 
temp[ temp == 9 ] = 'mCAF COL4A1+'
temp[ temp ==10 ] = 'ssCAF CXCL14+'
temp[ temp ==11 ] = 'iCAF ISG15+'
temp[ temp ==12 ] = 'apCAF CD74+'
temp[ temp ==13 ] = 'ssCAF CXCL14+'
temp[ temp ==14 ] = 'ssCAF PI16+'
temp[ temp ==15 ] = 'mCAF LRRC15+'
temp[ temp ==16 ] = 'ssCAF CXCL14+'
temp[ temp ==17 ] = 'mCAF LRRC15+'
temp[ temp ==18 ] = 'Low Quality'
temp[ temp ==19 ] = 'mCAF COL4A1+'
temp[ temp ==20 ] = 'mCAF LRRC15+'
temp[ temp ==21 ] = 'mCAF COL4A1+'
temp[ temp ==22 ] = 'ssCAF CXCL14+'
obj$cellType = temp
Idents(obj) = obj$cellType
saveRDS( obj, 'sketch_v5_all/caf_it6.rds', compress = F )
obj = subset( obj, cellType != 'Low Quality')
saveRDS( obj, 'sketch_v5_all/caf_final.rds', compress = F )



# Endothelial_it1
temp = as.character( obj$harmony.cluster.full )
temp[ temp == 0 ] = 'Endothelial'
temp[ temp == 1 ] = 'Epithelial'
temp[ temp == 2 ] = 'Lymphoid'
temp[ temp == 3 ] = 'Myeloid'
temp[ temp == 4 ] = 'Fibroblast'
temp[ temp == 5 ] = 'Mural'
temp[ temp == 6 ] = 'Cycling'
temp[ temp == 7 ] = 'Endothelial'
temp[ temp == 8 ] = 'Lymphoid'
temp[ temp == 9 ] = 'Lymphoid'
temp[ temp ==10 ] = 'Myeloid'
obj$cellType = temp
Idents(obj) = obj$cellType
saveRDS( obj, 'sketch_v5_all/cec_it1.rds', compress = F )
obj = subset( obj, cellType == 'Endothelial' )
obj = analyzeObj( obj, 'cec_it2' )

# Endothelial_it2
temp = as.character( obj$RNA_snn_res.0.2 )
temp[ temp == 0 ] = 'Endothelial'
temp[ temp == 1 ] = 'Low Quality'
temp[ temp == 2 ] = 'Endothelial'
temp[ temp == 3 ] = 'Endothelial'
temp[ temp == 4 ] = 'Endothelial'
temp[ temp == 5 ] = 'Endothelial'
temp[ temp == 6 ] = 'Epithelial'
temp[ temp == 7 ] = 'Fibroblast'
temp[ temp == 8 ] = 'Endothelial'
temp[ temp == 9 ] = 'Endothelial'
temp[ temp ==10 ] = 'Endothelial'
temp[ temp ==11 ] = 'Low Quality'
temp[ temp ==12 ] = 'Low Quality'
temp[ temp ==13 ] = 'Immune'
temp[ temp ==14 ] = 'Epithelial'
obj$cellType = temp
saveRDS( obj, 'sketch_v5_all/cec_it2.rds', compress = F )
obj = subset( obj, cellType == 'Endothelial' )
obj = analyzeObj( obj, 'cec_it3' )

# Endothelial_it3
temp = as.character( obj$RNA_snn_res.0.16 )
temp[ temp == 0 ] = 'Endothelial'
temp[ temp == 1 ] = 'Endothelial'
temp[ temp == 2 ] = 'Endothelial'
temp[ temp == 3 ] = 'Endothelial'
temp[ temp == 4 ] = 'Endothelial'
temp[ temp == 5 ] = 'Low Quality'
temp[ temp == 6 ] = 'Low Quality'
temp[ temp == 7 ] = 'Low Quality'
temp[ temp == 8 ] = 'Low Quality' 
temp[ temp == 9 ] = 'Low Quality' 
temp[ temp ==10 ] = 'Low Quality' 
temp[ temp ==11 ] = 'Low Quality'
temp[ temp ==12 ] = 'Low Quality'
obj$cellType = temp
saveRDS( obj, 'sketch_v5_all/cec_it3.rds', compress = F )
obj = subset( obj, cellType == 'Endothelial' )
obj = analyzeObj( obj, 'cec_it4' )

# Endothelial_it4
temp = as.character( obj$RNA_snn_res.0.09 )
temp[ temp == 0 ] = 'Vein ACKR1+'
temp[ temp == 1 ] = 'Tip ESM1+'
temp[ temp == 2 ] = 'Artery GJA5+'
temp[ temp == 3 ] = 'Lymphatic PROX1+'
temp[ temp == 4 ] = 'Capillary CA4+'
obj$cellType = temp
Idents(obj) = obj$cellType
Idents(obj) <- factor(x = Idents(obj), levels = sort(levels(obj), decreasing=F))
saveRDS( obj, 'sketch_v5_all/cec_it4.rds', compress = F )
saveRDS( obj, 'sketch_v5_all/cec_final.rds', compress = F )



# Mural_it1
temp = as.character( obj$harmony.cluster.full )
temp[ temp == 0 ] = 'Low Quality'
temp[ temp == 1 ] = 'Mural'
temp[ temp == 2 ] = 'Mural'
temp[ temp == 3 ] = 'Fibroblast'
temp[ temp == 4 ] = 'Mural'
temp[ temp == 5 ] = 'Lymphoid'
temp[ temp == 6 ] = 'Epithelial'
temp[ temp == 7 ] = 'Endothelial'
temp[ temp == 8 ] = 'Myeloid'
temp[ temp == 9 ] = 'Epithelial'
temp[ temp ==10 ] = 'Mural'
temp[ temp ==11 ] = 'Mural'
temp[ temp ==12 ] = 'Fibroblast'
temp[ temp ==13 ] = 'Mural'
temp[ temp ==14 ] = 'Mural'
temp[ temp ==15 ] = 'Lymphoid'
temp[ temp ==16 ] = 'Endothelial'
temp[ temp ==17 ] = 'Mural'
temp[ temp ==18 ] = 'Low Quality'
temp[ temp ==19 ] = 'Low Quality'
temp[ temp ==20 ] = 'Low Quality'
temp[ temp ==21 ] = 'Low Quality'
temp[ temp ==22 ] = 'Low Quality'
temp[ temp ==23 ] = 'Low Quality'
temp[ temp ==24 ] = 'Low Quality'
temp[ temp ==25 ] = 'Low Quality'
temp[ temp ==26 ] = 'Low Quality'
temp[ temp ==27 ] = 'Low Quality'
temp[ temp ==28 ] = 'Low Quality'
temp[ temp ==29 ] = 'Low Quality'
temp[ temp ==30 ] = 'Low Quality'
obj$cellType = temp
Idents(obj) = obj$cellType
saveRDS( obj, 'sketch_v5_all/cmc_it1.rds', compress = F )
obj = subset( obj, cellType == 'Mural' )
obj = analyzeObj( obj, 'cmc_it2' )


# Mural_it2
temp = as.character( obj$harmony.cluster.full )
temp[ temp == 0 ] = 'Mural'
temp[ temp == 1 ] = 'Mural'
temp[ temp == 2 ] = 'Mural' 
temp[ temp == 3 ] = 'Mural'
temp[ temp == 4 ] = 'Mural'
temp[ temp == 5 ] = 'Cycling'
temp[ temp == 6 ] = 'Lymphoid'
temp[ temp == 7 ] = 'Nerve'
temp[ temp == 8 ] = 'Cycling'
temp[ temp == 9 ] = 'Low Quality'
temp[ temp ==10 ] = 'Low Quality'
obj$cellType = temp
Idents(obj) = obj$cellType
saveRDS( obj, 'sketch_v5_all/cmc_it2.rds', compress = F )
obj = subset( obj, cellType == 'Mural' )
obj = analyzeObj( obj, 'cmc_it3' )

# Mural_it3
temp = as.character( obj$harmony.cluster.full )
temp[ temp == 0 ] = 'Mural'
temp[ temp == 1 ] = 'Mural'
temp[ temp == 2 ] = 'Mural'
temp[ temp == 3 ] = 'Mural' 
temp[ temp == 4 ] = 'Mural'
temp[ temp == 5 ] = 'Mural'
temp[ temp == 6 ] = 'Mural'
temp[ temp == 7 ] = 'Mural'
temp[ temp == 8 ] = 'Mural'
temp[ temp == 9 ] = 'Mural'
temp[ temp ==10 ] = 'Low Quality'
temp[ temp ==11 ] = 'Low Quality'
temp[ temp ==12 ] = 'Mural'
temp[ temp ==13 ] = 'Mural'
temp[ temp ==14 ] = 'Mural'
temp[ temp ==15 ] = 'Low Quality'
temp[ temp ==16 ] = 'Mural'
temp[ temp ==17 ] = 'Mural'
temp[ temp ==18 ] = 'Low Quality'
temp[ temp ==19 ] = 'Low Quality'
temp[ temp ==20 ] = 'Low Quality'
temp[ temp ==21 ] = 'Low Quality'
temp[ temp ==22 ] = 'Low Quality'
obj$cellType = temp
Idents(obj) = obj$cellType
saveRDS( obj, 'sketch_v5_all/cmc_it3.rds', compress = F )
obj = subset( obj, cellType == 'Mural' )
obj = analyzeObj( obj, 'cmc_it4' )

# Mural_it4
temp = as.character( obj$harmony.cluster.full )
temp[ temp == 0 ] = 'Mural'
temp[ temp == 1 ] = 'Mural'
temp[ temp == 2 ] = 'Low Quality'
temp[ temp == 3 ] = 'Mural'
temp[ temp == 4 ] = 'Mural'
temp[ temp == 5 ] = 'Mural'
temp[ temp == 6 ] = 'Mural'
temp[ temp == 7 ] = 'Mural'
temp[ temp == 8 ] = 'Mural'
temp[ temp == 9 ] = 'Mural'
temp[ temp ==10 ] = 'Mural'
temp[ temp ==11 ] = 'Mural'
temp[ temp ==12 ] = 'Mural'
temp[ temp ==13 ] = 'Mural'
temp[ temp ==14 ] = 'Mural'
temp[ temp ==15 ] = 'Mural'
temp[ temp ==16 ] = 'Low Quality'
obj$cellType = temp
Idents(obj) = obj$cellType
saveRDS( obj, 'sketch_v5_all/cmc_it4.rds', compress = F )
obj = subset( obj, cellType == 'Mural' )
obj = analyzeObj( obj, 'cmc_it5' )

# Mural_it5
temp = as.character( obj$harmony.cluster.full )
temp[ temp == 0 ] = 'Pericyte CD248+'
temp[ temp == 1 ] = 'SMC vascular RERGL+'
temp[ temp == 2 ] = 'Pericyte CD248+'
temp[ temp == 3 ] = 'SMC vascular RERGL+'
temp[ temp == 4 ] = 'SMC vascular RERGL+'
temp[ temp == 5 ] = 'SMC vascular RERGL+'
temp[ temp == 6 ] = 'Pericyte CD248+'
temp[ temp == 7 ] = 'Pericyte CCL2+'
temp[ temp == 8 ] = 'Pericyte CD248+'
temp[ temp == 9 ] = 'SMC vascular WFDC1+'
temp[ temp ==10 ] = 'Pericyte CD248+'
temp[ temp ==11 ] = 'Pericyte ISG15+'
temp[ temp ==12 ] = 'Pericyte CD248+'
temp[ temp ==13 ] = 'Pericyte CCL2+'
temp[ temp ==14 ] = 'Pericyte CCL2+'
temp[ temp ==15 ] = 'Pericyte CD248+'
temp[ temp ==16 ] = 'SMC vascular WFDC1+'
temp[ temp ==17 ] = 'Pericyte CD248+'
temp[ temp ==18 ] = 'Pericyte CCL2+'
temp[ temp ==19 ] = 'SMC vascular WFDC1+'
temp[ temp ==20 ] = 'Low Quality' 
temp[ temp ==21 ] = 'SMC visceral DES+'
temp[ temp ==22 ] = 'SMC vascular WFDC1+'
temp[ temp ==23 ] = 'SMC visceral DES+'
temp[ temp ==24 ] = 'Pericyte CCL2+'
temp[ temp ==25 ] = 'Pericyte CCL2+'
temp[ temp ==26 ] = 'Pericyte CCL2+'
temp[ temp ==27 ] = 'Pericyte CD248+'
temp[ temp ==28 ] = 'Pericyte CCL2+'
temp[ temp ==29 ] = 'SMC vascular WFDC1+'
temp[ temp ==30 ] = 'SMC vascular WFDC1+'
temp[ temp ==31 ] = 'Low Quality'
temp[ temp ==32 ] = 'Low Quality'
temp[ temp ==33 ] = 'Low Quality'
obj$cellType = temp
Idents(obj) = obj$cellType
saveRDS( obj, 'sketch_v5_all/cmc_it5.rds', compress = F )
obj = subset( obj, cellType != 'Low Quality' )
saveRDS( obj, 'sketch_v5_all/cmc_final.rds', compress = F )



# Nerve_it1
temp = as.character( obj$harmony.cluster.full )
temp[ temp == 0 ] = 'Epithelial'
temp[ temp == 1 ] = 'Epithelial'
temp[ temp == 2 ] = 'Nerve'
temp[ temp == 3 ] = 'Lymphoid'
temp[ temp == 4 ] = 'Fibroblast'
temp[ temp == 5 ] = 'Myeloid'
temp[ temp == 6 ] = 'Fibroblast'
temp[ temp == 7 ] = 'Endothelial'
temp[ temp == 8 ] = 'Epithelial'
temp[ temp == 9 ] = 'Cycling'
temp[ temp ==10 ] = 'Nerve'
temp[ temp ==11 ] = 'Mural'
temp[ temp ==12 ] = 'Nerve'
temp[ temp ==13 ] = 'Fibroblast'
temp[ temp ==14 ] = 'Low Quality'
temp[ temp ==15 ] = 'Mural'
temp[ temp ==16 ] = 'Mural' 
temp[ temp ==17 ] = 'Epithelial'
temp[ temp ==18 ] = 'Epithelial'
temp[ temp ==19 ] = 'Low Quality' 
temp[ temp ==20 ] = 'Low Quality'
obj$cellType = temp
Idents(obj) = obj$cellType
saveRDS( obj, 'sketch_v5_all/can_it1.rds', compress = F )
obj = subset( obj, cellType == 'Nerve' )
obj = analyzeObj( obj, 'can_it2' )

# Nerve_it2
temp = as.character( obj$RNA_snn_res.0.2 )
temp[ temp == 0 ] = 'Nerve'
temp[ temp == 1 ] = 'Nerve'
temp[ temp == 2 ] = 'Low Quality' 
temp[ temp == 3 ] = 'Nerve' 
temp[ temp == 4 ] = 'Epithelial'
temp[ temp == 5 ] = 'Nerve' 
temp[ temp == 6 ] = 'Lymphoid'
temp[ temp == 7 ] = 'Nerve'
temp[ temp == 8 ] = 'Nerve'
temp[ temp == 9 ] = 'Myeloid'
temp[ temp ==10 ] = 'Nerve'
temp[ temp ==11 ] = 'Nerve'
temp[ temp ==12 ] = 'Myeloid'
temp[ temp ==13 ] = 'Myeloid'
temp[ temp ==14 ] = 'Low Quality' 
temp[ temp ==15 ] = 'Nerve'
temp[ temp ==16 ] = 'Low Quality' 
temp[ temp ==17 ] = 'Fibroblast'
temp[ temp ==18 ] = 'Low Quality'
temp[ temp ==19 ] = 'Low Quality' 
temp[ temp ==20 ] = 'Low Quality'
temp[ temp ==21 ] = 'Low Quality'
temp[ temp ==22 ] = 'Low Quality'
temp[ temp ==23 ] = 'Low Quality'
temp[ temp ==24 ] = 'Low Quality'
temp[ temp ==25 ] = 'Low Quality'
temp[ temp ==26 ] = 'Low Quality' 
temp[ temp ==27 ] = 'Low Quality'
temp[ temp ==28 ] = 'Low Quality'
temp[ temp ==29 ] = 'Low Quality' 
temp[ temp ==30 ] = 'Low Quality'
obj$cellType = temp
Idents(obj) = obj$cellType
saveRDS( obj, 'sketch_v5_all/can_it2.rds', compress = F )
obj = subset( obj, cellType == 'Nerve' )
obj = analyzeObj( obj, 'can_it3' )

# Nerve_it3
temp = as.character( obj$RNA_snn_res.0.1 )
temp[ temp == 0 ] = 'Low Quality'
temp[ temp == 1 ] = 'Low Quality'
temp[ temp == 2 ] = 'Low Quality'
temp[ temp == 3 ] = 'Nerve'
temp[ temp == 4 ] = 'Low Quality'
temp[ temp == 5 ] = 'Nerve'
temp[ temp == 6 ] = 'Lymphoid'
temp[ temp == 7 ] = 'Nerve'
temp[ temp == 8 ] = 'Myeloid'
temp[ temp == 9 ] = 'Low Quality'
temp[ temp ==10 ] = 'Nerve'
temp[ temp ==11 ] = 'Low Quality'
temp[ temp ==12 ] = 'Low Quality'
temp[ temp ==13 ] = 'Fibroblast'
temp[ temp ==14 ] = 'Fibroblast'
temp[ temp ==15 ] = 'Low Quality' 
temp[ temp ==16 ] = 'Low Quality'
obj$cellType = temp
Idents(obj) = obj$cellType
saveRDS( obj, 'sketch_v5_all/can_it3.rds', compress = F )
obj = subset( obj, cellType == 'Nerve' )
obj = analyzeObj( obj, 'can_it4' )

# Nerve_it4
temp = as.character( obj$RNA_snn_res.0.4 )
temp[ temp == 0 ] = 'Nerve' 
temp[ temp == 1 ] = 'Immune'
temp[ temp == 2 ] = 'Nerve'
temp[ temp == 3 ] = 'Nerve'
temp[ temp == 4 ] = 'Low Quality'
temp[ temp == 5 ] = 'Low Quality'
temp[ temp == 6 ] = 'Low Quality'
temp[ temp == 7 ] = 'Low Quality'
temp[ temp == 8 ] = 'Nerve'
temp[ temp == 9 ] = 'Low Quality'
temp[ temp ==10 ] = 'Low Quality'
temp[ temp ==11 ] = 'Low Quality'
temp[ temp ==12 ] = 'Nerve'
temp[ temp ==13 ] = 'Low Quality'
temp[ temp ==14 ] = 'Low Quality'
temp[ temp ==15 ] = 'Low Quality'
temp[ temp ==16 ] = 'Low Quality'
temp[ temp ==17 ] = 'Epithelial'
temp[ temp ==18 ] = 'Immune'
temp[ temp ==19 ] = 'Low Quality'
temp[ temp ==20 ] = 'Low Quality'
temp[ temp ==21 ] = 'Low Quality'
temp[ temp ==22 ] = 'Low Quality'
obj$cellType = temp
Idents(obj) = obj$cellType
saveRDS( obj, 'sketch_v5_all/can_it4.rds', compress = F )
obj = subset( obj, cellType == 'Nerve' )
obj = analyzeObj( obj, 'can_it5' )

# Nerve_it5
temp = as.character( obj$RNA_snn_res.0.4 )
temp[ temp == 0 ] = 'Nerve' 
temp[ temp == 1 ] = 'Low Quality'
temp[ temp == 2 ] = 'Nerve' 
temp[ temp == 3 ] = 'Nerve'
temp[ temp == 4 ] = 'Lymphoid' 
temp[ temp == 5 ] = 'Low Quality' 
temp[ temp == 6 ] = 'Nerve'
temp[ temp == 7 ] = 'Low Quality'
temp[ temp == 8 ] = 'Lymphoid'
temp[ temp == 9 ] = 'Nerve'
temp[ temp ==10 ] = 'Low Quality'
temp[ temp ==11 ] = 'Nerve'
temp[ temp ==12 ] = 'Mural'
temp[ temp ==13 ] = 'Low Quality'
obj$cellType = temp
Idents(obj) = obj$cellType
saveRDS( obj, 'sketch_v5_all/can_it5.rds', compress = F )
obj = subset( obj, cellType == 'Nerve' )
obj = analyzeObj( obj, 'can_it6' )

# Nerve_it6
temp = as.character( obj$RNA_snn_res.0.06 )
temp[ temp == 0 ] = 'Nerve' # MLANA 
temp[ temp == 1 ] = 'Nerve' # Schwann
temp[ temp == 2 ] = 'Low Quality'
obj$cellType = temp
Idents(obj) = obj$cellType
saveRDS( obj, 'sketch_v5_all/can_it6.rds', compress = F )
obj = subset( obj, cellType == 'Nerve' )
obj = analyzeObj( obj, 'can_it7' )

# Nerve_it7
temp = as.character( obj$RNA_snn_res.0.2 )
temp[ temp == 0 ] = 'Melanocyte MLANA+'
temp[ temp == 1 ] = 'Melanocyte MLANA+'
temp[ temp == 2 ] = 'Schwann myel MPZ+'
temp[ temp == 3 ] = 'Schwann non-myel NGFR+'
temp[ temp == 4 ] = 'Schwann myel MPZ+'
obj$cellType = temp
Idents(obj) = obj$cellType
saveRDS( obj, 'sketch_v5_all/can_it7.rds', compress = F )
saveRDS( obj, 'sketch_v5_all/can_final.rds', compress = F )



# Lymphoid_it1
temp = as.character( obj$harmony.cluster.full )
temp[ temp == 0 ] = 'Lymphoid'
temp[ temp == 1 ] = 'Epithelial'
temp[ temp == 2 ] = 'Myeloid'
temp[ temp == 3 ] = 'Fibroblast'
temp[ temp == 4 ] = 'Lymphoid'
temp[ temp == 5 ] = 'Lymphoid'
temp[ temp == 6 ] = 'Endothelial'
temp[ temp == 7 ] = 'Cycling'
temp[ temp == 8 ] = 'Myeloid'
temp[ temp == 9 ] = 'Low Quality'
temp[ temp ==10 ] = 'Nerve'
temp[ temp ==11 ] = 'Endothelial'
temp[ temp ==12 ] = 'Epithelial'
temp[ temp ==13 ] = 'Low Quality'
obj$cellType = temp
Idents(obj) = obj$cellType
saveRDS( obj, 'sketch_v5_all/cal_it1.rds', compress = F )
obj = subset( obj, cellType == 'Lymphoid' )
obj = analyzeObj( obj, 'cal_it2' )

# Lymphoid_it2
temp = as.character( obj$harmony.cluster.full )
temp[ temp == 0 ] = 'Lymphoid'
temp[ temp == 1 ] = 'Epithelial'
temp[ temp == 2 ] = 'Lymphoid'
temp[ temp == 3 ] = 'Lymphoid'
temp[ temp == 4 ] = 'Lymphoid'
temp[ temp == 5 ] = 'Myeloid'
temp[ temp == 6 ] = 'Low Quality' # Doublet?
temp[ temp == 7 ] = 'Low Quality'
obj$cellType = temp
Idents(obj) = obj$cellType
saveRDS( obj, 'sketch_v5_all/cal_it2.rds', compress = F )
obj = subset( obj, cellType == 'Lymphoid' )
obj = analyzeObj( obj, 'cal_it3' )

# Lymphoid_it3
temp = as.character( obj$harmony.cluster.full )
temp[ temp == 0 ] = 'Low Quality'
temp[ temp == 1 ] = 'Low Quality'
temp[ temp == 2 ] = 'Lymphoid'
temp[ temp == 3 ] = 'Lymphoid'
temp[ temp == 4 ] = 'Low Quality'
temp[ temp == 5 ] = 'Lymphoid'
temp[ temp == 6 ] = 'Lymphoid'
temp[ temp == 7 ] = 'Lymphoid'
temp[ temp == 8 ] = 'Lymphoid'
temp[ temp == 9 ] = 'Epithelial'
temp[ temp ==10 ] = 'Lymphoid'
temp[ temp ==11 ] = 'Lymphoid'
temp[ temp ==12 ] = 'Low Quality'
temp[ temp ==13 ] = 'Lymphoid'
temp[ temp ==14 ] = 'Lymphoid'
temp[ temp ==15 ] = 'Lymphoid'
temp[ temp ==16 ] = 'Lymphoid'
temp[ temp ==17 ] = 'Lymphoid'
temp[ temp ==18 ] = 'Lymphoid'
temp[ temp ==19 ] = 'Lymphoid'
temp[ temp ==20 ] = 'Low Quality'
temp[ temp ==21 ] = 'Lymphoid'
temp[ temp ==22 ] = 'Low Quality'
temp[ temp ==23 ] = 'Lymphoid'
temp[ temp ==24 ] = 'Lymphoid'
temp[ temp ==25 ] = 'Low Quality'
temp[ temp ==26 ] = 'Low Quality'
temp[ temp ==27 ] = 'Lymphoid'
temp[ temp ==28 ] = 'Low Quality'
temp[ temp ==29 ] = 'Low Quality'
obj$cellType = temp
Idents(obj) = obj$cellType
saveRDS( obj, 'sketch_v5_all/cal_it3.rds', compress = F )
obj = subset( obj, cellType == 'Lymphoid' )
obj = analyzeObj( obj, 'cal_it4' )

# Lymphoid_it4
temp = as.character( obj$harmony.cluster.full )
temp[ temp == 0 ] = 'Lymphoid'
temp[ temp == 1 ] = 'Lymphoid'
temp[ temp == 2 ] = 'Lymphoid'
temp[ temp == 3 ] = 'Low Quality'
temp[ temp == 4 ] = 'Lymphoid'
temp[ temp == 5 ] = 'Lymphoid'
temp[ temp == 6 ] = 'Lymphoid'
temp[ temp == 7 ] = 'Lymphoid'
temp[ temp == 8 ] = 'Lymphoid'
temp[ temp == 9 ] = 'Lymphoid'
temp[ temp ==10 ] = 'Lymphoid'
temp[ temp ==11 ] = 'Lymphoid'
temp[ temp ==12 ] = 'Low Quality'
temp[ temp ==13 ] = 'Lymphoid'
temp[ temp ==14 ] = 'Lymphoid'
temp[ temp ==15 ] = 'Lymphoid'
temp[ temp ==16 ] = 'Cycling'
temp[ temp ==17 ] = 'Lymphoid'
temp[ temp ==18 ] = 'Lymphoid'
temp[ temp ==19 ] = 'Low Quality'
temp[ temp ==20 ] = 'Low Quality'
temp[ temp ==21 ] = 'Low Quality'
obj$cellType = temp
Idents(obj) = obj$cellType
saveRDS( obj, 'sketch_v5_all/cal_it4.rds', compress = F )
obj = subset( obj, cellType == 'Lymphoid' )
obj = analyzeObj( obj, 'cal_it5' )

# Lymphoid_it5
temp = as.character( obj$harmony.cluster.full )
temp[ temp == 0 ] = 'T CD8+/GZMK+'
temp[ temp == 1 ] = 'Treg FOXP3+'
temp[ temp == 2 ] = 'T CD8+/GZMK+'
temp[ temp == 3 ] = 'NK KLRD1+'
temp[ temp == 4 ] = 'Plasma JCHAIN+'
temp[ temp == 5 ] = 'T CD4+/IL7R+'
temp[ temp == 6 ] = 'B mem BANK1'
temp[ temp == 7 ] = 'B naive IGHM+'
temp[ temp == 8 ] = 'NK KLRD1+'
temp[ temp == 9 ] = 'T CD8+/GZMB+'
temp[ temp ==10 ] = 'T CD4+/IL7R+'
temp[ temp ==11 ] = 'Plasma JCHAIN+'
temp[ temp ==12 ] = 'T CD8+/GZMB+'
temp[ temp ==13 ] = 'Plasma JCHAIN+'
temp[ temp ==14 ] = 'T CD8+/GZMK+'
temp[ temp ==15 ] = 'Tfh CXCL13+'
temp[ temp ==16 ] = 'Plasma JCHAIN+'
temp[ temp ==17 ] = 'T CD8+/ISG15+'
temp[ temp ==18 ] = 'B germinal RGS13+'
temp[ temp ==19 ] = 'Low Quality'
temp[ temp ==20 ] = 'Low Quality'
temp[ temp ==21 ] = 'Low Quality'
temp[ temp ==22 ] = 'Low Quality'
temp[ temp ==23 ] = 'Low Quality'
obj$cellType = temp
Idents(obj) = obj$cellType
saveRDS( obj, 'sketch_v5_all/cal_it5.rds', compress = F )
obj = subset( obj, cellType != 'Low Quality' )
saveRDS( obj, 'sketch_v5_all/cal_final.rds', compress = F )



# Myeloid_it1
temp = as.character( obj$harmony.cluster.full )
temp[ temp == 0 ] = 'Myeloid'
temp[ temp == 1 ] = 'Low Quality' # LQ
temp[ temp == 2 ] = 'Epithelial'
temp[ temp == 3 ] = 'Lymphoid'
temp[ temp == 4 ] = 'Myeloid'
temp[ temp == 5 ] = 'Fibroblast'
temp[ temp == 6 ] = 'Myeloid'
temp[ temp == 7 ] = 'Cycling'
temp[ temp == 8 ] = 'Myeloid'
temp[ temp == 9 ] = 'Endothelial'
temp[ temp ==10 ] = 'Lymphoid'
temp[ temp ==11 ] = 'Lympohid'
temp[ temp ==12 ] = 'Low Quality' # Doublet?
obj$cellType = temp
Idents(obj) = obj$cellType
saveRDS( obj, 'sketch_v5_all/cam_it1.rds', compress = F )
obj = subset( obj, cellType == 'Myeloid' )
obj = analyzeObj( obj, 'cam_it2' )

# Myeloid_it2
temp = as.character( obj$harmony.cluster.full )
temp[ temp == 0 ] = 'Myeloid'
temp[ temp == 1 ] = 'Myeloid'
temp[ temp == 2 ] = 'Myeloid'
temp[ temp == 3 ] = 'Myeloid'
temp[ temp == 4 ] = 'Myeloid'
temp[ temp == 5 ] = 'Low Quality'
temp[ temp == 6 ] = 'Myeloid'
temp[ temp == 7 ] = 'Myeloid'
temp[ temp == 8 ] = 'Lymphoid'
temp[ temp == 9 ] = 'Myeloid'
temp[ temp ==10 ] = 'Myeloid'
temp[ temp ==11 ] = 'Myeloid'
temp[ temp ==12 ] = 'Myeloid'
temp[ temp ==13 ] = 'Epithelial'
temp[ temp ==14 ] = 'Myeloid'
temp[ temp ==15 ] = 'Myeloid'
temp[ temp ==16 ] = 'Lymphoid'
temp[ temp ==17 ] = 'Myeloid'
temp[ temp ==18 ] = 'Myeloid'
temp[ temp ==19 ] = 'Myeloid'
temp[ temp ==20 ] = 'Myeloid'
temp[ temp ==21 ] = 'Myeloid'
temp[ temp ==22 ] = 'Fibroblast'
temp[ temp ==23 ] = 'Myeloid'
temp[ temp ==24 ] = 'Myeloid'
temp[ temp ==25 ] = 'Low Quality'
temp[ temp ==26 ] = 'Low Quality'
temp[ temp ==27 ] = 'Low Quality'
obj$cellType = temp
Idents(obj) = obj$cellType
saveRDS( obj, 'sketch_v5_all/cam_it2.rds', compress = F )
obj = subset( obj, cellType == 'Myeloid' )
obj = analyzeObj( obj, 'cam_it3' )

# Myeloid_it3
temp = as.character( obj$harmony.cluster.full )
temp[ temp == 0 ] = 'Neut CSF3R+'
temp[ temp == 1 ] = 'Mac C1QC+'
temp[ temp == 2 ] = 'Mast CPA3+'
temp[ temp == 3 ] = 'Mono CD14+'
temp[ temp == 4 ] = 'DC cDC2 CD1C+'
temp[ temp == 5 ] = 'Mac MARCO+'
temp[ temp == 6 ] = 'Low Quality'
temp[ temp == 7 ] = 'Mac ISG15+'
temp[ temp == 8 ] = 'Mac SPP1+'
temp[ temp == 9 ] = 'Mac C1QC+'
temp[ temp ==10 ] = 'Low Quality'
temp[ temp ==11 ] = 'Mac C1QC+' 
temp[ temp ==12 ] = 'Mono CD16+'
temp[ temp ==13 ] = 'DC mregDC LAMP3+'
temp[ temp ==14 ] = 'Mac C1QC+'
temp[ temp ==15 ] = 'Mac ISG15+'
temp[ temp ==16 ] = 'Mac C1QC+'
temp[ temp ==17 ] = 'DC cDC1 CLEC9A+'
temp[ temp ==18 ] = 'Osteoclast CTSK+'
temp[ temp ==19 ] = 'Mac C1QC+'
temp[ temp ==20 ] = 'Low Quality'
temp[ temp ==21 ] = 'Mac SPP1+'
temp[ temp ==22 ] = 'Mac C1QC+'
temp[ temp ==23 ] = 'Mac C1QC+'
temp[ temp ==24 ] = 'Low Quality'
temp[ temp ==25 ] = 'Low Quality' 
temp[ temp ==26 ] = 'Low Quality' 
temp[ temp ==27 ] = 'Low Quality' 
obj$cellType = temp
Idents(obj) = obj$cellType
saveRDS( obj, 'sketch_v5_all/cam_it3.rds', compress = F )
obj = subset( obj, cellType != 'Low Quality' )
saveRDS( obj, 'sketch_v5_all/cam_final.rds', compress = F )






makeFinalFigures( obj_caf, 'caf_final', genes = c('CD74','CXCL8','ISG15','COL4A1','LRRC15','CXCL14','PI16') )
makeFinalFigures( obj_cec, 'cec_final', genes = c('GJA5','CA4','PROX1','ESM1','ACKR1') )
makeFinalFigures( obj_cmc, 'cmc_final', genes = c('CCL2','CD248','ISG15','RERGL','WFDC1','DES') )
makeFinalFigures( obj_can, 'can_final', genes = c('MLANA','MPZ','NGFR') )
makeFinalFigures( obj_cal, 'cal_final', genes = c('RGS13','BANK1','IGHM','KLRD1','JCHAIN','IL7R','GZMB','GZMK','ISG15','CXCL13','FOXP3') )
makeFinalFigures( obj_cam, 'cam_final', genes = c('CLEC9A','CD1C','LAMP3','C1QC','ISG15','MARCO','SPP1','CPA3','VCAN','CD52','CSF3R','CTSK') )



############ CAF ############
obj = readRDS('sketch_v5_all/caf_final_joined.rds')
genes <- c("PDGFRA", "LUM", "DCN",
           "CD74","HLA-DRA", "HLA-DRB1",
           "CXCL8", "CXCL1", "CXCL3",
           "ISG15", "CXCL10", "IFIT3",
           "COL4A1","COL4A2","ITGA1",
           "LRRC15","SDC1","CTHRC1",
           "COL14A1", "IGF1", "NFIB",
           "PI16", "CFD", "CD34")
neg <- c("PECAM1", "RGS5", "EPCAM",
           "PLP1","MSLN")
makeSuppPanel( obj, 'caf', genes, neg, dotColor='black' )

############ CEC ############
obj = readRDS('sketch_v5_all/cec_final_joined.rds')
genes <- c("PECAM1", "VWF", "CDH5",
           "GJA5","HEY1", "SEMA3G",
           "CA4", "CD36","MT1E",
           #"ISG15", "CXCL10", "IFIT3",
           "PROX1","PDPN","CCL21",
           "ESM1","KDR","FLT1",
           "ACKR1", "SELE", "CLU",
           "PLVAP", "IGFBP3", "HSPG2",
           "COL4A1", "SELP", "POSTN")
neg <- c("EPCAM", "PDGFRA", "RGS5",
           "PLP1","PTPRC")
makeSuppPanel( obj, 'cec', genes, neg, dotColor = 'black' )

############ CMC ############
obj = readRDS('sketch_v5_all/cmc_final_joined.rds')
genes <- c("RGS5", "PDGFRB", "KCNJ8",
           "CSPG4","PROCR", "MYH11",
           "CCL2", "CCL19", "CCL21",
           "CD248","THY1","PLXDC1",
           "ISG15","CXCL10","IFIT3",
           "RERGL", "NET1", "MT1A",
           "WFDC1", "LTBP1", "SLIT3",
           "DES", "ACTG2", "PRUNE2")
neg <- c("EPCAM", "PDGFRA", "PECAM1",
           "PLP1","PTPRC")
makeSuppPanel( obj, 'cmc', genes, neg, dotColor = 'black' )

############ CAN ############
obj = readRDS('sketch_v5_all/can_final_joined.rds')
genes <- c(#"PLP1", "GPM6B", "S100B",
           "MLANA","TYRP1", "PMEL",
           "MPZ", "NRXN1", "NTM",
           "NGFR","TNC","GBP1")
neg <- c("EPCAM", "PDGFRA", "PECAM1",
           "RGS5","PTPRC")
makeSuppPanel( obj, 'can', genes, neg, dotColor = 'black', rasterFeaturePlot = F )

############ CAL ############
obj = readRDS('sketch_v5_all/cal_final_joined.rds')
genes <- c("MS4A1", "CD79A","HLA-DRA", # All B cell pops
           "RGS13","BANK1", "IGHM",    # Specific B cell pops
           "KLRD1","GNLY","NKG7",      # NK
           "JCHAIN", "MZB1", "DERL3",  # Plasma
           "CD4","CD2","CCR7",         # CD4 all
           "IL7R", "CXCL13", "FOXP3",  # CD4 subpops
           "CD8A", "CD8B", "CCL5",     # CD8 all
           "GZMB", "GZMK", "ISG15")    # CD8 subpops
neg <- c("EPCAM", "PDGFRA", "PECAM1",
           "RGS5","PLP1")
makeSuppPanel( obj, 'cal', genes, neg, dotColor = 'black' )

############ CAM ############
obj = readRDS('sketch_v5_all/cam_final.rds')
genes <- c("HLA-DRA", "FLT3", "ITGAX", # DC all
           "CLEC9A","CD1C", "LAMP3",   # DC subpops
           "CD68", "ITGAM", "C1QC",    # Mac all (1-2), subpop (3)
           "ISG15","MARCO","SPP1",     # Mac subpops
           "FCN1","VCAN","CD52",       # Mono all (1), subpop (2-3)
           "CPA3", "CSF3R", "CTSK")    # Mast, Neut, Osteo
neg <- c("EPCAM", "PDGFRA", "PECAM1",
           "RGS5","PLP1")
makeSuppPanel( obj, 'cam', genes, neg, dotColor = 'black' )




