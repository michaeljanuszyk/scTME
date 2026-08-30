# ==============================================================================
# Script: tmeMetaAnalysis_human.R
# Description: Full pipeline for integrating and analyzing human scTME datasets.
# ==============================================================================

# ------------------------------------------------------------------------------
# Required Libraries & Core Functions
# ------------------------------------------------------------------------------
library(Seurat)
library(SeuratWrappers)
library(ggplot2)
library(stringr)
library(grid)
library(gridExtra)
library(dplyr)
library(BPCells)
library(Matrix)
library(zellkonverter)
library(SingleCellExperiment)
library(reticulate)

# Source core analytical functions (Ensure this path is correct for your repo structure)
source("../tmeCore.R")

options(Seurat.object.assay.version = "v3")

# ------------------------------------------------------------------------------
# Configuration (Paths and Directories)
# ------------------------------------------------------------------------------
DIR_INPUT     <- "readyForSeurat"
DIR_QC        <- "01_qc"
DIR_PROCESSED <- "processed"
DIR_SKETCH    <- "sketch_v5_all"
DIR_H5AD      <- file.path(DIR_SKETCH, "H5AD")
DIR_BP        <- "BP"
DIR_BP_OBJ    <- "BP_object"

# Define the relative path to the reference data directory
DIR_REF       <- file.path("..", "reference")

# Define full paths for the specific CSVs
FILE_SYNONYMS <- file.path(DIR_REF, "GRCh38.p14_geneNames_geneSynonyms.csv")
FILE_METADATA <- file.path(DIR_REF, "human.meta.data.csv")

# Ensure primary output directories exist
invisible(lapply(c(DIR_QC, DIR_PROCESSED, DIR_SKETCH, DIR_H5AD, DIR_BP, DIR_BP_OBJ), function(d) {
  if (!dir.exists(d)) dir.create(d, recursive = TRUE, showWarnings = FALSE)
}))


# ==============================================================================
# 1. Data Import and QC
# ==============================================================================

obj_list <- list.files(DIR_INPUT, pattern = "*object.rds", recursive = TRUE)
obj_list <- substring(obj_list, 1, nchar(obj_list) - 11)
write.csv(obj_list, file.path(DIR_QC, "obj_list_01_preQC.csv"), row.names = FALSE)

pb <- txtProgressBar(min = 1, max = length(obj_list), style = 3)
count <- 0

for (obj_name in obj_list) {
  count <- count + 1
  setTxtProgressBar(pb, count)
  # Standard assignment replacing eval(call(...))
  assign(obj_name, readRDS(file.path(DIR_INPUT, paste0(obj_name, ".object.rds"))))
}

# Confirm col = cells; row = features
out <- data.frame(matrix(ncol = 4, nrow = 0))
for (i in seq_along(obj_list)) {
  print(i)
  obj_name <- obj_list[i]
  x <- get(obj_name)
  x[['RNA']] <- as(x[['RNA']], Class = 'Assay')
  
  dim1_1 <- x@assays[["RNA"]]@counts@Dimnames[[1]][1]
  dim1_10k <- x@assays[["RNA"]]@counts@Dimnames[[1]][10000]
  dim2_1 <- x@assays[["RNA"]]@counts@Dimnames[[2]][1]
  
  out[nrow(out) + 1, ] <- c(obj_name, dim1_1, dim1_10k, dim2_1) 
}
write.csv(out, file.path(DIR_QC, "head.csv"), row.names = FALSE)

# Confirm data has no non-integer values
out_counts <- data.frame(matrix(ncol = 2, nrow = 0))
for (i in seq_along(obj_list)) {
  print(i)
  obj_name <- obj_list[i]
  x <- get(obj_name)
  is_integer <- any(GetAssayData(object = x, slot = "counts") %% 1 != 0)
  out_counts[nrow(out_counts) + 1, ] <- c(obj_name, is_integer)
}
write.csv(out_counts, file.path(DIR_QC, "max_counts.csv"), row.names = FALSE)

# Get assay name
sink(file = file.path(DIR_QC, "assays.txt"))
for (obj_name in obj_list) {
  x <- get(obj_name)
  print(obj_name)
  show(x)
}
sink(file = NULL)

# QC Execution
for (i in seq_along(obj_list)) {
  print(i)
  objName <- obj_list[i]
  setTxtProgressBar(pb, i)
  obj <- get(objName)
  DefaultAssay(obj) <- 'RNA'
  
  if (!('nFeature_RNA' %in% colnames(obj@meta.data))) {
    counts <- obj[['RNA']]$counts
    obj$nCount_RNA <- colSums(counts)
    obj$nFeature_RNA <- colSums(counts > 0)
  }
  
  if (max(obj$nFeature_RNA) <= 200) next
  obj <- subset(obj, subset = nFeature_RNA > 200)
  if (ncol(obj) <= 1) next
  
  obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^MT-")
  VlnPlot(obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
  ggsave(file.path(DIR_QC, paste0(gsub("/", "_", objName), "_qc-01-pre-filt.jpg")))
  
  obj <- subset(obj, subset = nFeature_RNA > 200 & percent.mt < 20)
  if (ncol(obj) <= 1) next
  
  VlnPlot(obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
  ggsave(file.path(DIR_QC, paste0(gsub("/", "_", objName), "_qc-02-post-filt.jpg")))
  
  assign(objName, obj)
  
  dir_path <- file.path(DIR_PROCESSED, strsplit(objName, "/")[[1]][1])
  if (!dir.exists(dir_path)) dir.create(dir_path, recursive = TRUE, showWarnings = FALSE)
  saveRDS(obj, file.path(DIR_PROCESSED, paste0(objName, "_01_qc.rds")), compress = FALSE)
}

# Output number of cells in file
out_cells <- data.frame(matrix(ncol = 2, nrow = 0))
for (obj_name in obj_list) {
  x <- get(obj_name)
  out_cells[nrow(out_cells) + 1, ] <- c(obj_name, ncol(x))
}
write.csv(out_cells, file.path(DIR_QC, "cell_number_01_qc.csv"), row.names = FALSE)


# ==============================================================================
# 2. Pre-processing and Gene Standardization
# ==============================================================================

obj_list <- as.vector(read.csv(file.path(DIR_QC, "obj_list_01_preQC.csv"))$x)
pb <- txtProgressBar(min = 1, max = length(obj_list), style = 3)
count <- 0

# Format individual post-QC objects (streamline metadata, rename gene names)
for (i in seq_along(obj_list)) {
  print(i)
  count <- count + 1
  setTxtProgressBar(pb, count)
  obj <- obj_list[i]
  
  file_path <- file.path(DIR_PROCESSED, paste0(obj, "_01_qc.rds"))
  if (!file.exists(file_path)) next
  
  x <- readRDS(file_path)
  if (ncol(x) <= 2) next
  
  x <- DietSeurat(x)
  x@meta.data <- x@meta.data %>%
    select(any_of(c("orig.ident", "nCount_RNA", "nFeature_RNA", "study",
                    "organ", "cancer", "model", "sorting", "site", "percent.mt")))
  x <- RenameCells(x, new.names = paste(x$orig.ident, Cells(x), sep = "_"))
  
  # Deal with duplicate gene names and non-standard gene names
  counts.mat <- GetAssayData(x, slot = 'counts')
  geneNames <- rownames(counts.mat)
  synToName <- read.csv(FILE_SYNONYMS, header = TRUE)
  
  for (j in seq_along(geneNames)) {
    name <- geneNames[j]
    if (name %in% synToName$Gene.name) {
      next
    } else if (name %in% synToName$Gene.Synonym) {
      geneNames[j] <- synToName$Gene.name[which(synToName$Gene.Synonym == name)[1]]
    }
  }
  rownames(counts.mat) <- geneNames
  
  # Sum duplicate rows
  geneNamesDup <- unique(geneNames[which(duplicated(geneNames))])
  if (length(geneNamesDup) > 0) {
    counts.mat <- convert_matrix_type(counts.mat, type = "uint32_t")
    dup.counts.mat <- as(matrix(0, nrow = length(geneNamesDup), ncol = ncol(counts.mat)), Class = "dgCMatrix")
    colnames(dup.counts.mat) <- colnames(counts.mat)
    rownames(dup.counts.mat) <- geneNamesDup
    
    for (j in seq_along(geneNamesDup)) {
      gene <- geneNamesDup[j]
      rows <- which(rownames(counts.mat) == gene)
      temp <- as(counts.mat[rows, ], Class = "dgCMatrix")
      temp <- as(t(colSums(temp)), Class = "dgCMatrix")
      dup.counts.mat[j, ] <- temp
    }
    duprows <- which(duplicated(geneNames) | duplicated(geneNames, fromLast = TRUE))
    counts.mat <- rbind(as(counts.mat[-duprows, ], Class = "dgCMatrix"), dup.counts.mat)
  }
  
  x <- CreateSeuratObject(counts = counts.mat, meta.data = x@meta.data)
  saveRDS(x, file.path(DIR_PROCESSED, paste0(obj, "_01b_correctGeneName.rds")), compress = FALSE)
}


# ==============================================================================
# 3. Doublet Detection (Scrublet via Reticulate)
# ==============================================================================

# Export sample files to Python format
obj_list <- list.files(DIR_PROCESSED, pattern = "*_01b_correctGeneName.rds", recursive = TRUE)
for (i in seq_along(obj_list)) {
  print(i)
  obj <- readRDS(file.path(DIR_PROCESSED, obj_list[i]))
  sce <- as.SingleCellExperiment(obj)
  writeH5AD(sce, file.path(DIR_H5AD, paste0("forScrublet_", unique(obj$study), '_', unique(obj$orig.ident), '.h5ad')))
}

# Write python scrublet script to file securely
py_script <- "
import scanpy as sc
import scib
import anndata2ri
import anndata
import os
import scrublet as scr
import pickle
from rpy2.robjects import r
from rpy2.robjects.conversion import localconverter

anndata2ri.activate()
os.chdir('sketch_v5_all/H5AD')
dbList = []
fileList = [f for f in os.listdir() if f.endswith('.h5ad')]
badCount = 0

for i, file in enumerate(fileList):
    print(i, file)
    obj = anndata.read_h5ad(file)
    if obj.n_obs <= 30:
        badCount += 1
        print(f'bad: {badCount}')
    else:
        scrub = scr.Scrublet(obj.X)
        doublet_scores, predicted_doublets = scrub.scrub_doublets()
        obj.obs['doublet_scores'] = doublet_scores
        obj.obs['predicted_doublets'] = predicted_doublets
        dbList.append(obj.obs)

with open('scrublets.pkl', 'wb') as file:
    pickle.dump(dbList, file)
"
writeLines(py_script, "run_scrublet.py")

# Execute python script via reticulate (assumes conda env is available)
message("Running scrublet in python...")
use_condaenv("scib/anndata2ri", required = FALSE)
source_python("run_scrublet.py")

# Filter for doublets using scrublet results
source_python("pickle_reader.py") # Assuming this exists locally
pickle_data <- read_pickle_file(file.path(DIR_H5AD, "scrublets.pkl"))

obj_list <- list.files(DIR_PROCESSED, pattern = "*_01b_correctGeneName.rds", recursive = TRUE)
obj_list <- substring(obj_list, 1, nchar(obj_list) - 24)

temp <- sapply(pickle_data, function(pd) paste0(as.character(unique(pd$study)), '_', as.character(unique(pd$orig.ident))))

for (i in seq_along(obj_list)) {
  print(i)
  objName <- obj_list[i]
  obj <- readRDS(file.path(DIR_PROCESSED, paste0(objName, "_01b_correctGeneName.rds")))
  denom <- ncol(obj)

  j <- which(paste0(unique(obj$study), '_', unique(obj$orig.ident)) == temp)
  if (length(j) > 0) {
    obj$scrublet_doublet_scores <- pickle_data[[j]]$doublet_scores
    obj$scrublet_predicted_doublets <- pickle_data[[j]]$predicted_doublets
    obj <- subset(obj, scrublet_predicted_doublets == FALSE)
  }
  
  num <- ncol(obj)
  print(c(100 * (denom - num) / denom, ncol(obj)))
  saveRDS(obj, file.path(DIR_PROCESSED, paste0(objName, "_01c_singlets.rds")), compress = FALSE)
}


# ==============================================================================
# 4. Standardization and BPCells Conversion
# ==============================================================================

# Standardize rows into common gene space
obj_list <- list.files(DIR_PROCESSED, pattern = "*_01c_singlets.rds", recursive = TRUE)
obj_list <- substring(obj_list, 1, nchar(obj_list) - 17)
synToName <- read.csv(FILE_SYNONYMS, header = TRUE)

keepGenes <- sort(unique(synToName$Gene.name))
keepGenes <- keepGenes[-1] # remove empty "" gene

for (i in seq_along(obj_list)) {
  print(i)
  obj <- obj_list[i]
  x <- readRDS(file.path(DIR_PROCESSED, paste0(obj, "_01c_singlets.rds")))
  if (ncol(x) <= 1) next
  
  counts.mat <- GetAssayData(x, slot = 'counts')
  counts.mat <- convert_matrix_type(counts.mat, type = "uint32_t")
  counts.mat <- counts.mat[rownames(counts.mat) %in% synToName$Gene.name, ]
  counts.mat <- as(counts.mat, Class = "dgCMatrix")
  
  mat1 <- counts.mat
  ii1 <- match(keepGenes, rownames(mat1))
  ii1[is.na(ii1)] <- nrow(mat1) + 1
  mat1 <- rbind(mat1, rep(0, ncol(mat1)))
  new1 <- mat1[ii1, ]
  rownames(new1) <- keepGenes
  
  x <- CreateSeuratObject(counts = new1, meta.data = x@meta.data)
  dir_path <- file.path(DIR_PROCESSED, strsplit(obj, "/")[[1]][1])
  if (!dir.exists(dir_path)) dir.create(dir_path, recursive = TRUE, showWarnings = FALSE)
  
  saveRDS(x, file.path(DIR_PROCESSED, paste0(obj, "_01d_stddim.rds")), compress = FALSE)
}

# Merge samples for each study
study_list <- list.files(DIR_PROCESSED, recursive = FALSE)
for (i in seq_along(study_list)) {
  print(i)
  study <- study_list[i]
  study_list_samples <- list.files(file.path(DIR_PROCESSED, study), pattern = "*_01d_stddim.rds")
  if (length(study_list_samples) == 0) next
  
  study_list_obj <- list()
  for (sample in study_list_samples) {
    x <- readRDS(file.path(DIR_PROCESSED, study, sample))
    study_list_obj <- append(study_list_obj, x)
  }
  
  if (length(study_list_obj) > 1) {
    obj <- merge(study_list_obj[[1]], study_list_obj[2:length(study_list_obj)])
    if (("data" %in% Layers(obj)) && length(Layers(obj)) == 2) {
      warning('Layer mismatch detected.')
    } else {
      obj <- JoinLayers(obj)
    }
  } else {
    obj <- study_list_obj[[1]]
  }
  
  if (ncol(obj) > 1) {
    saveRDS(obj, file.path(DIR_PROCESSED, study, paste0(study, "_01e_study.rds")), compress = FALSE)
  }
}

# Find minimal expression genes
obj_list <- list.files(DIR_PROCESSED, pattern = "*_01e_study.rds", recursive = TRUE)
cell_count <- 0
obj <- readRDS(file.path(DIR_PROCESSED, obj_list[1]))
gene_counts <- as.data.frame(matrix(0, ncol = ncol(obj), nrow = 0))
colnames(gene_counts) <- ncol(obj)

for (i in seq_along(obj_list)) {
  print(i)
  objName <- obj_list[i]
  file_path <- file.path(DIR_PROCESSED, objName)
  if (!file.exists(file_path)) next
  
  obj <- readRDS(file_path)
  counts.mat <- obj[["RNA"]]$counts
  cell_count <- cell_count + ncol(counts.mat)
  gene_counts <- rbind(gene_counts, rowSums(counts.mat > 0))
}

nCount_Feature <- colSums(gene_counts, na.rm = TRUE)
keepGenesIdx <- c(nCount_Feature >= 0.001 * cell_count)
keepGenesNames <- colnames(gene_counts)[keepGenesIdx]
write.csv(keepGenesIdx, "keepGenesIdx.csv", quote = FALSE)
write.csv(keepGenesNames, "keepGenesNames.csv", quote = FALSE, row.names = FALSE)

# Convert to BP cell space
for (i in seq_along(obj_list)) {
  print(i)
  obj <- readRDS(file.path(DIR_PROCESSED, obj_list[i]))
  counts.mat <- obj[["RNA"]]$counts[keepGenesIdx, ]
  
  base_name <- substring(obj_list[i], 1, nchar(obj_list[i]) - 14)
  bp_path <- file.path(DIR_BP, base_name)
  
  write_matrix_dir(mat = counts.mat, dir = bp_path)
  counts.mat <- open_matrix_dir(dir = bp_path)
  
  obj <- CreateSeuratObject(counts = counts.mat, meta.data = obj@meta.data)
  
  dir_path <- file.path(DIR_BP_OBJ, strsplit(obj_list[i], "/")[[1]][1])
  if (!dir.exists(dir_path)) dir.create(dir_path, recursive = TRUE, showWarnings = FALSE)
  
  saveRDS(obj, file.path(DIR_BP_OBJ, paste0(base_name, "_01f_prunedGenes.rds")), compress = FALSE)
}


# ==============================================================================
# 5. Core Integrations (Harmony, MNN, scVI)
# ==============================================================================

options(Seurat.object.assay.version = "v5")
obj_list <- list.files(DIR_BP_OBJ, pattern = "*_01f_prunedGenes.rds", recursive = TRUE)
first_elements <- sapply(strsplit(obj_list, "/"), function(x) x[1])

caf_list <- list()
cell_count <- 0

for (i in seq_along(obj_list)) {
  print(i)
  obj <- readRDS(file.path(DIR_BP_OBJ, obj_list[i]))
  obj <- NormalizeData(obj)
  obj$study <- first_elements[i]
  caf_list[[i]] <- obj
  cell_count <- cell_count + ncol(obj)
}

merged <- merge(caf_list[[1]], caf_list[2:length(caf_list)])
saveRDS(merged, file.path(DIR_SKETCH, "humanMerged.rds"), compress = FALSE)

# Attach clinical metadata
clinical <- read.csv(FILE_METADATA, header = TRUE)
clinical <- clinical[, -1]
mapping <- match(merged$orig.ident, clinical$orig.ident)

for (col in colnames(clinical)) {
  merged[[col]] <- clinical[[col]][mapping]
}

# --- Harmony Integration ---
obj <- NormalizeData(merged)
saveRDS(obj, file.path(DIR_SKETCH, "unfiltered_counts_obj_norm.rds"), compress = FALSE)
obj <- FindVariableFeatures(obj, verbose = TRUE)
saveRDS(obj, file.path(DIR_SKETCH, "unfiltered_counts_obj_findvar.rds"), compress = FALSE)

obj <- SketchData(obj, ncells = 1000, method = "LeverageScore", sketched.assay = "sketch")
saveRDS(obj, file.path(DIR_SKETCH, "unfiltered_counts_obj_findvar_sketch1000.rds"), compress = FALSE)

DefaultAssay(obj) <- "sketch"
obj <- FindVariableFeatures(obj, verbose = TRUE)
obj <- ScaleData(obj, verbose = FALSE)
obj <- RunPCA(obj, npcs = 30, verbose = TRUE)
saveRDS(obj, file.path(DIR_SKETCH, "unfiltered_counts_obj_findvar_sketch1000_pca.rds"), compress = FALSE)

obj <- IntegrateLayers(obj, method = HarmonyIntegration, orig.reduction = "pca", new.reduction = "harmony", kmeans_init_nstart = 20, kmeans_init_iter_max = 5000, verbose = TRUE)
saveRDS(obj, file.path(DIR_SKETCH, "unfiltered_counts_obj_findvar_sketch1000_harmony.rds"), compress = FALSE)

obj <- RunUMAP(obj, reduction = "harmony", dims = 1:30, reduction.name = "umap.harmony", return.model = TRUE, min.dist = 0.001)
saveRDS(obj, file.path(DIR_SKETCH, "unfiltered_counts_obj_findvar_sketch1000_harmony_umap.rds"), compress = FALSE)

obj <- FindNeighbors(obj, reduction = "harmony", dims = 1:30)
for (i in c(0.06, 0.07, 0.08, 0.09, 0.1, 0.12, 0.14, 0.16, 0.018, 0.2)) {
  obj <- FindClusters(obj, resolution = i, graph.name = 'sketch_snn')
  DimPlot(obj, raster = TRUE, label = TRUE)
  ggsave(file.path(DIR_SKETCH, paste0("unfiltered_counts_obj_findvar_sketch1000_harmony_umap_res", i, ".jpg")), width = 5, height = 5, units = "in", limitsize = FALSE)
  saveRDS(obj, file.path(DIR_SKETCH, paste0("unfiltered_counts_obj_findvar_sketch1000_harmony_umap_res", i, ".rds")), compress = FALSE)
}

# Project to full dataset
obj$seurat_clusters <- obj$sketch_snn_res.0.06
obj <- ProjectIntegration(object = obj, reduction = "harmony")
options(future.globals.maxSize = 8000 * 1024^2)
obj <- ProjectData(object = obj, sketched.reduction = "harmony.full", full.reduction = "harmony.full", umap.model = "umap.harmony", dims = 1:30, refdata = list(harmony.cluster.full = "seurat_clusters"))
saveRDS(obj, file.path(DIR_SKETCH, "unfiltered_counts_obj_findvar_sketch1000_harmony_umap_project.rds"), compress = FALSE)

# --- MNN Integration ---
obj_mnn <- readRDS(file.path(DIR_SKETCH, "unfiltered_counts_obj_findvar_sketch1000_pca.rds"))
obj_mnn <- IntegrateLayers(object = obj_mnn, method = FastMNNIntegration, orig.reduction = "pca", new.reduction = "integrated.mnn", verbose = TRUE)
saveRDS(obj_mnn, file.path(DIR_SKETCH, "unfiltered_counts_obj_findvar_sketch1000_mnn.rds"), compress = FALSE)

obj_mnn <- RunUMAP(obj_mnn, reduction = "integrated.mnn", dims = 1:30, reduction.name = "umap.mnn", return.model = TRUE, min.dist = 0.001)
saveRDS(obj_mnn, file.path(DIR_SKETCH, "unfiltered_counts_obj_findvar_sketch1000_mnn_umap.rds"), compress = FALSE)

obj_mnn <- FindNeighbors(obj_mnn, reduction = "integrated.mnn", dims = 1:30)
for (i in c(0.05, 0.06, 0.07, 0.09, 0.1, 0.15, 0.2)) {
  obj_mnn <- FindClusters(obj_mnn, resolution = i, graph.name = 'RNA_snn')
  DimPlot(obj_mnn, raster = TRUE, label = TRUE, reduction = 'umap.mnn')
  ggsave(file.path(DIR_SKETCH, paste0("unfiltered_counts_obj_findvar_sketch1000_mnn_umap_res", i, ".jpg")), width = 5, height = 5, units = "in", limitsize = FALSE)
  saveRDS(obj_mnn, file.path(DIR_SKETCH, paste0("unfiltered_counts_obj_findvar_sketch1000_mnn_umap_res", i, ".rds")), compress = FALSE)
}

obj_mnn <- readRDS(file.path(DIR_SKETCH, "unfiltered_counts_obj_findvar_sketch1000_mnn_umap_res0.5.rds"))
obj_mnn$seurat_clusters <- obj_mnn$RNA_snn_res.0.3
Idents(obj_mnn) <- obj_mnn$seurat_clusters 
obj_mnn <- ProjectIntegration(object = obj_mnn, reduction = "integrated.mnn")
obj_mnn <- ProjectData(object = obj_mnn, sketched.reduction = "integrated.mnn.full", full.reduction = "integrated.mnn.full", umap.model = "umap.mnn", dims = 1:30, refdata = list(mnn.cluster.full = "seurat_clusters"))
saveRDS(obj_mnn, file.path(DIR_SKETCH, "unfiltered_counts_obj_findvar_sketch1000_mnn_umap_project.rds"), compress = FALSE)

# --- scVI Integration ---
Sys.setenv(RETICULATE_PYTHON = "scvi-env/bin/python")
obj_scvi <- readRDS(file.path(DIR_SKETCH, "unfiltered_counts_obj_findvar_sketch1000_pca.rds"))
obj_scvi <- IntegrateLayers(object = obj_scvi, method = scVIIntegration, orig.reduction = "pca", new.reduction = "integrated.scvi", conda_env = '~/.conda/envs/scvi-env/', verbose = TRUE)
saveRDS(obj_scvi, file.path(DIR_SKETCH, "unfiltered_counts_obj_findvar_sketch1000_scvi.rds"), compress = FALSE)

obj_scvi <- RunUMAP(obj_scvi, reduction = "integrated.scvi", dims = 1:30, reduction.name = "umap.scvi", return.model = TRUE, min.dist = 0.001)
ggsave(file.path(DIR_SKETCH, "unfiltered_obj_sketch1000_scvi_featureplot.jpg"), width = 8, height = 10, units = "in", limitsize = FALSE)
saveRDS(obj_scvi, file.path(DIR_SKETCH, "unfiltered_counts_obj_findvar_sketch1000_scvi_umap.rds"), compress = FALSE)

obj_scvi <- FindNeighbors(obj_scvi, reduction = "integrated.scvi", dims = 1:30)
for (i in c(0.05, 0.06, 0.07, 0.09, 0.1, 0.15, 0.2)) {
  obj_scvi <- FindClusters(obj_scvi, resolution = i, graph.name = 'RNA_snn')
  DimPlot(obj_scvi, raster = TRUE, reduction = 'umap.scvi', label = TRUE, group.by = paste0('RNA_snn_res.', i))
  ggsave(file.path(DIR_SKETCH, paste0("unfiltered_counts_obj_findvar_sketch1000_scvi_umap_res", i, ".jpg")), width = 5, height = 5, units = "in", limitsize = FALSE)
  saveRDS(obj_scvi, file.path(DIR_SKETCH, paste0("unfiltered_counts_obj_findvar_sketch1000_scvi_umap_res", i, ".rds")), compress = FALSE)
}

obj_scvi$seurat_clusters <- obj_scvi$RNA_snn_res.0.2
Idents(obj_scvi) <- obj_scvi$seurat_clusters
obj_scvi <- ProjectIntegration(object = obj_scvi, reduction = "integrated.scvi")
obj_scvi <- ProjectData(object = obj_scvi, sketched.reduction = "integrated.scvi.full", full.reduction = "integrated.scvi.full", umap.model = "umap.scvi", dims = 1:30, refdata = list(scvi.cluster.full = "seurat_clusters"))
saveRDS(obj_scvi, file.path(DIR_SKETCH, "unfiltered_counts_obj_findvar_sketch1000_scvi_umap_project.rds"), compress = FALSE)


# ==============================================================================
# 6. Cluster Annotations (Dictionary Mapping)
# ==============================================================================

# --- Harmony Annotations ---
harmony_map <- c(
  "0" = "Epithelial", "1" = "Lymphoid", "2" = "Myeloid", "3" = "Fibroblast",
  "4" = "Endothelial", "5" = "Cycling", "6" = "Mural", "7" = "Lymphoid",
  "8" = "Lymphoid", "9" = "Nerve", "10" = "Myeloid"
)
obj$cellType <- unname(harmony_map[as.character(obj$harmony.cluster.full)])
DefaultAssay(obj) <- obj$cellType
saveRDS(obj, file.path(DIR_SKETCH, 'allCells_harmony.rds'), compress = FALSE)

# --- MNN Annotations ---
mnn_map <- c(
  "0" = "Epithelial", "1" = "Lymphoid", "2" = "Myeloid", "3" = "Fibroblast",
  "4" = "Epithelial", "5" = "Endothelial", "6" = "Epithelial", "7" = "Epithelial",
  "8" = "Mural", "9" = "Nerve", "10" = "Lymphoid", "11" = "Epithelial",
  "12" = "Lymphoid", "13" = "Myeloid", "14" = "Cycling", "15" = "Myeloid"
)
obj_mnn$cellType <- unname(mnn_map[as.character(obj_mnn$mnn.cluster.full)])
Idents(obj_mnn) <- obj_mnn$cellType
saveRDS(obj_mnn, file.path(DIR_SKETCH, 'allCells_mnn.rds'), compress = FALSE)

# --- scVI Annotations ---
scvi_map <- c(
  "0" = "Myeloid", "1" = "Fibroblast", "2" = "Lymphoid", "3" = "Epithelial",
  "4" = "Endothelial", "5" = "Epithelial", "6" = "Epithelial", "7" = "Mural",
  "8" = "Lymphoid", "9" = "Epithelial", "10" = "Lymphoid", "11" = "Epithelial",
  "12" = "Epithelial", "13" = "Nerve", "14" = "Lymphoid", "15" = "Epithelial",
  "16" = "Lymphoid", "17" = "Myeloid", "18" = "Epithelial", "19" = "Epithelial",
  "20" = "Fibroblast", "21" = "Endothelial", "22" = "Epithelial", "23" = "Epithelial",
  "24" = "Epithelial", "25" = "Lymphoid", "26" = "Epithelial", "27" = "Epithelial",
  "28" = "Epithelial", "29" = "Epithelial", "30" = "Lymphoid", "31" = "Myeloid"
)
obj_scvi$cellType <- unname(scvi_map[as.character(obj_scvi$scvi.cluster.full)])
Idents(obj_scvi) <- obj_scvi$cellType
saveRDS(obj_scvi, file.path(DIR_SKETCH, 'allCells_scvi.rds'), compress = FALSE)


# ==============================================================================
# 7. Sub-Clustering and Fine Annotations
# ==============================================================================

# Core Integration Splits
obj_caf <- combineIntegration('Fibroblast')
obj_caf <- analyzeObj(obj_caf, 'caf_it1')

obj_cec <- combineIntegration('Endothelial')
obj_cec <- analyzeObj(obj_cec, 'cec_it1')

obj_cmc <- combineIntegration('Mural')
obj_cmc <- analyzeObj(obj_cmc, 'cmc_it1')

obj_can <- combineIntegration('Nerve')
obj_can <- analyzeObj(obj_can, 'can_it1')

obj_cal <- combineIntegration('Lymphoid')
obj_cal <- analyzeObj(obj_cal, 'cal_it1')

obj_cam <- combineIntegration('Myeloid')
obj_cam <- analyzeObj(obj_cam, 'cam_it1')

# --- Fibroblast Sub-Clustering Iterations ---
caf_it1_map <- c("0"="Epithelial", "1"="Fibroblast", "2"="Fibroblast", "3"="Fibroblast", "4"="Fibroblast", "5"="Lymphoid", "6"="Myeloid", "7"="Endothelial", "8"="Lymphoid", "9"="Nerve", "10"="Epithelial", "11"="Low Quality", "12"="Nerve", "13"="Low Quality", "14"="Low Quality")
obj_caf$cellType <- unname(caf_it1_map[as.character(obj_caf$harmony.cluster.full)])
Idents(obj_caf) <- obj_caf$cellType
saveRDS(obj_caf, file.path(DIR_SKETCH, 'caf_it1.rds'), compress = FALSE)

obj_caf <- buildFromSamples(target = 'Fibroblast', iteration = 2, useSketch = TRUE, addIdents = FALSE)
obj_caf <- iterativeFilter(obj = obj_caf, name = 'caf_it2', skipToScale = TRUE)
obj_caf <- projectSketchedData(obj = obj_caf, name = 'caf_it2', sketchedLabels = obj_caf$cellType)

caf_it2_map <- c("0"="Fibroblast", "1"="Low Quality", "2"="Fibroblast", "3"="Mural", "4"="Fibroblast", "5"="Low Quality", "6"="Fibroblast", "7"="Fibroblast", "8"="Cycling", "9"="Low Quality", "10"="Low Quality", "11"="Low Quality", "12"="Low Quality")
obj_caf$cellType <- unname(caf_it2_map[as.character(obj_caf$harmony.cluster.full)])
Idents(obj_caf) <- obj_caf$cellType
saveRDS(obj_caf, file.path(DIR_SKETCH, 'caf_it2.rds'), compress = FALSE)

obj_caf <- buildFromSamples(target = 'Fibroblast', iteration = 3, useSketch = TRUE, addIdents = FALSE)
obj_caf <- iterativeFilter(obj = obj_caf, name = 'caf_it3', skipToScale = TRUE)
obj_caf <- projectSketchedData(obj = obj_caf, name = 'caf_it3', sketchedLabels = obj_caf$cellType)

caf_it3_map <- c("0"="Fibroblast", "1"="Fibroblast", "2"="Low Quality", "3"="Fibroblast", "4"="Fibroblast", "5"="Fibroblast", "6"="Fibroblast", "7"="Fibroblast", "8"="Fibroblast", "9"="Fibroblast", "10"="Fibroblast", "11"="Fibroblast", "12"="Fibroblast", "13"="Low Quality", "14"="Fibroblast", "15"="Low Quality", "16"="Low Quality")
obj_caf$cellType <- unname(caf_it3_map[as.character(obj_caf$harmony.cluster.full)])
Idents(obj_caf) <- obj_caf$cellType
saveRDS(obj_caf, file.path(DIR_SKETCH, 'caf_it3.rds'), compress = FALSE)

obj_caf <- subset(obj_caf, cellType == 'Fibroblast')
obj_caf <- analyzeObj(obj_caf, 'caf_it4')

caf_it4_map <- c("0"="Fibroblast", "1"="Fibroblast", "2"="Low Quality", "3"="Low Quality", "4"="Fibroblast", "5"="Fibroblast", "6"="Fibroblast", "7"="Fibroblast", "8"="Fibroblast", "9"="Fibroblast", "10"="Fibroblast", "11"="Low Quality", "12"="Low Quality", "13"="Low Quality", "14"="Low Quality")
obj_caf$cellType <- unname(caf_it4_map[as.character(obj_caf$harmony.cluster.full)])
Idents(obj_caf) <- obj_caf$cellType
saveRDS(obj_caf, file.path(DIR_SKETCH, 'caf_it4.rds'), compress = FALSE)

obj_caf <- subset(obj_caf, cellType == 'Fibroblast')
obj_caf <- analyzeObj(obj_caf, 'caf_it5')

caf_it5_map <- c("0"="Fibroblast", "1"="Fibroblast", "2"="Fibroblast", "3"="Fibroblast", "4"="Fibroblast", "5"="Tumor", "6"="Fibroblast", "7"="Fibroblast", "8"="Fibroblast", "9"="Fibroblast", "10"="Fibroblast", "11"="Low Quality")
obj_caf$cellType <- unname(caf_it5_map[as.character(obj_caf$harmony.cluster.full)])
Idents(obj_caf) <- obj_caf$cellType
saveRDS(obj_caf, file.path(DIR_SKETCH, 'caf_it5.rds'), compress = FALSE)

obj_caf <- subset(obj_caf, cellType == 'Fibroblast')
obj_caf <- analyzeObj(obj_caf, 'caf_it6')

caf_it6_map <- c("0"="mCAF LRRC15+", "1"="ssCAF PI16+", "2"="iCAF CXCL8+", "3"="mCAF LRRC15+", "4"="mCAF COL4A1+", "5"="ssCAF PI16+", "6"="ssCAF PI16+", "7"="Low Quality", "8"="iCAF CXCL8+", "9"="mCAF COL4A1+", "10"="ssCAF CXCL14+", "11"="iCAF ISG15+", "12"="apCAF CD74+", "13"="ssCAF CXCL14+", "14"="ssCAF PI16+", "15"="mCAF LRRC15+", "16"="ssCAF CXCL14+", "17"="mCAF LRRC15+", "18"="Low Quality", "19"="mCAF COL4A1+", "20"="mCAF LRRC15+", "21"="mCAF COL4A1+", "22"="ssCAF CXCL14+")
obj_caf$cellType <- unname(caf_it6_map[as.character(obj_caf$harmony.cluster.full)])
Idents(obj_caf) <- obj_caf$cellType
saveRDS(obj_caf, file.path(DIR_SKETCH, 'caf_it6.rds'), compress = FALSE)

obj_caf <- subset(obj_caf, cellType != 'Low Quality')
saveRDS(obj_caf, file.path(DIR_SKETCH, 'caf_final.rds'), compress = FALSE)

# --- Endothelial Sub-Clustering ---
cec_it1_map <- c("0"="Endothelial", "1"="Epithelial", "2"="Lymphoid", "3"="Myeloid", "4"="Fibroblast", "5"="Mural", "6"="Cycling", "7"="Endothelial", "8"="Lymphoid", "9"="Lymphoid", "10"="Myeloid")
obj_cec$cellType <- unname(cec_it1_map[as.character(obj_cec$harmony.cluster.full)])
Idents(obj_cec) <- obj_cec$cellType
saveRDS(obj_cec, file.path(DIR_SKETCH, 'cec_it1.rds'), compress = FALSE)

obj_cec <- subset(obj_cec, cellType == 'Endothelial')
obj_cec <- analyzeObj(obj_cec, 'cec_it2')

cec_it2_map <- c("0"="Endothelial", "1"="Low Quality", "2"="Endothelial", "3"="Endothelial", "4"="Endothelial", "5"="Endothelial", "6"="Epithelial", "7"="Fibroblast", "8"="Endothelial", "9"="Endothelial", "10"="Endothelial", "11"="Low Quality", "12"="Low Quality", "13"="Immune", "14"="Epithelial")
obj_cec$cellType <- unname(cec_it2_map[as.character(obj_cec$RNA_snn_res.0.2)])
saveRDS(obj_cec, file.path(DIR_SKETCH, 'cec_it2.rds'), compress = FALSE)

obj_cec <- subset(obj_cec, cellType == 'Endothelial')
obj_cec <- analyzeObj(obj_cec, 'cec_it3')

cec_it3_map <- c("0"="Endothelial", "1"="Endothelial", "2"="Endothelial", "3"="Endothelial", "4"="Endothelial", "5"="Low Quality", "6"="Low Quality", "7"="Low Quality", "8"="Low Quality", "9"="Low Quality", "10"="Low Quality", "11"="Low Quality", "12"="Low Quality")
obj_cec$cellType <- unname(cec_it3_map[as.character(obj_cec$RNA_snn_res.0.16)])
saveRDS(obj_cec, file.path(DIR_SKETCH, 'cec_it3.rds'), compress = FALSE)

obj_cec <- subset(obj_cec, cellType == 'Endothelial')
obj_cec <- analyzeObj(obj_cec, 'cec_it4')

cec_it4_map <- c("0"="Vein ACKR1+", "1"="Tip ESM1+", "2"="Artery GJA5+", "3"="Lymphatic PROX1+", "4"="Capillary CA4+")
obj_cec$cellType <- unname(cec_it4_map[as.character(obj_cec$RNA_snn_res.0.09)])
Idents(obj_cec) <- obj_cec$cellType
Idents(obj_cec) <- factor(x = Idents(obj_cec), levels = sort(levels(obj_cec), decreasing = FALSE))
saveRDS(obj_cec, file.path(DIR_SKETCH, 'cec_it4.rds'), compress = FALSE)
saveRDS(obj_cec, file.path(DIR_SKETCH, 'cec_final.rds'), compress = FALSE)

# --- Mural Sub-Clustering ---
cmc_it1_map <- c("0"="Low Quality", "1"="Mural", "2"="Mural", "3"="Fibroblast", "4"="Mural", "5"="Lymphoid", "6"="Epithelial", "7"="Endothelial", "8"="Myeloid", "9"="Epithelial", "10"="Mural", "11"="Mural", "12"="Fibroblast", "13"="Mural", "14"="Mural", "15"="Lymphoid", "16"="Endothelial", "17"="Mural", "18"="Low Quality", "19"="Low Quality", "20"="Low Quality", "21"="Low Quality", "22"="Low Quality", "23"="Low Quality", "24"="Low Quality", "25"="Low Quality", "26"="Low Quality", "27"="Low Quality", "28"="Low Quality", "29"="Low Quality", "30"="Low Quality")
obj_cmc$cellType <- unname(cmc_it1_map[as.character(obj_cmc$harmony.cluster.full)])
Idents(obj_cmc) <- obj_cmc$cellType
saveRDS(obj_cmc, file.path(DIR_SKETCH, 'cmc_it1.rds'), compress = FALSE)

obj_cmc <- subset(obj_cmc, cellType == 'Mural')
obj_cmc <- analyzeObj(obj_cmc, 'cmc_it2')

cmc_it2_map <- c("0"="Mural", "1"="Mural", "2"="Mural", "3"="Mural", "4"="Mural", "5"="Cycling", "6"="Lymphoid", "7"="Nerve", "8"="Cycling", "9"="Low Quality", "10"="Low Quality")
obj_cmc$cellType <- unname(cmc_it2_map[as.character(obj_cmc$harmony.cluster.full)])
Idents(obj_cmc) <- obj_cmc$cellType
saveRDS(obj_cmc, file.path(DIR_SKETCH, 'cmc_it2.rds'), compress = FALSE)

obj_cmc <- subset(obj_cmc, cellType == 'Mural')
obj_cmc <- analyzeObj(obj_cmc, 'cmc_it3')

cmc_it3_map <- c("0"="Mural", "1"="Mural", "2"="Mural", "3"="Mural", "4"="Mural", "5"="Mural", "6"="Mural", "7"="Mural", "8"="Mural", "9"="Mural", "10"="Low Quality", "11"="Low Quality", "12"="Mural", "13"="Mural", "14"="Mural", "15"="Low Quality", "16"="Mural", "17"="Mural", "18"="Low Quality", "19"="Low Quality", "20"="Low Quality", "21"="Low Quality", "22"="Low Quality")
obj_cmc$cellType <- unname(cmc_it3_map[as.character(obj_cmc$harmony.cluster.full)])
Idents(obj_cmc) <- obj_cmc$cellType
saveRDS(obj_cmc, file.path(DIR_SKETCH, 'cmc_it3.rds'), compress = FALSE)

obj_cmc <- subset(obj_cmc, cellType == 'Mural')
obj_cmc <- analyzeObj(obj_cmc, 'cmc_it4')

cmc_it4_map <- c("0"="Mural", "1"="Mural", "2"="Low Quality", "3"="Mural", "4"="Mural", "5"="Mural", "6"="Mural", "7"="Mural", "8"="Mural", "9"="Mural", "10"="Mural", "11"="Mural", "12"="Mural", "13"="Mural", "14"="Mural", "15"="Mural", "16"="Low Quality")
obj_cmc$cellType <- unname(cmc_it4_map[as.character(obj_cmc$harmony.cluster.full)])
Idents(obj_cmc) <- obj_cmc$cellType
saveRDS(obj_cmc, file.path(DIR_SKETCH, 'cmc_it4.rds'), compress = FALSE)

obj_cmc <- subset(obj_cmc, cellType == 'Mural')
obj_cmc <- analyzeObj(obj_cmc, 'cmc_it5')

cmc_it5_map <- c("0"="Pericyte CD248+", "1"="SMC vascular RERGL+", "2"="Pericyte CD248+", "3"="SMC vascular RERGL+", "4"="SMC vascular RERGL+", "5"="SMC vascular RERGL+", "6"="Pericyte CD248+", "7"="Pericyte CCL2+", "8"="Pericyte CD248+", "9"="SMC vascular WFDC1+", "10"="Pericyte CD248+", "11"="Pericyte ISG15+", "12"="Pericyte CD248+", "13"="Pericyte CCL2+", "14"="Pericyte CCL2+", "15"="Pericyte CD248+", "16"="SMC vascular WFDC1+", "17"="Pericyte CD248+", "18"="Pericyte CCL2+", "19"="SMC vascular WFDC1+", "20"="Low Quality", "21"="SMC visceral DES+", "22"="SMC vascular WFDC1+", "23"="SMC visceral DES+", "24"="Pericyte CCL2+", "25"="Pericyte CCL2+", "26"="Pericyte CCL2+", "27"="Pericyte CD248+", "28"="Pericyte CCL2+", "29"="SMC vascular WFDC1+", "30"="SMC vascular WFDC1+", "31"="Low Quality", "32"="Low Quality", "33"="Low Quality")
obj_cmc$cellType <- unname(cmc_it5_map[as.character(obj_cmc$harmony.cluster.full)])
Idents(obj_cmc) <- obj_cmc$cellType
saveRDS(obj_cmc, file.path(DIR_SKETCH, 'cmc_it5.rds'), compress = FALSE)

obj_cmc <- subset(obj_cmc, cellType != 'Low Quality')
saveRDS(obj_cmc, file.path(DIR_SKETCH, 'cmc_final.rds'), compress = FALSE)

# --- Nerve Sub-Clustering ---
can_it1_map <- c("0"="Epithelial", "1"="Epithelial", "2"="Nerve", "3"="Lymphoid", "4"="Fibroblast", "5"="Myeloid", "6"="Fibroblast", "7"="Endothelial", "8"="Epithelial", "9"="Cycling", "10"="Nerve", "11"="Mural", "12"="Nerve", "13"="Fibroblast", "14"="Low Quality", "15"="Mural", "16"="Mural", "17"="Epithelial", "18"="Epithelial", "19"="Low Quality", "20"="Low Quality")
obj_can$cellType <- unname(can_it1_map[as.character(obj_can$harmony.cluster.full)])
Idents(obj_can) <- obj_can$cellType
saveRDS(obj_can, file.path(DIR_SKETCH, 'can_it1.rds'), compress = FALSE)

obj_can <- subset(obj_can, cellType == 'Nerve')
obj_can <- analyzeObj(obj_can, 'can_it2')

can_it2_map <- c("0"="Nerve", "1"="Nerve", "2"="Low Quality", "3"="Nerve", "4"="Epithelial", "5"="Nerve", "6"="Lymphoid", "7"="Nerve", "8"="Nerve", "9"="Myeloid", "10"="Nerve", "11"="Nerve", "12"="Myeloid", "13"="Myeloid", "14"="Low Quality", "15"="Nerve", "16"="Low Quality", "17"="Fibroblast", "18"="Low Quality", "19"="Low Quality", "20"="Low Quality", "21"="Low Quality", "22"="Low Quality", "23"="Low Quality", "24"="Low Quality", "25"="Low Quality", "26"="Low Quality", "27"="Low Quality", "28"="Low Quality", "29"="Low Quality", "30"="Low Quality")
obj_can$cellType <- unname(can_it2_map[as.character(obj_can$RNA_snn_res.0.2)])
Idents(obj_can) <- obj_can$cellType
saveRDS(obj_can, file.path(DIR_SKETCH, 'can_it2.rds'), compress = FALSE)

obj_can <- subset(obj_can, cellType == 'Nerve')
obj_can <- analyzeObj(obj_can, 'can_it3')

can_it3_map <- c("0"="Low Quality", "1"="Low Quality", "2"="Low Quality", "3"="Nerve", "4"="Low Quality", "5"="Nerve", "6"="Lymphoid", "7"="Nerve", "8"="Myeloid", "9"="Low Quality", "10"="Nerve", "11"="Low Quality", "12"="Low Quality", "13"="Fibroblast", "14"="Fibroblast", "15"="Low Quality", "16"="Low Quality")
obj_can$cellType <- unname(can_it3_map[as.character(obj_can$RNA_snn_res.0.1)])
Idents(obj_can) <- obj_can$cellType
saveRDS(obj_can, file.path(DIR_SKETCH, 'can_it3.rds'), compress = FALSE)

obj_can <- subset(obj_can, cellType == 'Nerve')
obj_can <- analyzeObj(obj_can, 'can_it4')

can_it4_map <- c("0"="Nerve", "1"="Immune", "2"="Nerve", "3"="Nerve", "4"="Low Quality", "5"="Low Quality", "6"="Low Quality", "7"="Low Quality", "8"="Nerve", "9"="Low Quality", "10"="Low Quality", "11"="Low Quality", "12"="Nerve", "13"="Low Quality", "14"="Low Quality", "15"="Low Quality", "16"="Low Quality", "17"="Epithelial", "18"="Immune", "19"="Low Quality", "20"="Low Quality", "21"="Low Quality", "22"="Low Quality")
obj_can$cellType <- unname(can_it4_map[as.character(obj_can$RNA_snn_res.0.4)])
Idents(obj_can) <- obj_can$cellType
saveRDS(obj_can, file.path(DIR_SKETCH, 'can_it4.rds'), compress = FALSE)

obj_can <- subset(obj_can, cellType == 'Nerve')
obj_can <- analyzeObj(obj_can, 'can_it5')

can_it5_map <- c("0"="Nerve", "1"="Low Quality", "2"="Nerve", "3"="Nerve", "4"="Lymphoid", "5"="Low Quality", "6"="Nerve", "7"="Low Quality", "8"="Lymphoid", "9"="Nerve", "10"="Low Quality", "11"="Nerve", "12"="Mural", "13"="Low Quality")
obj_can$cellType <- unname(can_it5_map[as.character(obj_can$RNA_snn_res.0.4)])
Idents(obj_can) <- obj_can$cellType
saveRDS(obj_can, file.path(DIR_SKETCH, 'can_it5.rds'), compress = FALSE)

obj_can <- subset(obj_can, cellType == 'Nerve')
obj_can <- analyzeObj(obj_can, 'can_it6')

can_it6_map <- c("0"="Nerve", "1"="Nerve", "2"="Low Quality")
obj_can$cellType <- unname(can_it6_map[as.character(obj_can$RNA_snn_res.0.06)])
Idents(obj_can) <- obj_can$cellType
saveRDS(obj_can, file.path(DIR_SKETCH, 'can_it6.rds'), compress = FALSE)

obj_can <- subset(obj_can, cellType == 'Nerve')
obj_can <- analyzeObj(obj_can, 'can_it7')

can_it7_map <- c("0"="Melanocyte MLANA+", "1"="Melanocyte MLANA+", "2"="Schwann myel MPZ+", "3"="Schwann non-myel NGFR+", "4"="Schwann myel MPZ+")
obj_can$cellType <- unname(can_it7_map[as.character(obj_can$RNA_snn_res.0.2)])
Idents(obj_can) <- obj_can$cellType
saveRDS(obj_can, file.path(DIR_SKETCH, 'can_it7.rds'), compress = FALSE)
saveRDS(obj_can, file.path(DIR_SKETCH, 'can_final.rds'), compress = FALSE)

# --- Lymphoid Sub-Clustering ---
cal_it1_map <- c("0"="Lymphoid", "1"="Epithelial", "2"="Myeloid", "3"="Fibroblast", "4"="Lymphoid", "5"="Lymphoid", "6"="Endothelial", "7"="Cycling", "8"="Myeloid", "9"="Low Quality", "10"="Nerve", "11"="Endothelial", "12"="Epithelial", "13"="Low Quality")
obj_cal$cellType <- unname(cal_it1_map[as.character(obj_cal$harmony.cluster.full)])
Idents(obj_cal) <- obj_cal$cellType
saveRDS(obj_cal, file.path(DIR_SKETCH, 'cal_it1.rds'), compress = FALSE)

obj_cal <- subset(obj_cal, cellType == 'Lymphoid')
obj_cal <- analyzeObj(obj_cal, 'cal_it2')

cal_it2_map <- c("0"="Lymphoid", "1"="Epithelial", "2"="Lymphoid", "3"="Lymphoid", "4"="Lymphoid", "5"="Myeloid", "6"="Low Quality", "7"="Low Quality")
obj_cal$cellType <- unname(cal_it2_map[as.character(obj_cal$harmony.cluster.full)])
Idents(obj_cal) <- obj_cal$cellType
saveRDS(obj_cal, file.path(DIR_SKETCH, 'cal_it2.rds'), compress = FALSE)

obj_cal <- subset(obj_cal, cellType == 'Lymphoid')
obj_cal <- analyzeObj(obj_cal, 'cal_it3')

cal_it3_map <- c("0"="Low Quality", "1"="Low Quality", "2"="Lymphoid", "3"="Lymphoid", "4"="Low Quality", "5"="Lymphoid", "6"="Lymphoid", "7"="Lymphoid", "8"="Lymphoid", "9"="Epithelial", "10"="Lymphoid", "11"="Lymphoid", "12"="Low Quality", "13"="Lymphoid", "14"="Lymphoid", "15"="Lymphoid", "16"="Lymphoid", "17"="Lymphoid", "18"="Lymphoid", "19"="Lymphoid", "20"="Low Quality", "21"="Lymphoid", "22"="Low Quality", "23"="Lymphoid", "24"="Lymphoid", "25"="Low Quality", "26"="Low Quality", "27"="Lymphoid", "28"="Low Quality", "29"="Low Quality")
obj_cal$cellType <- unname(cal_it3_map[as.character(obj_cal$harmony.cluster.full)])
Idents(obj_cal) <- obj_cal$cellType
saveRDS(obj_cal, file.path(DIR_SKETCH, 'cal_it3.rds'), compress = FALSE)

obj_cal <- subset(obj_cal, cellType == 'Lymphoid')
obj_cal <- analyzeObj(obj_cal, 'cal_it4')

cal_it4_map <- c("0"="Lymphoid", "1"="Lymphoid", "2"="Lymphoid", "3"="Low Quality", "4"="Lymphoid", "5"="Lymphoid", "6"="Lymphoid", "7"="Lymphoid", "8"="Lymphoid", "9"="Lymphoid", "10"="Lymphoid", "11"="Lymphoid", "12"="Low Quality", "13"="Lymphoid", "14"="Lymphoid", "15"="Lymphoid", "16"="Cycling", "17"="Lymphoid", "18"="Lymphoid", "19"="Low Quality", "20"="Low Quality", "21"="Low Quality")
obj_cal$cellType <- unname(cal_it4_map[as.character(obj_cal$harmony.cluster.full)])
Idents(obj_cal) <- obj_cal$cellType
saveRDS(obj_cal, file.path(DIR_SKETCH, 'cal_it4.rds'), compress = FALSE)

obj_cal <- subset(obj_cal, cellType == 'Lymphoid')
obj_cal <- analyzeObj(obj_cal, 'cal_it5')

cal_it5_map <- c("0"="T CD8+/GZMK+", "1"="Treg FOXP3+", "2"="T CD8+/GZMK+", "3"="NK KLRD1+", "4"="Plasma JCHAIN+", "5"="T CD4+/IL7R+", "6"="B mem BANK1", "7"="B naive IGHM+", "8"="NK KLRD1+", "9"="T CD8+/GZMB+", "10"="T CD4+/IL7R+", "11"="Plasma JCHAIN+", "12"="T CD8+/GZMB+", "13"="Plasma JCHAIN+", "14"="T CD8+/GZMK+", "15"="Tfh CXCL13+", "16"="Plasma JCHAIN+", "17"="T CD8+/ISG15+", "18"="B germinal RGS13+", "19"="Low Quality", "20"="Low Quality", "21"="Low Quality", "22"="Low Quality", "23"="Low Quality")
obj_cal$cellType <- unname(cal_it5_map[as.character(obj_cal$harmony.cluster.full)])
Idents(obj_cal) <- obj_cal$cellType
saveRDS(obj_cal, file.path(DIR_SKETCH, 'cal_it5.rds'), compress = FALSE)

obj_cal <- subset(obj_cal, cellType != 'Low Quality')
saveRDS(obj_cal, file.path(DIR_SKETCH, 'cal_final.rds'), compress = FALSE)

# --- Myeloid Sub-Clustering ---
cam_it1_map <- c("0"="Myeloid", "1"="Low Quality", "2"="Epithelial", "3"="Lymphoid", "4"="Myeloid", "5"="Fibroblast", "6"="Myeloid", "7"="Cycling", "8"="Myeloid", "9"="Endothelial", "10"="Lymphoid", "11"="Lympohid", "12"="Low Quality")
obj_cam$cellType <- unname(cam_it1_map[as.character(obj_cam$harmony.cluster.full)])
Idents(obj_cam) <- obj_cam$cellType
saveRDS(obj_cam, file.path(DIR_SKETCH, 'cam_it1.rds'), compress = FALSE)

obj_cam <- subset(obj_cam, cellType == 'Myeloid')
obj_cam <- analyzeObj(obj_cam, 'cam_it2')

cam_it2_map <- c("0"="Myeloid", "1"="Myeloid", "2"="Myeloid", "3"="Myeloid", "4"="Myeloid", "5"="Low Quality", "6"="Myeloid", "7"="Myeloid", "8"="Lymphoid", "9"="Myeloid", "10"="Myeloid", "11"="Myeloid", "12"="Myeloid", "13"="Epithelial", "14"="Myeloid", "15"="Myeloid", "16"="Lymphoid", "17"="Myeloid", "18"="Myeloid", "19"="Myeloid", "20"="Myeloid", "21"="Myeloid", "22"="Fibroblast", "23"="Myeloid", "24"="Myeloid", "25"="Low Quality", "26"="Low Quality", "27"="Low Quality")
obj_cam$cellType <- unname(cam_it2_map[as.character(obj_cam$harmony.cluster.full)])
Idents(obj_cam) <- obj_cam$cellType
saveRDS(obj_cam, file.path(DIR_SKETCH, 'cam_it2.rds'), compress = FALSE)

obj_cam <- subset(obj_cam, cellType == 'Myeloid')
obj_cam <- analyzeObj(obj_cam, 'cam_it3')

cam_it3_map <- c("0"="Neut CSF3R+", "1"="Mac C1QC+", "2"="Mast CPA3+", "3"="Mono CD14+", "4"="DC cDC2 CD1C+", "5"="Mac MARCO+", "6"="Low Quality", "7"="Mac ISG15+", "8"="Mac SPP1+", "9"="Mac C1QC+", "10"="Low Quality", "11"="Mac C1QC+", "12"="Mono CD16+", "13"="DC mregDC LAMP3+", "14"="Mac C1QC+", "15"="Mac ISG15+", "16"="Mac C1QC+", "17"="DC cDC1 CLEC9A+", "18"="Osteoclast CTSK+", "19"="Mac C1QC+", "20"="Low Quality", "21"="Mac SPP1+", "22"="Mac C1QC+", "23"="Mac C1QC+", "24"="Low Quality", "25"="Low Quality", "26"="Low Quality", "27"="Low Quality")
obj_cam$cellType <- unname(cam_it3_map[as.character(obj_cam$harmony.cluster.full)])
Idents(obj_cam) <- obj_cam$cellType
saveRDS(obj_cam, file.path(DIR_SKETCH, 'cam_it3.rds'), compress = FALSE)

obj_cam <- subset(obj_cam, cellType != 'Low Quality')
saveRDS(obj_cam, file.path(DIR_SKETCH, 'cam_final.rds'), compress = FALSE)


# ==============================================================================
# 8. Final Plotting Workflows
# ==============================================================================

# Ensure latest joined structures exist (Assumes these are available via upstream execution)
if (file.exists(file.path(DIR_SKETCH, 'caf_final_joined.rds'))) {
  obj_caf <- readRDS(file.path(DIR_SKETCH, 'caf_final_joined.rds'))
  makeFinalFigures(obj_caf, 'caf_final', genes = c('CD74','CXCL8','ISG15','COL4A1','LRRC15','CXCL14','PI16'))
  genes <- c("PDGFRA", "LUM", "DCN", "CD74","HLA-DRA", "HLA-DRB1", "CXCL8", "CXCL1", "CXCL3", "ISG15", "CXCL10", "IFIT3", "COL4A1","COL4A2","ITGA1", "LRRC15","SDC1","CTHRC1", "COL14A1", "IGF1", "NFIB", "PI16", "CFD", "CD34")
  neg <- c("PECAM1", "RGS5", "EPCAM", "PLP1","MSLN")
  makeSuppPanel(obj_caf, 'caf', genes, neg, dotColor = 'black')
}

if (file.exists(file.path(DIR_SKETCH, 'cec_final_joined.rds'))) {
  obj_cec <- readRDS(file.path(DIR_SKETCH, 'cec_final_joined.rds'))
  makeFinalFigures(obj_cec, 'cec_final', genes = c('GJA5','CA4','PROX1','ESM1','ACKR1'))
  genes <- c("PECAM1", "VWF", "CDH5", "GJA5","HEY1", "SEMA3G", "CA4", "CD36","MT1E", "PROX1","PDPN","CCL21", "ESM1","KDR","FLT1", "ACKR1", "SELE", "CLU", "PLVAP", "IGFBP3", "HSPG2", "COL4A1", "SELP", "POSTN")
  neg <- c("EPCAM", "PDGFRA", "RGS5", "PLP1","PTPRC")
  makeSuppPanel(obj_cec, 'cec', genes, neg, dotColor = 'black')
}

if (file.exists(file.path(DIR_SKETCH, 'cmc_final_joined.rds'))) {
  obj_cmc <- readRDS(file.path(DIR_SKETCH, 'cmc_final_joined.rds'))
  makeFinalFigures(obj_cmc, 'cmc_final', genes = c('CCL2','CD248','ISG15','RERGL','WFDC1','DES'))
  genes <- c("RGS5", "PDGFRB", "KCNJ8", "CSPG4","PROCR", "MYH11", "CCL2", "CCL19", "CCL21", "CD248","THY1","PLXDC1", "ISG15","CXCL10","IFIT3", "RERGL", "NET1", "MT1A", "WFDC1", "LTBP1", "SLIT3", "DES", "ACTG2", "PRUNE2")
  neg <- c("EPCAM", "PDGFRA", "PECAM1", "PLP1","PTPRC")
  makeSuppPanel(obj_cmc, 'cmc', genes, neg, dotColor = 'black')
}

if (file.exists(file.path(DIR_SKETCH, 'can_final_joined.rds'))) {
  obj_can <- readRDS(file.path(DIR_SKETCH, 'can_final_joined.rds'))
  makeFinalFigures(obj_can, 'can_final', genes = c('MLANA','MPZ','NGFR'))
  genes <- c("MLANA","TYRP1", "PMEL", "MPZ", "NRXN1", "NTM", "NGFR","TNC","GBP1")
  neg <- c("EPCAM", "PDGFRA", "PECAM1", "RGS5","PTPRC")
  makeSuppPanel(obj_can, 'can', genes, neg, dotColor = 'black', rasterFeaturePlot = FALSE)
}

if (file.exists(file.path(DIR_SKETCH, 'cal_final_joined.rds'))) {
  obj_cal <- readRDS(file.path(DIR_SKETCH, 'cal_final_joined.rds'))
  makeFinalFigures(obj_cal, 'cal_final', genes = c('RGS13','BANK1','IGHM','KLRD1','JCHAIN','IL7R','GZMB','GZMK','ISG15','CXCL13','FOXP3'))
  genes <- c("MS4A1", "CD79A","HLA-DRA", "RGS13","BANK1", "IGHM", "KLRD1","GNLY","NKG7", "JCHAIN", "MZB1", "DERL3", "CD4","CD2","CCR7", "IL7R", "CXCL13", "FOXP3", "CD8A", "CD8B", "CCL5", "GZMB", "GZMK", "ISG15")
  neg <- c("EPCAM", "PDGFRA", "PECAM1", "RGS5","PLP1")
  makeSuppPanel(obj_cal, 'cal', genes, neg, dotColor = 'black')
}

if (file.exists(file.path(DIR_SKETCH, 'cam_final.rds'))) {
  obj_cam <- readRDS(file.path(DIR_SKETCH, 'cam_final.rds'))
  makeFinalFigures(obj_cam, 'cam_final', genes = c('CLEC9A','CD1C','LAMP3','C1QC','ISG15','MARCO','SPP1','CPA3','VCAN','CD52','CSF3R','CTSK'))
  genes <- c("HLA-DRA", "FLT3", "ITGAX", "CLEC9A","CD1C", "LAMP3", "CD68", "ITGAM", "C1QC", "ISG15","MARCO","SPP1", "FCN1","VCAN","CD52", "CPA3", "CSF3R", "CTSK")
  neg <- c("EPCAM", "PDGFRA", "PECAM1", "RGS5","PLP1")
  makeSuppPanel(obj_cam, 'cam', genes, neg, dotColor = 'black')
}
