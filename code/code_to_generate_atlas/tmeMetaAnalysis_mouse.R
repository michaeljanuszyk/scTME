# ==============================================================================
# Script: tmeMetaAnalysis_mouse.R
# Description: Full pipeline for integrating and analyzing mouse scTME datasets.
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
source('tmeFunctions.R')

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
DIR_REF       <- file.path("..","reference")

# Define full paths for the specific CSVs
FILE_SYNONYMS <- "GRCm39_geneNames_geneSynonyms.csv"
FILE_METADATA <- "mouse.meta.data.csv"

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

# QC Execution (Mouse uses lowercase "^mt-" prefix)
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
  
  obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^mt-")
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

obj_list <- as.vector(read.csv(file.path(DIR_QC, "obj_list_01_preQC.csv"), colClasses = c("NULL", NA))$x)
pb <- txtProgressBar(min = 1, max = length(obj_list), style = 3)
count <- 0

# Format individual post-QC objects (streamline mouse metadata, rename gene names)
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
                    "POD", "size", "age", "sequencingmethod", "treated", 
                    "model", "region", "sorting", "site", "percent.mt", "type", "strain")))
  x <- RenameCells(x, new.names = paste(x$orig.ident, Cells(x), sep = "_"))
  
  # Deal with duplicate gene names and non-standard gene names using GRCm39 synonyms
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
source_python("pickle_reader.py")
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

# Standardize rows into common mouse gene space
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
saveRDS(merged, file.path(DIR_SKETCH, "mouseMerged.rds"), compress = FALSE)

# Attach clinical/experimental metadata
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
obj$seurat_clusters <- obj$sketch_snn_res.0.2
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
for (i in c(0.05, 0.06, 0.07, 0.09, 1.0, 1.5, 2.0)) {
  obj_mnn <- FindClusters(obj_mnn, resolution = i, graph.name = 'RNA_snn')
  DimPlot(obj_mnn, raster = TRUE, label = TRUE, reduction = 'umap.mnn')
  ggsave(file.path(DIR_SKETCH, paste0("unfiltered_counts_obj_findvar_sketch1000_mnn_umap_res", i, ".jpg")), width = 5, height = 5, units = "in", limitsize = FALSE)
  saveRDS(obj_mnn, file.path(DIR_SKETCH, paste0("unfiltered_counts_obj_findvar_sketch1000_mnn_umap_res", i, ".rds")), compress = FALSE)
}

obj_mnn$seurat_clusters <- obj_mnn$RNA_snn_res.0.15
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
for (i in c(0.05, 0.06, 0.07, 0.09, 1.0, 1.5, 2.0)) {
  obj_scvi <- FindClusters(obj_scvi, resolution = i, graph.name = 'RNA_snn')
  DimPlot(obj_scvi, raster = TRUE, reduction = 'umap.scvi', label = TRUE)
  ggsave(file.path(DIR_SKETCH, paste0("unfiltered_counts_obj_findvar_sketch1000_scvi_umap_res", i, ".jpg")), width = 5, height = 5, units = "in", limitsize = FALSE)
  saveRDS(obj_scvi, file.path(DIR_SKETCH, paste0("unfiltered_counts_obj_findvar_sketch1000_scvi_umap_res", i, ".rds")), compress = FALSE)
}

obj_scvi$seurat_clusters <- obj_scvi$RNA_snn_res.0.4
obj_scvi <- ProjectIntegration(object = obj_scvi, reduction = "integrated.scvi")
obj_scvi <- ProjectData(object = obj_scvi, sketched.reduction = "integrated.scvi.full", full.reduction = "integrated.scvi.full", umap.model = "umap.scvi", dims = 1:30, refdata = list(scvi.cluster.full = "seurat_clusters"))
saveRDS(obj_scvi, file.path(DIR_SKETCH, "unfiltered_counts_obj_findvar_sketch1000_scvi_umap_res0.4_project.rds"), compress = FALSE)


# ==============================================================================
# 6. Cluster Annotations (Dictionary Mapping)
# ==============================================================================

# --- Harmony Annotations ---
harmony_map <- c(
  "0" = "Epithelial", "1" = "Fibroblast", "2" = "Myeloid", "3" = "Lymphoid",
  "4" = "Endothelial", "5" = "Cycling", "6" = "Myeloid", "7" = "Lymphoid",
  "8" = "Fibroblast", "9" = "Mural", "10" = "Fibroblast", "11" = "Epithelial",
  "12" = "Epithelial", "13" = "Nerve", "14" = "Epithelial", "15" = "Lymphoid", "16" = "Lymphoid"
)
obj$cellType <- unname(harmony_map[as.character(obj$harmony.cluster.full)])

# --- MNN Annotations ---
mnn_map <- c(
  "0" = "Epithelial", "1" = "Fibroblast", "2" = "Myeloid", "3" = "Lymphoid",
  "4" = "Endothelial", "5" = "Myeloid", "6" = "Lymphoid", "7" = "Epithelial",
  "8" = "Epithelial", "9" = "Mural", "10" = "Fibroblast", "11" = "Epithelial",
  "12" = "Epithelial", "13" = "Epithelial", "14" = "Epithelial"
)
obj_mnn$cellType <- unname(mnn_map[as.character(obj_mnn$mnn.cluster.full)])
obj_mnn$cellType[obj_mnn$RNA_snn_res.4 == 42] <- 'Nerve'

# --- scVI Annotations ---
scvi_map <- c(
  "0" = "Fibroblast", "1" = "Myeloid", "2" = "Lymphoid", "3" = "Fibroblast",
  "4" = "Endothelial", "5" = "Myeloid", "6" = "Epithelial", "7" = "Epithelial",
  "8" = "Epithelial", "9" = "Epithelial", "10" = "Epithelial", "11" = "Lymphoid",
  "12" = "Myeloid", "13" = "Mural", "14" = "Nerve", "15" = "Fibroblast",
  "16" = "Epithelial", "17" = "Myeloid", "18" = "Epithelial", "19" = "Lymphoid",
  "20" = "Epithelial", "21" = "Endothelial", "22" = "Epithelial", "23" = "Nerve",
  "24" = "Lymphoid", "25" = "Epithelial", "26" = "Mural", "27" = "Myeloid",
  "28" = "Epithelial", "29" = "Lymphoid", "30" = "Myeloid", "31" = "Fibroblast",
  "32" = "Lymphoid", "33" = "Myeloid", "34" = "Epithelial", "35" = "Epithelial"
)
obj_scvi$cellType <- unname(scvi_map[as.character(obj_scvi$scvi.cluster.full)])


# ==============================================================================
# 7. Sub-Clustering and Fine Annotations
# ==============================================================================

obj <- combineIntegration('Fibroblast')
obj <- analyzeObj(obj, 'caf_it1')
obj <- combineIntegration('Endothelial')
obj <- analyzeObj(obj, 'cec_it1')
obj <- combineIntegration('Mural')
obj <- analyzeObj(obj, 'cmc_it1')
obj <- combineIntegration('Nerve')
obj <- analyzeObj(obj, 'can_it1')
obj <- combineIntegration('Lymphoid')
obj <- analyzeObj(obj, 'cal_it1')
obj <- combineIntegration('Myeloid')
obj <- analyzeObj(obj, 'cam_it1')

# --- Fibroblast Sub-Clustering Iterations ---
caf_it1_map <- c("0"="Fibroblast", "1"="Cycling", "2"="Epithelial", "3"="Myeloid", "4"="Fibroblast", "5"="Endothelial", "6"="Nerve", "7"="Fibroblast", "8"="Lymphoid", "9"="Low Quality", "10"="Fibroblast", "11"="Epithelial", "12"="Epithelial", "13"="Myeloid", "14"="Epithelial", "15"="Low Quality", "16"="Fibroblast", "17"="Low Quality")
obj$cellType <- unname(caf_it1_map[as.character(obj$RNA_snn_res.0.09)])
saveRDS(obj, file.path(DIR_SKETCH, 'caf_it1.rds'), compress = FALSE)
obj <- subset(obj, cellType == 'Fibroblast')
obj <- analyzeObj(obj, 'caf_it2')

caf_it2_map <- c("0"="Fibroblast", "1"="Fibroblast", "2"="Fibroblast", "3"="Cycling", "4"="Lymphoid", "5"="Lymphoid")
obj$cellType <- unname(caf_it2_map[as.character(obj$RNA_snn_res.0.06)])
saveRDS(obj, file.path(DIR_SKETCH, 'caf_it2.rds'), compress = FALSE)
obj <- subset(obj, cellType == 'Fibroblast')
obj <- analyzeObj(obj, 'caf_it3')

caf_it3_map <- c("0"="Fibroblast", "1"="Fibroblast", "2"="Fibroblast", "3"="Fibroblast", "4"="Fibroblast", "5"="Fibroblast", "6"="Epithelial", "7"="Low Quality", "8"="Epithelial", "9"="Fibroblast", "10"="Fibroblast", "11"="Fibroblast", "12"="Fibroblast", "13"="Low Quality")
obj$cellType <- unname(caf_it3_map[as.character(obj$RNA_snn_res.0.18)])
saveRDS(obj, file.path(DIR_SKETCH, 'caf_it3.rds'), compress = FALSE)
obj <- subset(obj, cellType == 'Fibroblast')
obj <- analyzeObj(obj, 'caf_it4')

caf_it4_map <- c("0"="mCAF Lrrc15+", "1"="ssCAF Cxcl12+", "2"="Low Quality", "3"="ssCAF Pi16+", "4"="apCAF Cd74+", "5"="iCAF Isg15+", "6"="ssCAF Pi16+", "7"="Pericyte", "8"="Low Quality", "9"="apCAF Cd74+", "10"="ssCAF Pi16+", "11"="mCAF Lrrc15+", "12"="ssCAF Cxcl12+", "13"="Epithelial")
obj$cellType <- unname(caf_it4_map[as.character(obj$RNA_snn_res.0.2)])
saveRDS(obj, file.path(DIR_SKETCH, 'caf_it4.rds'), compress = FALSE)
obj <- subset(obj, cellType != 'Low Quality' & cellType != 'Pericyte' & cellType != 'Epithelial')
obj <- analyzeObj(obj, 'caf_it5')

caf_it5_map <- c("0"="ssCAF Cxcl12+", "1"="mCAF Lrrc15+", "2"="ssCAF Pi16+", "3"="iCAF Spp1+", "4"="apCAF Cd74+", "5"="iCAF Isg15+", "6"="ssCAF Cxcl12+", "7"="apCAF Cd74+", "8"="apCAF Cd74+", "9"="ssCAF Pi16+", "10"="mCAF Lrrc15+", "11"="apCAF Cd74+")
obj$cellType <- unname(caf_it5_map[as.character(obj$RNA_snn_res.0.16)])
Idents(obj) <- obj$cellType
saveRDS(obj, file.path(DIR_SKETCH, 'caf_final.rds'), compress = FALSE)

# --- Endothelial Sub-Clustering ---
cec_it1_map <- c("0"="Endothelial", "1"="Low Quality", "2"="Fibroblast", "3"="Endothelial", "4"="Myeloid", "5"="Lymphoid", "6"="Epithelial", "7"="Endothelial", "8"="Lymphoid", "9"="Cycling", "10"="Myeloid")
obj$cellType <- unname(cec_it1_map[as.character(obj$RNA_snn_res.0.1)])
saveRDS(obj, file.path(DIR_SKETCH, 'cec_it1.rds'), compress = FALSE)
obj <- subset(obj, cellType == 'Endothelial')
obj <- analyzeObj(obj, 'cec_it2')

cec_it2_map <- c("0"="Endothelial", "1"="Endothelial", "2"="Endothelial", "3"="Endothelial", "4"="Lymphoid", "5"="Low Quality", "6"="Lymphoid")
obj$cellType <- unname(cec_it2_map[as.character(obj$RNA_snn_res.0.1)])
saveRDS(obj, file.path(DIR_SKETCH, 'cec_it2.rds'), compress = FALSE)
obj <- subset(obj, cellType == 'Endothelial')
obj <- analyzeObj(obj, 'cec_it3')

cec_it3_map <- c("0"="Capillary Car4+", "1"="Vein Ackr1+", "2"="Tip Esm1+", "3"="Lymphatic Prox1+", "4"="Artery Gja5+", "5"="Capillary Car4+")
obj$cellType <- unname(cec_it3_map[as.character(obj$RNA_snn_res.0.14)])
Idents(obj) <- obj$cellType
saveRDS(obj, file.path(DIR_SKETCH, 'cec_final.rds'), compress = FALSE)

# --- Mural Sub-Clustering ---
cmc_it1_map <- c("0"="Mural", "1"="Fibroblast", "2"="Low Quality", "3"="Fibroblast", "4"="Myeloid", "5"="Cycling", "6"="Endothelial", "7"="Fibroblast", "8"="Low Quality")
obj$cellType <- unname(cmc_it1_map[as.character(obj$RNA_snn_res.0.06)])
saveRDS(obj, file.path(DIR_SKETCH, 'cmc_it1.rds'), compress = FALSE)
obj <- subset(obj, cellType == 'Mural')
obj <- analyzeObj(obj, 'cmc_it2')

cmc_it2_map <- c("0"="Smooth muscle", "1"="Pericyte Cd248+", "2"="Pericyte Ccl2+", "3"="Pericyte Cd248+", "4"="Smooth muscle")
obj$cellType <- unname(cmc_it2_map[as.character(obj$RNA_snn_res.0.16)])
saveRDS(obj, file.path(DIR_SKETCH, 'cmc_final.rds'), compress = FALSE)


# ==============================================================================
# 8. Final Plotting Workflows
# ==============================================================================

if (file.exists(file.path(DIR_SKETCH, 'caf_final.rds'))) {
  obj_caf <- readRDS(file.path(DIR_SKETCH, 'caf_final.rds'))
  makeFinalFigures(obj_caf, 'mouse_Fibroblast', c('Cd74','Isg15','Spp1','Lrrc15','Cxcl12','Pi16'))
}

if (file.exists(file.path(DIR_SKETCH, 'cec_final.rds'))) {
  obj_cec <- readRDS(file.path(DIR_SKETCH, 'cec_final.rds'))
  makeFinalFigures(obj_cec, 'mouse_Endothelial', c('Gja5','Car4','Prox1','Esm1','Ackr1'))
}

if (file.exists(file.path(DIR_SKETCH, 'cmc_final.rds'))) {
  obj_cmc <- readRDS(file.path(DIR_SKETCH, 'cmc_final.rds'))
  makeFinalFigures(obj_cmc, 'mouse_Mural', c('Ccl2','Cd248'))
}
