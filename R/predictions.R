#' Run scTME Predictions
#' 
#' @importFrom SingleR SingleR
#' @importFrom Seurat as.SingleCellExperiment NormalizeData GetAssayData Idents
#' @export
#' @param obj The query dataset
#' @param ref The reference atlas
#' @param clusters Optional cluster labels
#' @param level Resolution level (1 or 2)
#' @param label.col Target metadata column
#' @param max.cells Downsampling threshold
#' @param fast Downsample the reference atlas
#' @param return.details Return prediction strength metrics
scTME <- function (obj, ref = NULL, clusters = NULL, level = 2, label.col = NULL, max.cells = 50000, fast = FALSE, return.details = FALSE) {

    # 0. Track if we are using the default atlas
    using_default_ref <- is.null(ref)
    if (using_default_ref) {
        message("No reference provided. Loading built-in scTME atlas...")
        utils::data("scTME_atlas", package = "scTME", envir = environment())
        ref <- scTME_atlas
    }

    # Helper function: safely check if a Seurat object needs normalization (v4 & v5 compatible)
    needs_norm <- function(seurat_obj) {
        tryCatch({
            suppressWarnings({
                counts_mat <- Seurat::GetAssayData(seurat_obj, layer = "counts")
                data_mat <- Seurat::GetAssayData(seurat_obj, layer = "data")
                identical(counts_mat, data_mat)
            })
        }, error = function(e) {
            tryCatch({
                suppressWarnings({
                    counts_mat <- Seurat::GetAssayData(seurat_obj, slot = "counts")
                    data_mat <- Seurat::GetAssayData(seurat_obj, slot = "data")
                    identical(counts_mat, data_mat)
                })
            }, error = function(e2) {
                TRUE 
            })
        })
    }

    # 1. Defensive Prep: Normalize and clean the Query (obj)
    if (inherits(obj, "Seurat")) {
        if (needs_norm(obj)) {
            message("Log-normalizing query object...")
            obj <- Seurat::NormalizeData(obj, verbose = FALSE)
        }
        obj@reductions <- list()
    }

    # 2. Defensive Prep: Normalize the Reference (ref)
    if (inherits(ref, "Seurat")) {
        if (needs_norm(ref)) {
            message("Log-normalizing reference object...")
            ref <- Seurat::NormalizeData(ref, verbose = FALSE)
        }
        ref@reductions <- list()
    }

    # 3. Downsample the query object FIRST
    if (ncol(obj) > max.cells) {
        keep <- sample(colnames(obj), max.cells)
        obj <- obj[, keep]
    }
    
    # 3.5 Downsample the reference for fast mode
    if (fast && ncol(ref) > 5000) {
        message("Fast mode enabled: downsampling reference to 5,000 cells...")
        keep_ref <- sample(colnames(ref), 5000)
        ref <- ref[, keep_ref]
    }

    # 4. Convert query to SingleCellExperiment
    if (!inherits(obj, "SingleCellExperiment")) {
        sce = Seurat::as.SingleCellExperiment(obj)
    } else {
        sce = obj
    }
    
    # 5. Convert reference to SingleCellExperiment
    if (!inherits(ref, "SingleCellExperiment")) {
        ref_sce = Seurat::as.SingleCellExperiment(ref)
    } else {
        ref_sce = ref
    }

    # 6. Extract query clusters dynamically
    cluster_vec <- NULL
    if (!is.null(clusters)) {
        if (length(clusters) == 1 && is.character(clusters)) {
            if (clusters == "auto") {
                cluster_vec <- if (inherits(obj, "Seurat")) Seurat::Idents(obj) else SingleCellExperiment::colLabels(sce)
            } else {
                if (inherits(obj, "Seurat")) {
                    if (!clusters %in% colnames(obj@meta.data)) stop(paste("Cluster column", clusters, "not found in query metadata."))
                    cluster_vec <- obj@meta.data[[clusters]]
                } else {
                    if (!clusters %in% colnames(sce@colData)) stop(paste("Cluster column", clusters, "not found in query metadata."))
                    cluster_vec <- sce@colData[[clusters]]
                }
            }
        } else {
            cluster_vec <- clusters 
        }
    }

    # 7. Extract reference labels dynamically
    if (using_default_ref && is.null(label.col)) {
        if (level == 1) {
            target_col <- "cellType_fine_abbr"
        } else if (level == 2) {
            target_col <- "cellType_fine"
        } else {
            stop("Error: 'level' must be 1 (coarse) or 2 (fine).")
        }
        
        if (!target_col %in% colnames(ref_sce@colData)) {
            stop(paste0("Error: Target atlas column '", target_col, "' not found. Available metadata columns are: ", paste(colnames(ref_sce@colData), collapse = ", ")))
        }
        labels <- ref_sce@colData[[target_col]]
        
    } else if (is.null(label.col)) {
        if (inherits(ref, "Seurat")) {
            labels <- Seurat::Idents(ref)
        } else {
            if (!is.null(SingleCellExperiment::colLabels(ref_sce))) {
                labels <- SingleCellExperiment::colLabels(ref_sce)
            } else if ("ident" %in% colnames(ref_sce@colData)) {
                labels <- ref_sce@colData[["ident"]]
            } else {
                stop("Error: No default labels found in reference. Please specify label.col.")
            }
        }
    } else {
        if (!label.col %in% colnames(ref_sce@colData)) {
            stop(paste0("Error: Column '", label.col, "' not found in reference metadata."))
        }
        labels <- ref_sce@colData[[label.col]]
    }

    # 8. Run the SingleR Pipeline
    common <- intersect(rownames(sce), rownames(ref_sce))
    sce <- sce[common, ]
    ref_sce <- ref_sce[common, ]
    
    if( !is.null(cluster_vec) ) {
      pred <- SingleR(test = sce, ref = ref_sce, labels = labels, de.method = "wilcox", clusters = cluster_vec)
      newLabels <- pred$labels
      names(newLabels) <- rownames(pred)
      tempLabels <- newLabels[as.character(cluster_vec)]
    } else {
      pred <- SingleR(test = sce, ref = ref_sce, labels = labels, de.method = "wilcox" )
      tempLabels <- pred$labels
    }
    
    # 9. Return Logic
    if (return.details) {
        returnList = list(tempLabels, pred)
        names(returnList) = c("cell.level.predictions", "cluster.level.details")
        return(returnList)
    } else {
        return(tempLabels)
    }
}

