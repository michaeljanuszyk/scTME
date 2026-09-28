###############################################################################
## subpopulationEnrichment.R
##
## Enrichment of cell subpopulations between groups of samples defined by
## sample-level meta-data, for an integrated single-cell atlas.
##
## Usage
## -----
## Assign the integrated Seurat object, then source this file:
##
##   obj <- readRDS("atlas.rds")
##   source("subpopulationEnrichment.R")
##
## ATLAS may be set instead of `obj`, either to a Seurat object or to the path
## of an .rds file holding one. OUTDIR and OUTFILE set the destination of the
## results table.
##
## Which comparisons run depends on which meta-data columns the object carries:
## tumor vs non-tumor tissue requires `site`; advanced vs early stage requires
## `clinical_stage`; treatment non-responders vs responders requires
## `treatment_response`. Comparisons whose columns are absent are skipped, so
## the same file serves atlases with different depths of clinical annotation.
##
## Method summary
## --------------
## For each parent cell lineage (fibroblast, endothelial, mural, lymphoid,
## myeloid) and each constituent subpopulation, the percentage of cells
## belonging to that subpopulation is computed per sample. The unit of
## replication is therefore the individual scRNA-seq sample, never the
## individual cell.
##
## Samples are weighted by the square root of their cell number within the
## lineage, which limits the influence of samples with very high capture rates
## without discarding information from smaller samples. Weighted means and
## standard errors use the Kish effective sample size,
## n_eff = (sum w)^2 / sum(w^2).
##
## Fold-enrichment is the ratio of the two groups' weighted mean proportions;
## its standard error is obtained by propagating the two weighted standard
## errors, and its 95% confidence interval by the delta method on the log
## scale (so the interval cannot extend below zero).
##
## Significance is assessed by a two-sided weighted Welch t-test for two
## independent samples on the per-sample proportions. Because subpopulations
## are nested within parent lineages, Benjamini-Hochberg adjustment is applied
## separately within each lineage; a global adjustment across all
## subpopulations is also reported for completeness.
##
## Dependencies: dplyr, tidyr, weights
###############################################################################

library(dplyr)
library(tidyr)
library(weights)

###############################################################################
## Analysis parameters
###############################################################################

## The five stromal and immune lineages that make up the shared tumor
## microenvironment. Only these are analyzed; epithelial, neural and
## unassigned cells are outside the scope of the comparison.
LINEAGES <- c("Fib", "End", "Mur", "Lym", "Myl")


###############################################################################
## Functions
###############################################################################

#' Per-cell metadata from a Seurat object
#'
#' Accepts a Seurat object, the path of an .rds file holding one, or a
#' metadata data.frame.
#'
#' @param x    Seurat object, file path, or data.frame
#' @param what label used in error messages
#' @return a per-cell metadata data.frame
seurat_meta <- function(x, what = "input") {
  if (is.character(x)) {
    if (!file.exists(x)) stop(what, " not found: ", x)
    x <- readRDS(x)
  }
  if (inherits(x, "Seurat")) return(x@meta.data)
  if (is.data.frame(x))      return(x)
  stop(what, " must be a Seurat object, a path to one, or a data.frame.")
}


#' Locate the lineage and subpopulation annotation columns
#'
#' Annotation column names differ between atlas releases, so detect whichever
#' pair is present rather than hardcoding one. If only the subpopulation
#' column exists, the parent lineage is derived from its label prefix by
#' `calculate_enrichment()`.
annotation_cols <- function(meta) {
  pairs <- list(
    c(lineage = "integrated.annotation.abbr", cell_type = "integrated.annotation"),
    c(lineage = "cellType_fine_abbr",         cell_type = "cellType_fine")
  )
  for (p in pairs) if (all(p %in% names(meta))) return(p)

  for (ct in c("integrated.annotation", "cellType_fine")) {
    if (ct %in% names(meta)) return(c(lineage = NA_character_, cell_type = ct))
  }
  stop("Could not find a subpopulation annotation column.\n",
       "Available columns: ", paste(names(meta), collapse = ", "))
}


#' Weighted mean and standard error using the Kish effective sample size
#'
#' @param x numeric vector of per-sample subpopulation percentages
#' @param w numeric vector of sample weights (sqrt of cells per sample)
#' @return named numeric vector: wavg, wse, wsd, n_eff
weighted_stats <- function(x, w) {
  if (length(x) < 2 || sum(w) == 0) {
    return(c(wavg = NA_real_, wse = NA_real_, wsd = NA_real_, n_eff = NA_real_))
  }
  wavg  <- sum(w * x) / sum(w)
  n_eff <- (sum(w))^2 / sum(w^2)
  if (n_eff <= 1) {
    return(c(wavg = wavg, wse = NA_real_, wsd = NA_real_, n_eff = n_eff))
  }
  ## Bias correction for the weighted variance, analogous to the n/(n-1)
  ## correction of the unweighted sample variance.
  correction <- n_eff / (n_eff - 1)
  wsd <- sqrt(correction * sum(w * (x - wavg)^2) / sum(w))
  c(wavg = wavg, wse = wsd / sqrt(n_eff), wsd = wsd, n_eff = n_eff)
}


#' Per-sample subpopulation percentages within one parent lineage
#'
#' Subpopulations absent from a sample are retained as zero rather than
#' missing, so every sample contributes a value for every subpopulation of
#' that lineage.
#'
#' @return tibble with columns sample, group, cell_type, proportion, weight
lineage_proportions <- function(meta, sample_col, group_col, cell_type_col) {

  if (nrow(meta) == 0) return(NULL)

  sample_totals <- meta %>%
    group_by(sample = .data[[sample_col]], group = .data[[group_col]]) %>%
    summarise(total_cells = n(), .groups = "drop")

  cell_counts <- meta %>%
    group_by(sample    = .data[[sample_col]],
             group     = .data[[group_col]],
             cell_type = .data[[cell_type_col]]) %>%
    summarise(cell_type_count = n(), .groups = "drop")

  cell_counts %>%
    right_join(expand_grid(sample_totals,
                           cell_type = unique(meta[[cell_type_col]])),
               by = c("sample", "group", "cell_type")) %>%
    mutate(cell_type_count = replace_na(cell_type_count, 0),
           proportion      = (cell_type_count / total_cells) * 100,
           weight          = sqrt(total_cells)) %>%
    filter(weight > 0)
}


#' Enrichment of every cell subpopulation between two sample groups
#'
#' @param meta          a Seurat object, or its per-cell metadata data.frame
#' @param group_col     metadata column defining the two groups
#' @param group1        reference group (denominator of the fold-enrichment)
#' @param group2        comparison group (numerator)
#' @param sample_col    metadata column identifying the biological sample
#' @param lineage_col   metadata column giving the parent lineage
#' @param cell_type_col metadata column giving the subpopulation
#'
#' @return One row per subpopulation, carrying the exact n per group, the
#'   weighted group means and standard errors, the fold-enrichment with its
#'   standard error and 95% confidence interval, the t statistic, the
#'   Welch-Satterthwaite degrees of freedom, and raw and adjusted P values.
calculate_enrichment <- function(meta,
                                 group_col, group1, group2,
                                 sample_col    = "orig.ident",
                                 lineage_col   = NULL,
                                 cell_type_col = "integrated.annotation") {

  if (inherits(meta, "Seurat")) meta <- meta@meta.data

  ## Some atlas releases ship only the subpopulation label. The parent lineage
  ## is encoded as its first token ("Fib iCAF ISG15+" -> "Fib"), so derive it
  ## when no separate lineage column is available.
  if (is.null(lineage_col) || is.na(lineage_col) ||
      !lineage_col %in% names(meta)) {
    if (!cell_type_col %in% names(meta)) {
      stop("Column not found: ", cell_type_col,
           "\nAvailable columns: ", paste(names(meta), collapse = ", "))
    }
    meta[[".lineage"]] <- sub(" .*$", "", as.character(meta[[cell_type_col]]))
    lineage_col <- ".lineage"
  }

  missing <- setdiff(c(sample_col, group_col, lineage_col, cell_type_col),
                     names(meta))
  if (length(missing)) {
    stop("Column(s) not found: ", paste(missing, collapse = ", "),
         "\nAvailable columns: ", paste(names(meta), collapse = ", "))
  }

  observed <- unique(as.character(meta[[lineage_col]]))

  meta <- meta %>%
    filter(.data[[lineage_col]] %in% LINEAGES,
           .data[[group_col]] %in% c(group1, group2))

  if (nrow(meta) == 0) {
    stop("No cells left after lineage filtering. Expected lineages: ",
         paste(LINEAGES, collapse = ", "),
         "\nFound in '", lineage_col, "': ",
         paste(sort(observed), collapse = ", "))
  }

  out <- list()

  for (lineage in unique(meta[[lineage_col]])) {

    props <- lineage_proportions(
      meta[meta[[lineage_col]] == lineage, , drop = FALSE],
      sample_col, group_col, cell_type_col)
    if (is.null(props)) next

    rows <- list()

    for (ct in unique(props$cell_type)) {

      g1 <- props[props$cell_type == ct & props$group == group1, ]
      g2 <- props[props$cell_type == ct & props$group == group2, ]

      s1 <- weighted_stats(g1$proportion, g1$weight)
      s2 <- weighted_stats(g2$proportion, g2$weight)

      ## Fold-enrichment and its propagated standard error. se_log is the
      ## standard error of log(ratio) and gives the delta-method interval.
      ratio  <- as.numeric(s2["wavg"] / s1["wavg"])
      se_log <- as.numeric(sqrt((s2["wse"] / s2["wavg"])^2 +
                                (s1["wse"] / s1["wavg"])^2))
      fse    <- abs(ratio) * se_log

      ## Two-sided weighted Welch t-test on the per-sample proportions.
      ## samedata = FALSE is required: the two groups are independent sets of
      ## samples of unequal size, not paired observations.
      t_value <- df <- p_raw <- NA_real_
      if (nrow(g1) >= 2 && nrow(g2) >= 2) {
        tt <- wtd.t.test(x = g1$proportion, weight  = g1$weight,
                         y = g2$proportion, weighty = g2$weight,
                         samedata = FALSE)
        t_value <- as.numeric(tt$coefficients["t.value"])
        df      <- as.numeric(tt$coefficients["df"])
        p_raw   <- as.numeric(tt$coefficients["p.value"])
      }

      rows[[ct]] <- data.frame(
        comparison       = paste0(group2, " vs ", group1),
        lineage          = lineage,
        cell_type        = ct,
        n_group1         = nrow(g1),
        n_group2         = nrow(g2),
        wmean_group1     = as.numeric(s1["wavg"]),
        wse_group1       = as.numeric(s1["wse"]),
        n_eff_group1     = as.numeric(s1["n_eff"]),
        wmean_group2     = as.numeric(s2["wavg"]),
        wse_group2       = as.numeric(s2["wse"]),
        n_eff_group2     = as.numeric(s2["n_eff"]),
        fold_enrichment  = ratio,
        fold_enrich_se   = fse,
        fold_ci_low      = exp(log(ratio) - 1.96 * se_log),
        fold_ci_high     = exp(log(ratio) + 1.96 * se_log),
        t_value          = t_value,
        df               = df,
        p_val_raw        = p_raw,
        stringsAsFactors = FALSE
      )
    }

    lineage_df <- bind_rows(rows) %>% filter(!is.na(p_val_raw))
    if (nrow(lineage_df) == 0) next

    ## Subpopulations partition their parent lineage, so each lineage is its
    ## own family of tests for multiple-comparison purposes.
    lineage_df$p_val_adj_lineage <- p.adjust(lineage_df$p_val_raw, method = "BH")
    out[[lineage]] <- lineage_df
  }

  res <- bind_rows(out)
  if (nrow(res) == 0) return(res)

  ## Reported for completeness alongside the per-lineage adjustment.
  res$p_val_adj_global <- p.adjust(res$p_val_raw, method = "BH")

  ## P values below double precision are reported at that floor. The raw
  ## values are retained above so nothing is lost.
  res$p_val_raw_reported         <- pmax(res$p_val_raw,         2.2e-16)
  res$p_val_adj_lineage_reported <- pmax(res$p_val_adj_lineage, 2.2e-16)

  res %>% arrange(lineage, cell_type)
}


#' Format an enrichment table for publication as a supplementary table
format_enrichment_table <- function(stats_df) {
  stats_df %>%
    transmute(
      Lineage                      = lineage,
      `Cell subpopulation`         = cell_type,
      Comparison                   = comparison,
      `n (group 1)`                = n_group1,
      `n (group 2)`                = n_group2,
      `Weighted mean % (group 1)`  = signif(wmean_group1, 4),
      `Weighted s.e.m. (group 1)`  = signif(wse_group1, 3),
      `Weighted mean % (group 2)`  = signif(wmean_group2, 4),
      `Weighted s.e.m. (group 2)`  = signif(wse_group2, 3),
      `Fold-enrichment`            = signif(fold_enrichment, 4),
      `Fold-enrichment 95% CI`     = paste0(signif(fold_ci_low, 3), " to ",
                                            signif(fold_ci_high, 3)),
      `t`                          = round(t_value, 3),
      `d.f.`                       = round(df, 1),
      `P (two-sided)`              = signif(p_val_raw_reported, 3),
      `P (BH, within lineage)`     = signif(p_val_adj_lineage_reported, 3),
      `P (BH, all subpopulations)` = signif(p_val_adj_global, 3)
    )
}


###############################################################################
## Input
###############################################################################

if (!exists("ATLAS")) {
  if (!exists("obj")) {
    stop("No atlas found. Load the integrated Seurat object as `obj`, or set ",
         "ATLAS to a Seurat object or an .rds path, then source this file ",
         "again.")
  }
  ATLAS <- obj
}

if (!exists("OUTDIR"))  OUTDIR  <- "."
if (!exists("OUTFILE")) OUTFILE <- "subpopulation_enrichment.csv"
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)

meta  <- seurat_meta(ATLAS, "Atlas")
acols <- annotation_cols(meta)

## Only unsorted samples contribute to the enrichment analysis. Samples that
## underwent FACS-based enrichment or depletion have compositions determined by
## the sorting strategy rather than by the tissue, which would bias any
## comparison of subpopulation proportions.
if ("sorting" %in% names(meta)) {
  meta <- meta[meta$sorting == "none", , drop = FALSE]
}

enrich <- function(m, group_col, group1, group2) {
  calculate_enrichment(m, group_col, group1, group2,
                       lineage_col   = acols["lineage"],
                       cell_type_col = acols["cell_type"])
}

results <- list()


###############################################################################
## Comparisons
###############################################################################

## --- Tumor vs non-tumor tissue ----------------------------------------------
## Adjacent and normal tissue are pooled as the non-tumor comparator; the group
## takes whichever of the two names the atlas actually uses.
if ("site" %in% names(meta)) {

  site <- tolower(as.character(meta$site))
  non_tumor <- if (any(site == "adjacent")) "Adjacent" else "Normal"

  meta$bin_site <- ""
  meta$bin_site[site %in% c("adjacent", "normal")] <- non_tumor
  meta$bin_site[site == "tumor"]                   <- "Tumor"

  results$site <- enrich(meta, "bin_site", "Tumor", non_tumor)
}

## --- Advanced vs early stage ------------------------------------------------
## Restricted to tumor samples with a recorded stage and, where treatment
## response is recorded, to treatment-naive samples, so that stage effects are
## not confounded by neoadjuvant treatment.
if ("clinical_stage" %in% names(meta) && "site" %in% names(meta)) {

  meta$bin_stage <- ""
  meta$bin_stage[meta$clinical_stage %in% c(1, 2)]                    <- "Stage 1/2"
  meta$bin_stage[meta$clinical_stage %in% c(3, 4, "3_4", "3_or_4")]   <- "Stage 3/4"

  m <- meta %>% filter(tolower(site) == "tumor", bin_stage != "")
  if ("treatment_response" %in% names(m)) {
    m <- m %>% filter(!treatment_response %in% c("Y", "N"))
  }

  results$stage <- enrich(m, "bin_stage", "Stage 3/4", "Stage 1/2")
}

## --- Treatment non-responders vs responders ---------------------------------
if ("treatment_response" %in% names(meta) && "site" %in% names(meta)) {

  meta$bin_treatment <- ""
  meta$bin_treatment[meta$treatment_response == "Y"] <- "Treatment responders"
  meta$bin_treatment[meta$treatment_response == "N"] <- "Treatment non-responders"

  m <- meta %>% filter(tolower(site) == "tumor", bin_treatment != "")

  results$treatment <- enrich(m, "bin_treatment",
                              "Treatment non-responders",
                              "Treatment responders")
}

if (length(results) == 0) {
  stop("None of the comparisons could be run. The atlas carries none of the ",
       "meta-data columns 'site', 'clinical_stage' or 'treatment_response'.")
}


###############################################################################
## Results table
###############################################################################

write.csv(bind_rows(lapply(results, format_enrichment_table)),
          file.path(OUTDIR, OUTFILE), row.names = FALSE)

message("Wrote ", file.path(OUTDIR, OUTFILE), " (",
        paste(names(results), collapse = ", "), ")")


###############################################################################
## Samples contributing to each lineage
###############################################################################

## Not every sample captures every lineage, so these differ between lineages.
for (nm in names(results)) {
  d <- results[[nm]]
  if (nrow(d) == 0) next
  cat("\n", nm, ":\n", sep = "")
  print(d %>%
          group_by(lineage) %>%
          summarise(n_group1 = unique(n_group1),
                    n_group2 = unique(n_group2), .groups = "drop") %>%
          as.data.frame())
}
