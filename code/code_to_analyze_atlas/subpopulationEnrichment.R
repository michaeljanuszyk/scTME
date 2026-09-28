###############################################################################
## subpopulationEnrichment.R
##
## Enrichment of cell subpopulations between groups of samples defined by
## sample-level meta-data, for an integrated single-cell atlas.
##
## This script reproduces Supplementary Table 8 exactly.
##
## Usage
## -----
##   obj <- readRDS("allCells_final_meta.rds")
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
## Fold-enrichment is the ratio of the two groups' weighted mean proportions,
## group 1 over group 2, matching the orientation of Supplementary Table 8 and
## of the published forest plots. Its standard error is obtained by propagating
## the two weighted standard errors, and its 95% confidence interval by the
## delta method on the log scale, so the interval cannot extend below zero.
##
## Significance is assessed by a two-sided weighted Welch t-test for two
## independent samples on the per-sample proportions, giving a t statistic and
## Welch-Satterthwaite degrees of freedom. Because subpopulations are nested
## within parent lineages, Benjamini-Hochberg adjustment is applied separately
## within each lineage; a global adjustment across all subpopulations is also
## reported for completeness.
##
## Dependencies: dplyr, tidyr, weights
##
## Januszyk M, Lu JM, Longaker MT. Conservation of the tumor microenvironment
## across species and organs.
###############################################################################

library(dplyr)
library(tidyr)
library(weights)


###############################################################################
## Analysis parameters
###############################################################################

## Lineages excluded from the analysis. Epithelial and nerve cells are outside
## the scope of the stromal/immune comparison.
EXCLUDED_LINEAGES <- c("Ner", "Epi")

## Subpopulations excluded within the retained lineages. Osteoclasts, Kupffer
## cells and proliferating cells are organ-restricted or cell-cycle defined
## rather than transcriptional TME states.
##
## Note for anyone editing this pattern: an empty alternative (a "||", or a
## leading or trailing "|") matches every string, which would silently discard
## the entire dataset. Remove a term completely rather than blanking it.
EXCLUDED_SUBPOP_REGEX <- "Osteoclast|Epithelial|Kupffer|Proliferating|Ner"

## Minimum cells a sample must contribute to a lineage to enter that lineage's
## analysis. Proportions estimated from very few cells are too noisy to be
## informative. Applied independently per lineage, so the number of
## contributing samples differs between lineages.
MIN_CELLS_PER_LINEAGE <- 30


###############################################################################
## Functions
###############################################################################

#' Per-cell metadata from a Seurat object
#'
#' @param x    Seurat object, path to an .rds holding one, or a data.frame
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
#' pair is present rather than hardcoding one.
annotation_cols <- function(meta) {
  pairs <- list(
    c(lineage = "integrated.annotation.abbr", cell_type = "integrated.annotation"),
    c(lineage = "cellType_fine_abbr",         cell_type = "cellType_fine")
  )
  for (p in pairs) if (all(p %in% names(meta))) return(p)
  stop("Could not find a lineage/subpopulation annotation column pair.\n",
       "Available columns: ", paste(names(meta), collapse = ", "))
}


#' Weighted mean and standard error using the Kish effective sample size
#'
#' @param x numeric vector of per-sample subpopulation percentages
#' @param w numeric vector of sample weights (sqrt of cells per sample)
#' @return named numeric vector: wavg, wse, n_eff
weighted_stats <- function(x, w) {
  if (length(x) < 2 || sum(w) == 0) {
    return(c(wavg = NA_real_, wse = NA_real_, n_eff = NA_real_))
  }
  wavg  <- sum(w * x) / sum(w)
  n_eff <- (sum(w))^2 / sum(w^2)
  if (n_eff <= 1) return(c(wavg = wavg, wse = NA_real_, n_eff = n_eff))
  ## Bias correction for the weighted variance, analogous to the n/(n-1)
  ## correction of the unweighted sample variance.
  correction <- n_eff / (n_eff - 1)
  wsd <- sqrt(correction * sum(w * (x - wavg)^2) / sum(w))
  c(wavg = wavg, wse = wsd / sqrt(n_eff), n_eff = n_eff)
}


#' Enrichment of every cell subpopulation between two sample groups
#'
#' @param meta          a Seurat object, or its per-cell metadata data.frame
#' @param group_col     metadata column defining the two groups
#' @param group1        numerator group of the fold-enrichment
#' @param group2        denominator group
#' @param sample_col    metadata column identifying the biological sample
#' @param lineage_col   metadata column giving the parent lineage
#' @param cell_type_col metadata column giving the subpopulation
#' @param min_cells     per-lineage inclusion threshold
#'
#' @return One row per subpopulation, carrying the exact n per group, the
#'   weighted group means and standard errors, the fold-enrichment (group 1 /
#'   group 2) with its standard error and 95% confidence interval, the t
#'   statistic, the Welch-Satterthwaite degrees of freedom, and raw and
#'   adjusted P values.
calculate_enrichment <- function(meta,
                                 group_col, group1, group2,
                                 sample_col    = "sample",
                                 lineage_col   = "integrated.annotation.abbr",
                                 cell_type_col = "integrated.annotation",
                                 min_cells     = MIN_CELLS_PER_LINEAGE) {

  if (inherits(meta, "Seurat")) meta <- meta@meta.data

  missing <- setdiff(c(sample_col, group_col, lineage_col, cell_type_col),
                     names(meta))
  if (length(missing)) {
    stop("Column(s) not found: ", paste(missing, collapse = ", "),
         "\nAvailable columns: ", paste(names(meta), collapse = ", "))
  }

  ## Cells without a lineage or subpopulation assignment are dropped.
  filtered <- meta %>%
    filter(!is.na(.data[[lineage_col]]),
           !is.na(.data[[cell_type_col]]),
           !(.data[[lineage_col]] %in% EXCLUDED_LINEAGES),
           !grepl(EXCLUDED_SUBPOP_REGEX, .data[[cell_type_col]]),
           .data[[group_col]] %in% c(group1, group2))

  if (nrow(filtered) == 0) {
    stop("No cells remain after filtering for '", group1, "' vs '", group2, "'.")
  }

  out <- list()

  for (lineage in unique(filtered[[lineage_col]])) {

    bc <- filtered %>% filter(.data[[lineage_col]] == lineage)

    ## Samples contributing `min_cells` or fewer cells to this lineage are
    ## dropped, independently for each lineage.
    counts <- table(bc[[sample_col]])
    bc <- bc %>% filter(.data[[sample_col]] %in%
                          names(counts)[counts > min_cells])
    if (nrow(bc) == 0) next

    sample_totals <- bc %>%
      group_by(sample = .data[[sample_col]], group = .data[[group_col]]) %>%
      summarise(total_cells = n(), .groups = "drop")

    cell_counts <- bc %>%
      group_by(sample    = .data[[sample_col]],
               group     = .data[[group_col]],
               cell_type = .data[[cell_type_col]]) %>%
      summarise(cell_type_count = n(), .groups = "drop")

    ## Subpopulations absent from a sample are retained as zero rather than
    ## missing, so every sample contributes a value for every subpopulation.
    props <- cell_counts %>%
      right_join(expand_grid(sample_totals,
                             cell_type = unique(bc[[cell_type_col]])),
                 by = c("sample", "group", "cell_type")) %>%
      mutate(cell_type_count = replace_na(cell_type_count, 0),
             proportion      = (cell_type_count / total_cells) * 100,
             weight          = sqrt(total_cells)) %>%
      filter(weight > 0)

    rows <- list()

    for (ct in unique(props$cell_type)) {

      g1 <- props %>% filter(cell_type == ct, group == group1)
      g2 <- props %>% filter(cell_type == ct, group == group2)

      s1 <- weighted_stats(g1$proportion, g1$weight)
      s2 <- weighted_stats(g2$proportion, g2$weight)

      ## Fold-enrichment, group 1 over group 2. se_log is the standard error
      ## of log(fold) and gives the delta-method confidence interval.
      fold   <- as.numeric(s1["wavg"] / s2["wavg"])
      se_log <- as.numeric(sqrt((s1["wse"] / s1["wavg"])^2 +
                                (s2["wse"] / s2["wavg"])^2))

      ## Two-sided weighted Welch t-test on the per-sample proportions. A
      ## positive t means group 1 exceeds group 2, matching fold > 1.
      t_value <- df <- p_raw <- NA_real_
      if (nrow(g1) >= 2 && nrow(g2) >= 2) {
        tt <- wtd.t.test(x = g1$proportion, weight  = g1$weight,
                         y = g2$proportion, weighty = g2$weight)
        t_value <- as.numeric(tt$coefficients["t.value"])
        df      <- as.numeric(tt$coefficients["df"])
        p_raw   <- as.numeric(tt$coefficients["p.value"])
      }

      rows[[ct]] <- data.frame(
        comparison       = paste0(group1, " vs ", group2),
        lineage          = lineage,
        cell_type        = ct,
        n_group1         = nrow(g1),
        n_group2         = nrow(g2),
        wmean_group1     = as.numeric(s1["wavg"]),
        wse_group1       = as.numeric(s1["wse"]),
        wmean_group2     = as.numeric(s2["wavg"]),
        wse_group2       = as.numeric(s2["wse"]),
        fold_enrichment  = fold,
        fold_enrich_se   = abs(fold) * se_log,
        fold_ci_low      = exp(log(fold) - 1.96 * se_log),
        fold_ci_high     = exp(log(fold) + 1.96 * se_log),
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

  res$p_val_adj_global <- p.adjust(res$p_val_raw, method = "BH")

  ## P values below double precision are reported at that floor.
  res$p_val_raw         <- pmax(res$p_val_raw,         2.2e-16)
  res$p_val_adj_lineage <- pmax(res$p_val_adj_lineage, 2.2e-16)
  res$p_val_adj_global  <- pmax(res$p_val_adj_global,  2.2e-16)

  res %>% arrange(lineage, cell_type)
}


#' Format an enrichment table for publication as a supplementary table
format_enrichment_table <- function(stats_df) {
  stats_df %>%
    transmute(
      Comparison                   = comparison,
      Lineage                      = lineage,
      `Cell subpopulation`         = cell_type,
      `n (group 1)`                = n_group1,
      `n (group 2)`                = n_group2,
      `Weighted mean % (group 1)`  = signif(wmean_group1, 4),
      `Weighted s.e.m. (group 1)`  = signif(wse_group1, 3),
      `Weighted mean % (group 2)`  = signif(wmean_group2, 4),
      `Weighted s.e.m. (group 2)`  = signif(wse_group2, 3),
      `Fold-enrichment`            = signif(fold_enrichment, 4),
      `95% CI (lower)`             = signif(fold_ci_low, 3),
      `95% CI (upper)`             = signif(fold_ci_high, 3),
      `t`                          = round(t_value, 3),
      `d.f.`                       = round(df, 1),
      `P (two-sided, unadjusted)`  = signif(p_val_raw, 3),
      `P (BH, within lineage)`     = signif(p_val_adj_lineage, 3),
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

enrich <- function(m, group_col, group1, group2) {
  calculate_enrichment(m, group_col, group1, group2,
                       lineage_col   = acols[["lineage"]],
                       cell_type_col = acols[["cell_type"]])
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
if (all(c("site", "clinical_stage") %in% names(meta))) {

  meta$bin_stage <- ""
  meta$bin_stage[meta$clinical_stage %in% c(1, 2)]                  <- "Stage 1/2"
  meta$bin_stage[meta$clinical_stage %in% c(3, 4, "3_4", "3_or_4")] <- "Stage 3/4"

  m <- meta %>% filter(tolower(site) == "tumor", bin_stage != "")
  if ("treatment_response" %in% names(m)) {
    m <- m %>% filter(!treatment_response %in% c("Y", "N"))
  }

  results$stage <- enrich(m, "bin_stage", "Stage 3/4", "Stage 1/2")
}

## --- Treatment non-responders vs responders ---------------------------------
if (all(c("site", "treatment_response") %in% names(meta))) {

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

## Not every sample captures every lineage, and the per-lineage cell threshold
## removes different samples in each, so these differ between lineages.
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
