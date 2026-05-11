## ===============================================================
## Shared helpers for EMM table scripts.
##
## - BH-FDR adjustment applied within (timepoint, outcome-class)
##   families, separately for stratum-specific RD p-values and
##   for interaction p-values.
## - Lookup tables for modifier-level labels so behavioral
##   modifiers (diarrhea, fever) render as "No"/"Yes" while
##   pathogen indicators render as "Negative"/"Positive".
## ===============================================================

## BH-FDR adjustment over a numeric vector, preserving NAs.
bh_fdr <- function(p) {
  out <- rep(NA_real_, length(p))
  idx <- !is.na(p)
  if (any(idx)) out[idx] <- stats::p.adjust(p[idx], method = "BH")
  out
}

## Modifier-level labels for the pathogen/symptom EMM tables.
## V is the column name in src/09 output; lev0/lev1 are the
## labels for the 0/1 levels of that indicator.
modifier_level_lut <- tibble::tribble(
  ~V,                        ~lev0,       ~lev1,
  "ch_pos_giardia",          "Negative",  "Positive",
  "ch_pos_entamoeba",        "Negative",  "Positive",
  "ch_pos_crypto",           "Negative",  "Positive",
  "ch_qpcr_pos_trichuris",   "Negative",  "Positive",
  "ch_qpcr_pos_ascaris",     "Negative",  "Positive",
  "diar7d_t2",               "No",        "Yes",
  "diar7d_t3",               "No",        "Yes",
  "fever7d_t2",              "No",        "Yes",
  "fever7d_t3",              "No",        "Yes"
)

## Look up the label for a given (V, raw-level) pair.
## Falls back to the raw level if V is not in the LUT.
modifier_level_label <- function(V, level) {
  level <- as.character(level)
  row <- modifier_level_lut[modifier_level_lut$V == V, , drop = FALSE]
  if (nrow(row) == 0) return(level)
  vapply(level, function(x) {
    if (x == "0") row$lev0
    else if (x == "1") row$lev1
    else x
  }, character(1))
}

## Format a p-value with the project's standard conventions.
fmt_p_emm <- function(p) {
  ifelse(is.na(p), "", format.pval(p, digits = 2, eps = 0.001))
}
