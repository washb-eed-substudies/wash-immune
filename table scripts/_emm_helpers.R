## ===============================================================
## Shared helpers for the EMM (effect-modification) tables.
##
## Both consumers source this file:
##   - table scripts/path_subgroup_tables.R  (standalone docx per timepoint)
##   - manuscript/supplement.qmd             (full assembled supplement)
## Keep all table-building logic here so the two can never drift.
##
## For each (timepoint, modifier) we build one wide flextable with a
## 3-row block per cytokine ratio:
##   1. Outcome header row carrying subgroup-specific treatment effects
##      (RD, 95% CI, p-value) and the interaction p-value.
##   2. Control row with N, mean, SD broken out by subgroup level.
##   3. Nutrition+WSH row with N, mean, SD broken out by subgroup level.
##
## MULTIPLE TESTING: both the raw p-value and a p-value adjusted to
## control the family-wise error rate are reported, as specified in the
## pre-analysis plan. The adjustment is the project's Bonferroni
## convention - raw p x 2, capped at 1 - identical to bonpval() in
## "table scripts/immune-main-tables.R", so the EMM tables, the main
## tables and S2-S8 all correct the same way. This replaced a BH-FDR
## adjustment that was inconsistent with the pre-registered plan.
## Both the stratum-specific (per-subgroup) RD p-values and the
## interaction p-values are reported raw and corrected; the hypothesis
## under test in an EMM analysis is the interaction, so the per-subgroup
## effects remain descriptive and should be read alongside their 95% CIs.
##
## SPECIMEN TIMING (see the header of src/09 for the full audit):
##   14-month outcomes -> 14-month TAC qPCR modifiers (same-day stool)
##   28-month outcomes -> BOTH families, labelled distinctly: the
##                        14-month TAC specimen (13.9 months before the
##                        draw) and the ~30-month specimen (2.0 months
##                        after it)
##   change score      -> 14-month TAC qPCR modifiers (precede T3)
## The ~30-month specimen post-dates the 14-month blood draw by a median
## of 16.0 months and must not be used as a modifier there.
##
## NOTE: these tables carry no explanatory preamble. Each one footnotes
## its own clustered standard errors, p-value correction and cytokine
## group definitions, and the heading plus the specimen note below it
## carry the specimen, the age and the timing relative to the outcome.
## Everything else - the minimum-cell exclusion rule and the
## post-randomization caveat - has to be stated in the manuscript
## Methods instead.
## ===============================================================

library(here)
library(dplyr)
library(stringr)
library(tidyr)
library(flextable)
library(officer)

## Must match MIN_CELL in src/09 - used for the footnote wording.
MIN_CELL_REPORTED <- 10

## ---------------------------------------------------------------
## Small formatters
## ---------------------------------------------------------------

## Bonferroni correction, project convention: the raw p-value is doubled
## and capped at 1. This is the same rule bonpval() applies in
## "table scripts/immune-main-tables.R" for the main tables and S2-S8.
BONFERRONI_FACTOR <- 2
bonferroni_p <- function(p) pmin(1, p * BONFERRONI_FACTOR)

## Display rule for p-values, also lifted from bonpval() so that every
## table in the paper renders p-values the same way. Vectorized, and NA
## renders blank.
fmt_p_proj <- function(p) {
  ifelse(is.na(p), "",
    ifelse(p < 0.001, "<0.001",
      ifelse(p < 0.01, "<0.01", as.character(round(p, 2)))))
}

fmt_num <- function(x) ifelse(is.na(x), "", sprintf("%.2f", x))
fmt_int <- function(x) ifelse(is.na(x), "", as.character(as.integer(x)))

## ---------------------------------------------------------------
## Modifier metadata
## ---------------------------------------------------------------

## Level labels: behavioural modifiers render as "No"/"Yes", pathogen
## indicators as "Negative"/"Positive".
modifier_level_lut <- tibble::tribble(
  ~V,                        ~lev0,       ~lev1,
  ## ~30-month year-2 parasite specimen (Kato-Katz + qPCR, same sample)
  "ch_pos_giardia",          "Negative",  "Positive",
  "ch_pos_entamoeba",        "Negative",  "Positive",
  "ch_pos_crypto",           "Negative",  "Positive",
  "ch_pos_trichuris",        "Negative",  "Positive",
  "ch_pos_ascaris",          "Negative",  "Positive",
  ## 14-month TAC multiplex qPCR
  "tac14_pos_giardia",       "Negative",  "Positive",
  "tac14_pos_crypto",        "Negative",  "Positive",
  "tac14_pos_entamoeba",     "Negative",  "Positive",
  "tac14_pos_ascaris",       "Negative",  "Positive",
  "tac14_pos_trichuris",     "Negative",  "Positive",
  "tac14_pos_hookworm",      "Negative",  "Positive",
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

## Pathogen name and specimen for each modifier. The TIMING WORDING is
## NOT stored here: the same 14-month specimen is concurrent with the
## 14-month outcome but prior to the 28-month outcome, so the framing
## has to be resolved per timepoint - see modifier_heading()/modifier_note().
## Row order here drives table order within a timepoint: the 14-month
## specimen family first, then the 30-month family, then symptoms.
modifier_lut <- tibble::tribble(
  ~V,                    ~pathogen,                ~specimen,                 ~specimen_age,
  ## --- 14-month TAC multiplex qPCR
  "tac14_pos_giardia",   "Giardia",                "14-month stool, TAC qPCR",   14L,
  "tac14_pos_crypto",    "Cryptosporidium",        "14-month stool, TAC qPCR",   14L,
  "tac14_pos_entamoeba", "Entamoeba",              "14-month stool, TAC qPCR",   14L,
  "tac14_pos_ascaris",   "Ascaris",                "14-month stool, TAC qPCR",   14L,
  "tac14_pos_trichuris", "Trichuris",              "14-month stool, TAC qPCR",   14L,
  "tac14_pos_hookworm",  "Hookworm",               "14-month stool, TAC qPCR",   14L,
  ## --- ~30-month year-2 parasite specimen (Kato-Katz + qPCR, one sample)
  "ch_pos_giardia",      "Giardia",                "30-month stool, qPCR",       30L,
  "ch_pos_entamoeba",    "Entamoeba histolytica",  "30-month stool, qPCR",       30L,
  "ch_pos_crypto",       "Cryptosporidium",        "30-month stool, qPCR",       30L,
  "ch_pos_trichuris",    "Trichuris",              "30-month stool, Kato-Katz",  30L,
  "ch_pos_ascaris",      "Ascaris",                "30-month stool, Kato-Katz",  30L,
  ## --- caregiver-reported symptoms
  "diar7d_t2",           "Diarrhea, past 7 days",  "14-month visit",             14L,
  "diar7d_t3",           "Diarrhea, past 7 days",  "28-month visit",             28L,
  "fever7d_t2",          "Fever, past 7 days",     "14-month visit",             14L,
  "fever7d_t3",          "Fever, past 7 days",     "28-month visit",             28L
)

## When the modifier was measured relative to the outcome it modifies.
## "prior" / "concurrent" / "subsequent" - the wording readers need to
## judge what the stratification can and cannot support.
modifier_timing <- function(V, timepoint) {
  age <- modifier_lut$specimen_age[match(V, modifier_lut$V)]
  if (is.na(age)) return(NA_character_)
  if (timepoint == "T2")    return(if (age <= 14) "concurrent" else "subsequent")
  if (timepoint == "T3")    return(if (age <= 14) "prior" else if (age <= 28) "concurrent" else "subsequent")
  if (timepoint == "Delta") return(if (age <= 14) "start of interval" else "end of interval")
  NA_character_
}

## Plain-English gloss of when the modifier sits relative to the outcome.
timing_phrase <- function(timing) {
  switch(timing,
    "prior"             = "measured before the outcome",
    "concurrent"        = "measured concurrently with the outcome",
    "subsequent"        = "measured after the outcome",
    "start of interval" = "measured at the start of the interval",
    "end of interval"   = "measured at the end of the interval",
    NA_character_)
}

## Table heading text: what was measured, in which specimen/visit, and
## where it sits relative to the outcome. Stool modifiers name the assay;
## symptom modifiers just name the visit, so they take a different shape
## to avoid reading as "Diarrhea ... concurrent infection".
modifier_heading <- function(V, timepoint) {
  i <- match(V, modifier_lut$V)
  if (is.na(i)) return(V)
  timing <- timing_phrase(modifier_timing(V, timepoint))
  if (grepl("stool", modifier_lut$specimen[i])) {
    paste0(modifier_lut$pathogen[i], " (", modifier_lut$specimen[i],
           if (!is.na(timing)) paste0("; ", timing) else "", ")")
  } else {
    paste0(modifier_lut$pathogen[i], " at the ", modifier_lut$specimen[i],
           if (!is.na(timing)) paste0(" (", timing, ")") else "")
  }
}

## Footnote under each table, stating the measured gap to the blood draw.
modifier_note <- function(V, timepoint) {
  i <- match(V, modifier_lut$V)
  if (is.na(i)) return("")
  age <- modifier_lut$specimen_age[i]
  is_stool <- grepl("stool", modifier_lut$specimen[i])
  if (!is_stool) return("")

  spec14 <- paste0(
    "Infection status from the 14-month stool specimen assayed by TaqMan ",
    "Array Card multiplex qPCR (median age 14.3 months, IQR 12.7-15.7). ")
  spec30 <- paste0(
    "Infection status from the year-2 parasite stool specimen (median age ",
    "30.3 months, IQR 29.2-31.6). A single specimen was assayed both by ",
    "Kato-Katz (Ascaris, Trichuris) and by qPCR (Giardia, Cryptosporidium, ",
    "Entamoeba histolytica); these are the same sample at the same visit, ",
    "not two collections. ")

  gap <- if (age == 14 && timepoint == "T2") {
    paste0("Collected on the same day as the 14-month blood draw (median ",
           "difference 0.0 months), so infection is concurrent with the outcome.")
  } else if (age == 14 && timepoint == "T3") {
    paste0("Collected a median of 13.9 months (IQR 13.4-14.6) BEFORE the ",
           "28-month blood draw, in 100% of children, so infection precedes ",
           "the outcome.")
  } else if (age == 14 && timepoint == "Delta") {
    paste0("Collected at the start of the 14-to-28-month interval, so ",
           "infection precedes the 28-month endpoint of the change score.")
  } else if (age == 30 && timepoint == "T3") {
    paste0("Collected a median of 2.0 months (IQR 1.2-2.7) AFTER the ",
           "28-month blood draw, so this is a cross-sectional association ",
           "with concurrent infection, not a prospective one.")
  } else ""

  paste0(if (age == 14) spec14 else spec30, gap)
}

## ---------------------------------------------------------------
## Which tables get a number
## ---------------------------------------------------------------
## At 28 months two pathogen families are reported for the same
## outcomes: the concurrent ~30-month specimen and the prospective
## 14-month TAC specimen. The concurrent family keeps the table
## numbers (it is the one the manuscript text refers to); the
## prospective family follows at the end of the section as unnumbered
## "Table SX" for comparison. Nothing is unnumbered at T2 or Delta,
## where only one family is reported.
modifier_is_unnumbered <- function(V, timepoint) {
  identical(timepoint, "T3") && grepl("^tac14_", V)
}

## Ordered plan of tables for one timepoint: numbered first (in
## modifier_lut order), then the unnumbered comparison tables.
emm_table_plan <- function(df_tp, timepoint) {
  present <- modifier_lut$V[modifier_lut$V %in% unique(df_tp$V)]
  if (length(present) == 0) {
    return(tibble::tibble(V = character(0), numbered = logical(0)))
  }
  unnum <- vapply(present, modifier_is_unnumbered, logical(1),
                  timepoint = timepoint)
  tibble::tibble(
    V        = c(present[!unnum], present[unnum]),
    numbered = c(rep(TRUE, sum(!unnum)), rep(FALSE, sum(unnum)))
  )
}

## ---------------------------------------------------------------
## Cytokine labels and ordering
## ---------------------------------------------------------------

cytokine_lut <- tibble::tribble(
  ~raw,    ~label,
  "il1",   "Interleukin-1\u03b2",
  "il2",   "Interleukin-2",
  "il4",   "Interleukin-4",
  "il5",   "Interleukin-5",
  "il6",   "Interleukin-6",
  "il10",  "Interleukin-10",
  "il12",  "Interleukin-12",
  "il13",  "Interleukin-13",
  "il17",  "Interleukin-17",
  "il21",  "Interleukin-21",
  "ifn",   "Interferon-\u03b3",
  "tnf",   "Tumor necrosis factor-\u03b1",
  "gmc",   "Granulocyte-macrophage colony-stimulating factor",
  "th1",   "Th1",
  "th2",   "Th2",
  "th17",  "Th17",
  "pro",   "Pro-inflammatory cytokines"
)

clean_ratio_label <- function(stem) {
  ratio <- stem %>% str_remove("^ratio_")
  parts <- str_split(ratio, "_", simplify = TRUE)
  paste0(
    cytokine_lut$label[match(parts[, 1], cytokine_lut$raw)],
    "/",
    cytokine_lut$label[match(parts[, 2], cytokine_lut$raw)]
  )
}

ratio_order <- c(
  "ratio_th1_th2","ratio_il12_il4","ratio_ifn_il4","ratio_il12_il5","ratio_ifn_il5",
  "ratio_il12_il13","ratio_ifn_il13",
  "ratio_pro_il10","ratio_il1_il10","ratio_il6_il10","ratio_tnf_il10","ratio_il2_il10",
  "ratio_th1_il10","ratio_th2_il10","ratio_il12_il10","ratio_ifn_il10","ratio_il4_il10",
  "ratio_il5_il10","ratio_il13_il10","ratio_th17_il10","ratio_il17_il10","ratio_il21_il10",
  "ratio_gmc_il10",
  "ratio_th1_th17","ratio_il12_il17","ratio_ifn_il17","ratio_il12_il21","ratio_ifn_il21"
)

## ---------------------------------------------------------------
## Load and prepare the EMM results
## ---------------------------------------------------------------

emm_result_files <- c(
  T2    = "results/bangladesh-immune-posthoc-subgroup-results-t2.csv",
  T3    = "results/bangladesh-immune-posthoc-subgroup-results-t3.csv",
  Delta = "results/bangladesh-immune-posthoc-subgroup-results-delta.csv"
)

## Returns one long data frame across all available timepoints, with
## cytokine labels attached and a Bonferroni-corrected copy of each
## interaction p-value alongside the raw one.
load_emm_results <- function(files = emm_result_files) {
  read_one <- function(path, tp) {
    fp <- here(path)
    if (!file.exists(fp)) {
      warning("Missing ", path, " - skipping ", tp)
      return(NULL)
    }
    read.csv(fp, stringsAsFactors = FALSE) %>% mutate(timepoint = tp)
  }

  dat <- bind_rows(lapply(names(files),
                          function(tp) read_one(files[[tp]], tp)))

  if (is.null(dat) || nrow(dat) == 0) {
    stop("No EMM result CSVs found. Run src/09 first.")
  }

  ## Normalize the p-value column name that read.csv turns into "P.value"
  if ("P.value" %in% names(dat)) names(dat)[names(dat) == "P.value"] <- "P_value"

  ## Strip the timepoint prefix so one ratio label applies across timepoints
  dat <- dat %>%
    mutate(ratio_stem = Y %>% str_remove("^t[0-9]_") %>% str_remove("^d23_"),
           Y_label    = clean_ratio_label(ratio_stem),
           ratio_stem = factor(ratio_stem, levels = ratio_order)) %>%
    arrange(timepoint, V, ratio_stem, subgroup)

  ## Bonferroni-corrected copy of each interaction p-value. The raw
  ## value stays in int_Pval and both are reported - see header.
  int_p <- dat %>%
    group_by(timepoint, Y, V) %>%
    summarise(int_Pval = first(na.omit(int_Pval)), .groups = "drop") %>%
    mutate(int_P_bonf = bonferroni_p(int_Pval)) %>%
    select(timepoint, Y, V, int_P_bonf)

  dat %>% left_join(int_p, by = c("timepoint", "Y", "V"))
}

## ---------------------------------------------------------------
## Table construction
## ---------------------------------------------------------------

## One 3-row block per cytokine ratio. df_y must have exactly 2 rows,
## one per subgroup level.
build_outcome_block <- function(df_y) {
  df_y <- df_y %>% arrange(subgroup)
  if (nrow(df_y) != 2) return(NULL)

  rd_ci_1  <- sprintf("%.3f (%.3f, %.3f)", df_y$RD[1], df_y$ci.lb[1], df_y$ci.ub[1])
  rd_ci_2  <- sprintf("%.3f (%.3f, %.3f)", df_y$RD[2], df_y$ci.lb[2], df_y$ci.ub[2])
  raw_int  <- first(na.omit(df_y$int_Pval))
  bonf_int <- first(na.omit(df_y$int_P_bonf))

  header_row <- tibble::tibble(
    Outcome = df_y$Y_label[1],
    L1_N = "", L1_Mean = "", L1_SD = "", L1_RD = rd_ci_1,
    L1_P = fmt_p_proj(df_y$P_value[1]),
    L1_Pbonf = fmt_p_proj(bonferroni_p(df_y$P_value[1])),
    L2_N = "", L2_Mean = "", L2_SD = "", L2_RD = rd_ci_2,
    L2_P = fmt_p_proj(df_y$P_value[2]),
    L2_Pbonf = fmt_p_proj(bonferroni_p(df_y$P_value[2])),
    Int_P = fmt_p_proj(raw_int),
    Int_Pbonf = fmt_p_proj(bonf_int)
  )

  control_row <- tibble::tibble(
    Outcome = "Control",
    L1_N = fmt_int(df_y$n_Control[1]),
    L1_Mean = fmt_num(df_y$mean_Control[1]),
    L1_SD = fmt_num(df_y$sd_Control[1]),
    L1_RD = "", L1_P = "", L1_Pbonf = "",
    L2_N = fmt_int(df_y$n_Control[2]),
    L2_Mean = fmt_num(df_y$mean_Control[2]),
    L2_SD = fmt_num(df_y$sd_Control[2]),
    L2_RD = "", L2_P = "", L2_Pbonf = "",
    Int_P = "", Int_Pbonf = ""
  )

  nwsh_row <- tibble::tibble(
    Outcome = "Nutrition + WSH",
    L1_N = fmt_int(df_y$n_NWSH[1]),
    L1_Mean = fmt_num(df_y$mean_NWSH[1]),
    L1_SD = fmt_num(df_y$sd_NWSH[1]),
    L1_RD = "", L1_P = "", L1_Pbonf = "",
    L2_N = fmt_int(df_y$n_NWSH[2]),
    L2_Mean = fmt_num(df_y$mean_NWSH[2]),
    L2_SD = fmt_num(df_y$sd_NWSH[2]),
    L2_RD = "", L2_P = "", L2_Pbonf = "",
    Int_P = "", Int_Pbonf = ""
  )

  bind_rows(header_row, control_row, nwsh_row)
}

## One flextable for a (timepoint, modifier) pair.
build_modifier_table <- function(df_mod) {
  V_code <- unique(df_mod$V)

  ## Subgroup labels in the order they appear after arrange()
  lev_raw <- df_mod %>%
    arrange(ratio_stem, subgroup) %>%
    pull(subgroup) %>%
    as.character() %>%
    unique()
  if (length(lev_raw) != 2) {
    warning("Modifier ", V_code, " has ", length(lev_raw),
            " subgroup levels - skipping.")
    return(NULL)
  }
  lev1_lab <- modifier_level_label(V_code, lev_raw[1])
  lev2_lab <- modifier_level_label(V_code, lev_raw[2])

  body <- df_mod %>%
    group_split(ratio_stem) %>%
    lapply(build_outcome_block) %>%
    bind_rows()

  if (nrow(body) == 0) return(NULL)

  ft <- flextable(body)

  ## Bottom header (column labels) - 6 cols per subgroup
  hdr_bot <- c(
    "Outcome",
    "N", "Mean", "SD", "Difference: Intervention vs. Control (95% CI)",
    "P-value", "Bonferroni P-value",
    "N", "Mean", "SD", "Difference: Intervention vs. Control (95% CI)",
    "P-value", "Bonferroni P-value",
    "P-value", "Bonferroni P-value"
  )
  ## Top header (group spanners) - each subgroup spans 6 cols, and the
  ## raw/corrected interaction p-values share a spanner of their own.
  ## The label is repeated across the cells it covers because merge_h()
  ## below merges runs of identical cells; a blank would not merge.
  hdr_top <- c(
    "Outcome / Arm",
    rep(lev1_lab, 6),
    rep(lev2_lab, 6),
    rep("P-value for interaction", 2)
  )

  ft <- ft %>%
    set_header_labels(values = setNames(as.list(hdr_bot), names(body))) %>%
    add_header_row(values = hdr_top, colwidths = rep(1, length(hdr_top))) %>%
    merge_h(part = "header") %>%
    bold(part = "header") %>%
    bold(j = "Outcome", i = seq(1, nrow(body), by = 3), part = "body") %>%
    align(align = "center", part = "header") %>%
    align(align = "center", part = "body") %>%
    align(align = "left",   j = "Outcome", part = "body") %>%
    valign(valign = "center", part = "all") %>%
    theme_booktabs() %>%
    font(fontname = "Times New Roman", part = "all") %>%
    fontsize(size = 6, part = "body") %>%
    fontsize(size = 7, part = "header") %>%
    padding(padding = 1, part = "all")

  ft <- block_structure(ft, body$Outcome)

  ## Same footnote block as the sex EMM tables, minus the abbreviation
  ## list - these outcome labels are spelled out rather than abbreviated.
  ft <- attach_notes(
    ft,
    list(lead  = NOTE_CLUSTER,
         keyed = list(list(text    = NOTE_P_BOTH,
                           header1 = "^P-value for interaction$",
                           header2 = "^(P-value|Bonferroni P-value)$")),
         groups = TRUE),
    hdr_top, hdr_bot, body$Outcome)

  fit_columns(ft, hdr_bot)
}


## ---------------------------------------------------------------
## Column layout
## ---------------------------------------------------------------
## Left to autofit, Word picks the column widths itself, and with 15-18
## columns it picks badly: the estimate and its CI get wrapped onto three
## lines while N and Mean keep room they never use. Widths are assigned
## here from what each column holds, and the layout is switched to fixed
## so Word honours them.

## Usable text width of reference-landscape.docx: an 11in page less the
## 0.6in margins on each side.
PAGE_WIDTH_IN <- 9.8

## Relative width per column, set so an estimate and its interval stay on
## one line and the outcome labels are what wrap instead.
col_weight <- function(label) {
  l <- trimws(label)
  if (identical(l, "N")) return(0.5)
  if (l %in% c("Mean", "SD")) return(0.7)
  ## Wide enough for "Absolute" on one line - it cannot be broken.
  if (l %in% c("Absolute Mean", "Absolute SD")) return(0.85)
  if (grepl("95% CI|Difference", l)) return(2.1)
  if (identical(l, "P-value")) return(0.8)
  if (identical(l, "Bonferroni P-value")) return(1.0)
  1.0
}

## `first_col` is a weight, not inches - the outcome column carries the
## longest strings in every one of these tables.
fit_columns <- function(ft, hdr_bot, first_col = 2.8) {
  w <- vapply(hdr_bot, col_weight, numeric(1), USE.NAMES = FALSE)
  w[1] <- first_col
  w <- w / sum(w) * PAGE_WIDTH_IN
  for (j in seq_along(w)) ft <- width(ft, j = j, width = w[j])
  set_table_properties(ft, layout = "fixed",
                       opts_word = list(split = FALSE))
}

## These tables run to 85 rows of near-identical numbers. A hairline
## under each outcome block, and keeping the three rows of a block on one
## page, is most of what makes them readable.
block_structure <- function(ft, labels) {
  ends <- which(labels == "Nutrition + WSH")
  if (!length(ends)) return(ft)
  ft %>%
    hline(i = ends, part = "body",
          border = fp_border(color = "#BFBFBF", width = 0.5)) %>%
    keep_with_next(i = setdiff(seq_along(labels), ends), part = "body",
                   value = TRUE) %>%
    keep_with_next(i = ends, part = "body", value = FALSE)
}

## ---------------------------------------------------------------
## Table footnotes
## ---------------------------------------------------------------
## The notes below are verbatim from the supplement circulated to
## coauthors, with one deliberate exception: the p-value note on the
## effect-modification tables. Those tables now show the raw and the
## Bonferroni-corrected p-value side by side, so the old wording
## ("P-values shown are Bonferroni corrected") would be wrong.
##
## Markers are real superscripts attached to the cells they annotate,
## and the cells are located by matching their text rather than by
## column index, so the notes survive columns being added or reordered -
## which is exactly what happened when the raw p-value columns went in.

NOTE_CLUSTER <- paste(
  "Confidence intervals were adjusted for clustered observations using",
  "robust standard errors.")

NOTE_ABBREV <- paste0(
  "IL = Interleukin; TNF-\u03b1 = Tumor necrosis factor-\u03b1; ",
  "CRP = C-reactive protein (CRP); IFN-\u03b3 = interferon-\u03b3; ",
  "GM-CSF = Granulocyte-macrophage colony-stimulating factor; ",
  "AGP = Alpha-1-acid glycoprotein; IGF-1 = Insulin-like growth factor-1")

NOTE_ABBREV_RATIO <- paste0(
  "IL = Interleukin; TNF-\u03b1 = Tumor necrosis factor-\u03b1; ",
  "IFN-\u03b3 = interferon-\u03b3; ",
  "GM-CSF = Granulocyte-macrophage colony-stimulating factor")

NOTE_COVARIATES <- paste0(
  "Adjusted for pre-specified covariates: Field staff who collected ",
  "data, month of measurement, household food insecurity, child age, ",
  "child sex, mother\u2019s age, mother\u2019s height, mother\u2019s ",
  "education level, number of children <18 years in the household, ",
  "number of individuals living in the compound, distance in minutes to ",
  "the primary water source, household floor materials, household wall ",
  "materials, household electricity, and household assets (wardrobe, ",
  "table, chair, clock, khat, chouki, radio, television, refrigerator, ",
  "bicycle, motorcycle, sewing machine, mobile phone, cattle, goats, ",
  "and chickens).")

NOTE_IPCW <- "Inverse probability of censoring weighting."

## Reworded - see the note at the top of this section. Every table that
## reports a p-value now reports the raw one beside the corrected one,
## as the pre-analysis plan calls for, so the old wording ("P-values
## shown are Bonferroni corrected") no longer applies anywhere.
NOTE_P_BOTH <- paste(
  "Both the raw p-value and a p-value Bonferroni corrected to control",
  "for familywise error rates are shown.")

## Cytokine-group definitions, in the order their markers are assigned.
## The tokens are matched by regular expression against the outcome
## labels, so the pattern is built longest-first to keep Th17 from being
## swallowed by Th1.
CYTOKINE_GROUPS <- tibble::tribble(
  ~token,                       ~note,
  "Pro-inflammatory cytokines", "Pro-inflammatory cytokines: IL-1\u03b2, IL-6, TNF-\u03b1",
  "Th1",                        "Th1 cytokines: IL-12, IFN-\u03b3",
  "Th2",                        "Th2 cytokines: IL-4, IL-5, IL-13",
  "Th17",                       "Th17 cytokines: IL-17A, IL-21"
)

## Legacy table keys, by shape. NOTE: these are the keys the old scripts
## use, not the published table numbers - legacy S6/S7 are swapped
## relative to the supplement, as documented in supplement.qmd.
RATIO_TABLES   <- c("S3", "S5", "S7", "S10", "S12")
EMM_SEX_TABLES <- c("S9", "S10", "S11", "S12")

## Break an outcome label into chunks with each cytokine-group token
## followed by its superscript marker. Returns NULL when the label
## mentions no group.
group_chunks <- function(label, syms) {
  toks <- names(syms)[order(nchar(names(syms)), decreasing = TRUE)]
  m <- gregexpr(paste(toks, collapse = "|"), label)[[1]]
  if (m[1] == -1) return(NULL)

  starts <- as.integer(m)
  lens   <- attr(m, "match.length")
  out <- list()
  pos <- 1L
  for (k in seq_along(starts)) {
    if (starts[k] > pos)
      out <- c(out, list(as_chunk(substr(label, pos, starts[k] - 1L))))
    tok <- substr(label, starts[k], starts[k] + lens[k] - 1L)
    out <- c(out, list(as_chunk(tok), as_sup(syms[[tok]])))
    pos <- starts[k] + lens[k]
  }
  if (pos <= nchar(label))
    out <- c(out, list(as_chunk(substr(label, pos, nchar(label)))))
  out
}

## Attach a footnote block to a flextable. hdr1/hdr2/labels are the
## character vectors the table was built from and are used only to find
## the cells that need a marker.
attach_notes <- function(ft, notes, hdr1, hdr2, labels) {
  if (is.null(notes)) return(ft)

  lines <- lapply(notes$lead, function(t) list(sym = NA_character_, text = t))
  sym_i <- 0L
  next_sym <- function() { sym_i <<- sym_i + 1L; letters[sym_i] }

  ## Cells matching a keyed note get its marker as a superscript.
  for (k in notes$keyed) {
    sym <- next_sym()
    for (spec in list(list(p = k$header1, row = 1L, txt = hdr1),
                      list(p = k$header2, row = 2L, txt = hdr2))) {
      if (is.null(spec$p)) next
      for (j in which(grepl(spec$p, spec$txt)))
        ft <- compose(ft, i = spec$row, j = j, part = "header",
                      value = as_paragraph(as_chunk(spec$txt[j]), as_sup(sym)))
    }
    if (!is.null(k$body)) {
      for (i in which(grepl(k$body, labels)))
        ft <- compose(ft, i = i, j = 1, part = "body",
                      value = as_paragraph(as_chunk(labels[i]), as_sup(sym)))
    }
    lines <- c(lines, list(list(sym = sym, text = k$text)))
  }

  ## Cytokine-group markers go inline wherever the group is named.
  if (isTRUE(notes$groups)) {
    present <- CYTOKINE_GROUPS[vapply(CYTOKINE_GROUPS$token, function(t)
      any(grepl(t, labels, fixed = TRUE)), logical(1)), , drop = FALSE]
    if (nrow(present)) {
      syms <- setNames(vapply(seq_len(nrow(present)),
                              function(i) next_sym(), character(1)),
                       present$token)
      for (i in seq_along(labels)) {
        ch <- group_chunks(labels[i], syms)
        if (!is.null(ch))
          ft <- compose(ft, i = i, j = 1, part = "body",
                        value = do.call(as_paragraph, ch))
      }
      for (i in seq_len(nrow(present)))
        lines <- c(lines, list(list(sym  = syms[[present$token[i]]],
                                    text = present$note[i])))
    }
  }

  for (ln in lines) {
    ft <- add_footer_lines(ft, values = if (is.na(ln$sym))
      as_paragraph(as_chunk(ln$text))
    else
      as_paragraph(as_sup(ln$sym), as_chunk(ln$text)))
  }

  ft %>%
    merge_h(part = "footer") %>%
    font(fontname = "Times New Roman", part = "footer") %>%
    fontsize(size = 6, part = "footer") %>%
    align(align = "left", part = "footer") %>%
    padding(padding = 1, part = "footer")
}

## The footnote block for one legacy table, keyed on the legacy name.
legacy_notes <- function(key) {
  if (key == "S1") {
    return(list(
      lead = paste0(
        "Data are n (%) or mean (SD). Percentages were estimated from ",
        "slightly smaller denominators than those shown at the top of the ",
        "table for the following variables due to missing values: ",
        "mother\u2019s age, father\u2019s education, father works in ",
        "agriculture, acres of land owned, open defecation, latrine has a ",
        "concrete slab, latrine has a functional water seal, visible stool ",
        "on latrine slab or floor, ownership of child potty, observed ",
        "feces in the house or child\u2019s play area, handwashing ",
        "variables."),
      keyed = list(
        list(text = "Open defecation does not include diaper disposal of feces.",
             body = "^Reported daily open defecation$"),
        list(text = paste("Households who do not own a latrine typically share",
                          "a latrine with extended family members who live in",
                          "the same compound."),
             body = "^Owned$"),
        list(text = "Assessed by the Household Food Insecurity Access Scale.",
             body = "^Household is food secure$")
      ),
      groups = FALSE))
  }

  is_ratio <- key %in% RATIO_TABLES
  is_sex   <- key %in% EMM_SEX_TABLES

  ## S8 is the sum score - it has no cytokine abbreviations to define.
  lead <- NOTE_CLUSTER
  if (key != "S8")
    lead <- c(lead, if (is_ratio) NOTE_ABBREV_RATIO else NOTE_ABBREV)

  keyed <- if (is_sex) {
    list(list(text    = NOTE_P_BOTH,
              header1 = "^P-value for interaction$",
              header2 = "^(P-value|Bonferroni P-value)$"))
  } else {
    list(
      list(text = NOTE_COVARIATES, header1 = "^Fully adjusted difference"),
      list(text = NOTE_IPCW,       header1 = "^IPCW adjusted difference"),
      list(text = NOTE_P_BOTH,     header2 = "^(P-value|Bonferroni P-value)$")
    )
  }

  list(lead = lead, keyed = keyed, groups = is_ratio)
}

## ---------------------------------------------------------------
## Legacy table conversion
## ---------------------------------------------------------------
## Tables S1-S12 are built by the older scripts as data.tables whose
## COLUMN NAMES are the spanner row and whose FIRST DATA ROW is the
## column-label row. Convert that shape into a two-row-header flextable.

legacy_flextable <- function(tbl, fontsize_body = 7, fontsize_header = 8,
                             notes = NULL) {
  df <- as.data.frame(tbl, stringsAsFactors = FALSE)
  df[] <- lapply(df, as.character)

  hdr_top <- names(df)
  hdr_bot <- as.character(unlist(df[1, ]))
  body    <- df[-1, , drop = FALSE]

  ## Blank spanner labels are placeholders meaning "merge with the cell
  ## to the left"; flextable's merge_h does that when the values match,
  ## so carry the previous non-blank label forward.
  hdr_top <- trimws(hdr_top)
  for (i in seq_along(hdr_top)) {
    if (i > 1 && (is.na(hdr_top[i]) || hdr_top[i] == "" ||
                  grepl("^V[0-9]+$", hdr_top[i]))) {
      hdr_top[i] <- hdr_top[i - 1]
    }
  }
  hdr_top[is.na(hdr_top)] <- ""
  hdr_bot[is.na(hdr_bot)] <- ""

  ## Unique internal names so set_header_labels can address each column
  names(body) <- paste0("c", seq_len(ncol(body)))
  rownames(body) <- NULL

  ft <- flextable(body) %>%
    set_header_labels(values = setNames(as.list(hdr_bot), names(body))) %>%
    add_header_row(values = hdr_top, colwidths = rep(1, length(hdr_top))) %>%
    merge_h(part = "header") %>%
    bold(part = "header") %>%
    align(align = "center", part = "header") %>%
    align(align = "center", part = "body") %>%
    align(align = "left", j = 1, part = "body") %>%
    valign(valign = "center", part = "all") %>%
    theme_booktabs() %>%
    font(fontname = "Times New Roman", part = "all") %>%
    fontsize(size = fontsize_body, part = "body") %>%
    fontsize(size = fontsize_header, part = "header") %>%
    padding(padding = 1, part = "all")

  ft <- block_structure(ft, body[[1]])
  ft <- attach_notes(ft, notes, hdr_top, hdr_bot, body[[1]])
  fit_columns(ft, hdr_bot)
}
