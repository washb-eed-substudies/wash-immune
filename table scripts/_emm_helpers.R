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
## Interaction p-values are Benjamini-Hochberg FDR-adjusted within each
## (timepoint, modifier) family - i.e. across the cytokine ratios for a
## given modifier. Stratum-specific (per-subgroup) RD p-values are
## presented raw/unadjusted: the hypothesis under test in an EMM analysis
## is the interaction, so the per-subgroup effects are descriptive and
## should be read alongside their 95% CIs.
##
## SPECIMEN TIMING (see the header of src/09 for the full audit):
##   14-month outcomes -> 14-month TAC qPCR modifiers (same-day stool)
##   28-month outcomes -> ~30-month year-2 parasite specimen (+2 months)
##   change score      -> 14-month TAC qPCR modifiers (precede T3)
## The ~30-month specimen post-dates the 14-month blood draw by a median
## of 16.0 months and must not be used as a modifier there.
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

## BH-FDR adjustment over a numeric vector, preserving NAs.
bh_fdr <- function(p) {
  out <- rep(NA_real_, length(p))
  idx <- !is.na(p)
  if (any(idx)) out[idx] <- stats::p.adjust(p[idx], method = "BH")
  out
}

## Format a p-value with the project's standard conventions.
fmt_p_emm <- function(p) {
  ifelse(is.na(p), "", format.pval(p, digits = 2, eps = 0.001))
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
## Cytokine labels and ordering
## ---------------------------------------------------------------

cytokine_lut <- tibble::tribble(
  ~raw,    ~label,
  "il1",   "Interleukin-1β",
  "il2",   "Interleukin-2",
  "il4",   "Interleukin-4",
  "il5",   "Interleukin-5",
  "il6",   "Interleukin-6",
  "il10",  "Interleukin-10",
  "il12",  "Interleukin-12",
  "il13",  "Interleukin-13",
  "il17",  "Interleukin-17",
  "il21",  "Interleukin-21",
  "ifn",   "Interferon-γ",
  "tnf",   "Tumor necrosis factor-α",
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
## cytokine labels attached and interaction p-values FDR-adjusted within
## each (timepoint, modifier) family.
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

  ## BH-FDR for interaction tests within (timepoint, modifier).
  ## Stratum-specific RD p-values are left raw - see header.
  int_p <- dat %>%
    group_by(timepoint, Y, V) %>%
    summarise(int_Pval = first(na.omit(int_Pval)), .groups = "drop") %>%
    group_by(timepoint, V) %>%
    mutate(int_P_FDR = bh_fdr(int_Pval)) %>%
    ungroup() %>%
    select(timepoint, Y, V, int_P_FDR)

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

  rd_ci_1 <- sprintf("%.3f (%.3f, %.3f)", df_y$RD[1], df_y$ci.lb[1], df_y$ci.ub[1])
  rd_ci_2 <- sprintf("%.3f (%.3f, %.3f)", df_y$RD[2], df_y$ci.lb[2], df_y$ci.ub[2])
  p1      <- fmt_p_emm(df_y$P_value[1])   # raw stratum-specific p (not FDR)
  p2      <- fmt_p_emm(df_y$P_value[2])
  int_p   <- fmt_p_emm(first(na.omit(df_y$int_P_FDR)))

  header_row <- tibble::tibble(
    Outcome = df_y$Y_label[1],
    L1_N = "", L1_Mean = "", L1_SD = "", L1_RD = rd_ci_1, L1_P = p1,
    L2_N = "", L2_Mean = "", L2_SD = "", L2_RD = rd_ci_2, L2_P = p2,
    Int_P = int_p
  )

  control_row <- tibble::tibble(
    Outcome = "Control",
    L1_N = fmt_int(df_y$n_Control[1]),
    L1_Mean = fmt_num(df_y$mean_Control[1]),
    L1_SD = fmt_num(df_y$sd_Control[1]),
    L1_RD = "", L1_P = "",
    L2_N = fmt_int(df_y$n_Control[2]),
    L2_Mean = fmt_num(df_y$mean_Control[2]),
    L2_SD = fmt_num(df_y$sd_Control[2]),
    L2_RD = "", L2_P = "",
    Int_P = ""
  )

  nwsh_row <- tibble::tibble(
    Outcome = "Nutrition + WSH",
    L1_N = fmt_int(df_y$n_NWSH[1]),
    L1_Mean = fmt_num(df_y$mean_NWSH[1]),
    L1_SD = fmt_num(df_y$sd_NWSH[1]),
    L1_RD = "", L1_P = "",
    L2_N = fmt_int(df_y$n_NWSH[2]),
    L2_Mean = fmt_num(df_y$mean_NWSH[2]),
    L2_SD = fmt_num(df_y$sd_NWSH[2]),
    L2_RD = "", L2_P = "",
    Int_P = ""
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

  ## Bottom header (column labels) - 5 cols per subgroup
  hdr_bot <- c(
    "Outcome",
    "N", "Mean", "SD", "Difference: Intervention vs. Control (95% CI)", "p",
    "N", "Mean", "SD", "Difference: Intervention vs. Control (95% CI)", "p",
    ""
  )
  ## Top header (group spanners) - each subgroup spans 5 cols
  hdr_top <- c(
    "Outcome / Arm",
    lev1_lab, "", "", "", "",
    lev2_lab, "", "", "", "",
    "Interaction p (FDR)"
  )

  ft %>%
    set_header_labels(values = setNames(as.list(hdr_bot), names(body))) %>%
    add_header_row(values = hdr_top, colwidths = rep(1, length(hdr_top))) %>%
    merge_h(part = "header") %>%
    bold(part = "header") %>%
    bold(j = "Outcome", i = seq(1, nrow(body), by = 3), part = "body") %>%
    align(align = "center", part = "header") %>%
    align(align = "left",   j = "Outcome", part = "body") %>%
    align(align = "right",  j = c("L1_N","L1_Mean","L1_SD",
                                  "L2_N","L2_Mean","L2_SD"), part = "body") %>%
    valign(valign = "center", part = "all") %>%
    theme_booktabs() %>%
    font(fontname = "Times New Roman", part = "all") %>%
    fontsize(size = 7, part = "body") %>%
    fontsize(size = 8, part = "header") %>%
    padding(padding = 1, part = "all") %>%
    set_table_properties(layout = "autofit", width = 1) %>%
    autofit()
}

## Standard preamble paragraph for an EMM table section.
emm_intro_text <- function(tp_label) {
  paste0("Each table corresponds to one effect modifier. ",
         "For each cytokine ratio: the outcome header row shows the ",
         "subgroup-specific difference (Nutrition + WSH vs. Control, ",
         "95% CI) and the interaction p-value; the Control and ",
         "Nutrition + WSH rows give arm-specific N, mean, and SD ",
         "stratified by modifier level. Interaction p-values are ",
         "Benjamini-Hochberg FDR-adjusted within each (", tp_label,
         ", modifier) family (across the cytokine ratios). ",
         "Stratum-specific p-values are raw (unadjusted) and should ",
         "be read alongside their 95% confidence intervals. ",
         "Each pathogen modifier is labelled with the stool specimen, the ",
         "child age it was measured at, and whether infection was measured ",
         "prior to, concurrently with, or after the outcome; modifiers with ",
         "fewer than ", MIN_CELL_REPORTED, " children in any arm-by-status ",
         "cell were not estimable and are not reported. ",
         "All of these modifiers are measured after randomization and are ",
         "plausibly affected by the intervention, so stratifying on them ",
         "does not preserve the randomized comparison; these analyses are ",
         "descriptive and exploratory.")
}

## ---------------------------------------------------------------
## Legacy table conversion
## ---------------------------------------------------------------
## Tables S1-S12 are built by the older scripts as data.tables whose
## COLUMN NAMES are the spanner row and whose FIRST DATA ROW is the
## column-label row. Convert that shape into a two-row-header flextable.

legacy_flextable <- function(tbl, fontsize_body = 7, fontsize_header = 8) {
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

  flextable(body) %>%
    set_header_labels(values = setNames(as.list(hdr_bot), names(body))) %>%
    add_header_row(values = hdr_top, colwidths = rep(1, length(hdr_top))) %>%
    merge_h(part = "header") %>%
    bold(part = "header") %>%
    align(align = "center", part = "header") %>%
    align(align = "left", j = 1, part = "body") %>%
    valign(valign = "center", part = "all") %>%
    theme_booktabs() %>%
    font(fontname = "Times New Roman", part = "all") %>%
    fontsize(size = fontsize_body, part = "body") %>%
    fontsize(size = fontsize_header, part = "header") %>%
    padding(padding = 1, part = "all") %>%
    set_table_properties(layout = "autofit", width = 1) %>%
    autofit()
}
