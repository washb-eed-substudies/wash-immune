## ===============================================================
## EMM supplementary tables — one table per effect modifier.
##
## For each (timepoint, modifier) we build one wide flextable with
## a 3-row block per cytokine ratio:
##   1. Outcome header row carrying subgroup-specific treatment
##      effects (RD, 95% CI, p-value) and the interaction p-value.
##   2. Control row with N, mean, SD broken out by subgroup level.
##   3. Nutrition+WSH row with N, mean, SD broken out by subgroup level.
##
## P-values are Benjamini-Hochberg FDR-adjusted within each
## (timepoint, ratios) family, separately for stratum-specific
## RD p-values and for interaction tests.
##
## Output: one docx per timepoint, each containing 9 tables
## (one per modifier).
## ===============================================================

rm(list = ls())

library(here)
library(dplyr)
library(stringr)
library(tidyr)
library(flextable)
library(officer)

source(here("table scripts/_emm_helpers.R"))

## ---------------------------------------------------------------
## 1. Read all three timepoints
## ---------------------------------------------------------------
result_files <- c(
  T2    = "results/bangladesh-immune-posthoc-subgroup-results-t2.csv",
  T3    = "results/bangladesh-immune-posthoc-subgroup-results-t3.csv",
  Delta = "results/bangladesh-immune-posthoc-subgroup-results-delta.csv"
)

read_one <- function(path, tp) {
  fp <- here(path)
  if (!file.exists(fp)) {
    warning("Missing ", path, " — skipping ", tp)
    return(NULL)
  }
  read.csv(fp, stringsAsFactors = FALSE) %>% mutate(timepoint = tp)
}

dat <- bind_rows(lapply(names(result_files),
                        function(tp) read_one(result_files[[tp]], tp)))

if (is.null(dat) || nrow(dat) == 0) {
  stop("No EMM result CSVs found. Run src/09 first.")
}

## Normalize p-value column names that R turns into "P.value"
if ("P.value" %in% names(dat)) names(dat)[names(dat) == "P.value"] <- "P_value"

## Strip timepoint prefix so the same ratio label applies across timepoints
dat <- dat %>%
  mutate(ratio_stem = Y %>% str_remove("^t[0-9]_") %>% str_remove("^d23_"))

## ---------------------------------------------------------------
## 2. Cytokine labels and ordering
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

dat <- dat %>% mutate(Y_label = clean_ratio_label(ratio_stem))

ratio_order <- c(
  "ratio_th1_th2","ratio_il12_il4","ratio_ifn_il4","ratio_il12_il5","ratio_ifn_il5",
  "ratio_il12_il13","ratio_ifn_il13",
  "ratio_pro_il10","ratio_il1_il10","ratio_il6_il10","ratio_tnf_il10","ratio_il2_il10",
  "ratio_th1_il10","ratio_th2_il10","ratio_il12_il10","ratio_ifn_il10","ratio_il4_il10",
  "ratio_il5_il10","ratio_il13_il10","ratio_th17_il10","ratio_il17_il10","ratio_il21_il10",
  "ratio_gmc_il10",
  "ratio_th1_th17","ratio_il12_il17","ratio_ifn_il17","ratio_il12_il21","ratio_ifn_il21"
)
dat <- dat %>%
  mutate(ratio_stem = factor(ratio_stem, levels = ratio_order)) %>%
  arrange(timepoint, V, ratio_stem, subgroup)

## ---------------------------------------------------------------
## 3. BH-FDR within (timepoint, ratios)
## ---------------------------------------------------------------
dat <- dat %>%
  group_by(timepoint) %>%
  mutate(P_FDR = bh_fdr(P_value)) %>%
  ungroup()

int_p <- dat %>%
  group_by(timepoint, Y, V) %>%
  summarise(int_Pval = first(na.omit(int_Pval)), .groups = "drop") %>%
  group_by(timepoint) %>%
  mutate(int_P_FDR = bh_fdr(int_Pval)) %>%
  ungroup() %>%
  select(timepoint, Y, V, int_P_FDR)

dat <- dat %>% left_join(int_p, by = c("timepoint", "Y", "V"))

## ---------------------------------------------------------------
## 4. Modifier metadata
## ---------------------------------------------------------------
modifier_lut <- tibble::tribble(
  ~V,                      ~Modifier_name,
  "ch_pos_giardia",        "Giardia",
  "ch_pos_entamoeba",      "Entamoeba",
  "ch_pos_crypto",         "Cryptosporidium",
  "ch_qpcr_pos_trichuris", "Trichuris (qPCR)",
  "ch_qpcr_pos_ascaris",   "Ascaris (qPCR)",
  "diar7d_t2",             "Diarrhea (past 7 days at T2)",
  "diar7d_t3",             "Diarrhea (past 7 days at T3)",
  "fever7d_t2",            "Fever (past 7 days at T2)",
  "fever7d_t3",            "Fever (past 7 days at T3)"
)

fmt_num <- function(x) ifelse(is.na(x), "", sprintf("%.2f", x))
fmt_int <- function(x) ifelse(is.na(x), "", as.character(as.integer(x)))

## ---------------------------------------------------------------
## 5. Build one 3-row block per cytokine ratio
## ---------------------------------------------------------------
build_outcome_block <- function(df_y) {
  ## df_y is expected to have exactly 2 rows, one per subgroup level
  df_y <- df_y %>% arrange(subgroup)
  if (nrow(df_y) != 2) return(NULL)

  V_code <- unique(df_y$V)
  lev_raw <- as.character(df_y$subgroup)
  ## lev1 = first level (typically "0"), lev2 = second level (typically "1")

  rd_ci_1 <- sprintf("%.3f (%.3f, %.3f)", df_y$RD[1], df_y$ci.lb[1], df_y$ci.ub[1])
  rd_ci_2 <- sprintf("%.3f (%.3f, %.3f)", df_y$RD[2], df_y$ci.lb[2], df_y$ci.ub[2])
  p1      <- fmt_p_emm(df_y$P_FDR[1])
  p2      <- fmt_p_emm(df_y$P_FDR[2])
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

## ---------------------------------------------------------------
## 6. Assemble one flextable for a (timepoint, modifier) pair
## ---------------------------------------------------------------
build_modifier_table <- function(df_mod) {
  V_code <- unique(df_mod$V)

  ## Determine subgroup labels in the order they appear after arrange()
  lev_raw <- df_mod %>%
    arrange(ratio_stem, subgroup) %>%
    pull(subgroup) %>%
    as.character() %>%
    unique()
  if (length(lev_raw) != 2) {
    warning("Modifier ", V_code, " has ", length(lev_raw),
            " subgroup levels — skipping.")
    return(NULL)
  }
  lev1_lab <- modifier_level_label(V_code, lev_raw[1])
  lev2_lab <- modifier_level_label(V_code, lev_raw[2])

  body <- df_mod %>%
    group_split(ratio_stem) %>%
    lapply(build_outcome_block) %>%
    bind_rows()

  if (nrow(body) == 0) return(NULL)

  ## Build flextable with two-row spanner header
  ft <- flextable(body)

  ## Bottom header (column labels) — now 5 cols per subgroup
  hdr_bot <- c(
    "Outcome",
    "N", "Mean", "SD", "Difference: Intervention vs. Control (95% CI)", "p (FDR)",
    "N", "Mean", "SD", "Difference: Intervention vs. Control (95% CI)", "p (FDR)",
    ""
  )
  ## Top header (group spanners) — each subgroup spans 5 cols
  hdr_top <- c(
    "Outcome / Arm",
    lev1_lab, "", "", "", "",
    lev2_lab, "", "", "", "",
    "Interaction p (FDR)"
  )

  ft <- ft %>%
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

  ft
}

## ---------------------------------------------------------------
## 7. One docx per timepoint
## ---------------------------------------------------------------
for (tp in unique(dat$timepoint)) {
  df_tp <- dat %>% filter(timepoint == tp)
  if (nrow(df_tp) == 0) next

  doc <- read_docx() %>%
    body_set_default_section(prop_section(
      page_size = page_size(orient = "landscape", width = 11, height = 8.5)
    )) %>%
    body_add_par(
      paste0("Supplementary tables: Effect modification of intervention effects on cytokine ratios (", tp, ")"),
      style = "heading 1"
    ) %>%
    body_add_par(
      paste0("Each table corresponds to one effect modifier. ",
             "For each cytokine ratio: the outcome header row shows the ",
             "subgroup-specific difference (Nutrition + WSH vs. Control, ",
             "95% CI) and the interaction p-value; the Control and ",
             "Nutrition + WSH rows give arm-specific N, mean, and SD ",
             "stratified by modifier level. P-values are Benjamini-Hochberg ",
             "FDR-adjusted within the (", tp, ", ratios) family, separately ",
             "for stratum-specific risk differences and interaction tests."),
      style = "Normal"
    )

  table_idx <- 0
  for (i in seq_len(nrow(modifier_lut))) {
    V_code  <- modifier_lut$V[i]
    mod_name <- modifier_lut$Modifier_name[i]
    df_mod <- df_tp %>% filter(V == V_code)
    if (nrow(df_mod) == 0) next

    ft <- build_modifier_table(df_mod)
    if (is.null(ft)) next

    table_idx <- table_idx + 1
    doc <- doc %>%
      body_add_par(
        paste0("Table ", tp, "-", table_idx,
               ". Effect modification by ", mod_name, "."),
        style = "heading 2"
      ) %>%
      body_add_flextable(ft) %>%
      body_add_par("", style = "Normal")
  }

  out_path <- here(sprintf(
    "tables/bangladesh-immune-posthoc-EMM-supp-tables-%s.docx",
    tolower(tp)
  ))
  print(doc, target = out_path)
  message("Wrote ", out_path)
}
