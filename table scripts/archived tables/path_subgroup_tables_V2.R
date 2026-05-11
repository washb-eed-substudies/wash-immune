## ===============================================================
## Table S12-style pathogen effect modification tables (T3)
## Wide format with spanner headers like Supplementary Table S12
## ===============================================================

rm(list = ls())

library(here)
library(dplyr)
library(stringr)
library(tidyr)
library(flextable)
library(officer)

## ---------------------------------------------------------------
## 1. Read data
## ---------------------------------------------------------------
dat <- read.csv(
  here("results/bangladesh-immune-posthoc-subgroup-results-t3.csv"),
  stringsAsFactors = FALSE
)

## Keep only T3 cytokine ratios (as requested)
dat <- dat %>% filter(str_detect(Y, "^t3_"))

## ---------------------------------------------------------------
## 2. Standardize key column names (handles minor variations)
## ---------------------------------------------------------------

## Required columns: Y, V, RD, ci.lb, ci.ub, p-value, interaction p-value, modifier level
required_base <- c("Y", "V", "RD", "ci.lb", "ci.ub")

missing_base <- setdiff(required_base, names(dat))
if (length(missing_base) > 0) {
  stop(paste0("Missing required columns in CSV: ", paste(missing_base, collapse = ", ")))
}

## P-value column
p_candidates <- c("P.value", "P-value", "p_value", "p.value", "pval", "p_val")
p_col <- p_candidates[p_candidates %in% names(dat)][1]
if (is.na(p_col)) {
  stop("Could not find a p-value column. Expected one of: P.value, P-value, p_value, p.value, pval, p_val.")
}

## Interaction p-value column
int_candidates <- c("int_Pval", "int_pval", "int.pval", "p_int", "pint", "P.int", "Pint")
int_col <- int_candidates[int_candidates %in% names(dat)][1]
if (is.na(int_col)) {
  stop("Could not find an interaction p-value column. Expected one of: int_Pval, int_pval, int.pval, p_int, pint, P.int, Pint.")
}

## ---------------------------------------------------------------
## 3. Cytokine ratio labels (Y -> human readable)
## ---------------------------------------------------------------
cytokine_lut <- tibble::tribble(
  ~raw,    ~label,
  "il1",   "Interleukin-1\u03B2",
  "il2",   "Interleukin-2",
  "il4",   "Interleukin-4",
  "il5",   "Interleukin-5",
  "il6",   "Interleukin-6",
  "il10",  "Interleukin-10",
  "il12",  "Interleukin-12",
  "il13",  "Interleukin-13",
  "il17",  "Interleukin-17",
  "il21",  "Interleukin-21",
  "ifn",   "Interferon-\u03B3",
  "tnf",   "Tumor necrosis factor-\u03B1",
  "gmc",   "Granulocyte-macrophage colony-stimulating factor",
  "th1",   "Th1",
  "th2",   "Th2",
  "th17",  "Th17",
  "pro",   "Pro-inflammatory cytokines"
)

clean_Y_label <- function(Y) {
  ratio <- Y %>%
    str_remove("^t[0-9]_") %>%
    str_remove("^ratio_")
  parts <- str_split(ratio, "_", simplify = TRUE)
  paste0(
    cytokine_lut$label[match(parts[, 1], cytokine_lut$raw)],
    "/",
    cytokine_lut$label[match(parts[, 2], cytokine_lut$raw)]
  )
}

dat <- dat %>% mutate(Y_label = clean_Y_label(Y))

## Order ratios to match Figure 2 / Table S12 ordering
ratio_order <- paste0("t3_", c(
  "ratio_th1_th2","ratio_il12_il4","ratio_ifn_il4","ratio_il12_il5","ratio_ifn_il5",
  "ratio_il12_il13","ratio_ifn_il13",
  "ratio_pro_il10","ratio_il1_il10","ratio_il6_il10","ratio_tnf_il10","ratio_il2_il10",
  "ratio_th1_il10","ratio_th2_il10","ratio_il12_il10","ratio_ifn_il10","ratio_il4_il10",
  "ratio_il5_il10","ratio_il13_il10","ratio_th17_il10","ratio_il17_il10","ratio_il21_il10",
  "ratio_gmc_il10",
  "ratio_th1_th17","ratio_il12_il17","ratio_ifn_il17","ratio_il12_il21","ratio_ifn_il21"
))
dat <- dat %>%
  mutate(Y = factor(Y, levels = ratio_order)) %>%
  arrange(Y)

## ---------------------------------------------------------------
## 4. Pathogen definitions (V codes)
## ---------------------------------------------------------------
pathogen_lut <- tibble::tribble(
  ~V,                      ~Pathogen,
  "ch_pos_giardia",        "Giardia",
  "ch_pos_entamoeba",      "Entamoeba",
  "ch_pos_crypto",         "Cryptosporidium",
  "ch_qpcr_pos_trichuris", "Trichuris (qPCR)",
  "ch_qpcr_pos_ascaris",   "Ascaris (qPCR)"
)

## ---------------------------------------------------------------
## 5. Detect the modifier-level column (two-level strata)
##    If this fails, set modifier_col manually.
## ---------------------------------------------------------------

modifier_col = "subgroup"
message(paste0("Using modifier_col = ", modifier_col))

## ---------------------------------------------------------------
## 6. Formatting helpers
## ---------------------------------------------------------------
fmt_p <- function(p) format.pval(p, digits = 2, eps = 0.001)
fmt_rdci <- function(est, lo, hi) sprintf("%.3f (%.3f, %.3f)", est, lo, hi)

## ---------------------------------------------------------------
## 7. Build a Table S12-style wide table for one pathogen
##    Two strata levels become the two big column blocks (like Female/Male)
## ---------------------------------------------------------------
build_pathogen_table_S12_style <- function(df_one_pathogen) {
  
  ## Ensure subgroup is character
  df_one_pathogen <- df_one_pathogen %>%
    mutate(subgroup = as.character(subgroup))
  
  levs <- sort(unique(df_one_pathogen$subgroup))
  if (length(levs) != 2) {
    stop("Expected exactly 2 subgroup levels, found: ", paste(levs, collapse = ", "))
  }
  
  ## Labels for subgroup levels
  subgroup_label <- function(x) {
    if (x == "0") return("Negative")
    if (x == "1") return("Positive")
    x
  }
  
  lev1 <- levs[1]; lev2 <- levs[2]
  lev1_lab <- subgroup_label(lev1)
  lev2_lab <- subgroup_label(lev2)
  
  ## Format estimates
  df_fmt <- df_one_pathogen %>%
    mutate(
      RD_CI = sprintf("%.3f (%.3f, %.3f)", RD, ci.lb, ci.ub),
      p_fmt = format.pval(P.value, digits = 2, eps = 0.001)
    )
  
  ## One interaction p-value per cytokine ratio
  int_p <- df_fmt %>%
    group_by(Y, Y_label) %>%
    summarise(
      int_p = first(na.omit(int_Pval)),
      .groups = "drop"
    ) %>%
    mutate(int_p = format.pval(int_p, digits = 2, eps = 0.001))
  
  ## Build subgroup blocks aligned by Y / Y_label
  block <- function(level, prefix) {
    df_fmt %>%
      filter(subgroup == level) %>%
      select(Y, Y_label, RD_CI, p_fmt) %>%
      arrange(Y) %>%
      transmute(
        Y,
        Outcome = Y_label,
        !!paste0(prefix, "_N") := "",
        !!paste0(prefix, "_Mean") := "",
        !!paste0(prefix, "_SD") := "",
        !!paste0(prefix, "_RD") := RD_CI,
        !!paste0(prefix, "_P") := p_fmt
      )
  }
  
  b1 <- block(lev1, "L1")
  b2 <- block(lev2, "L2")
  
  ## Join blocks by cytokine ratio (NOT cross join)
  out <- b1 %>%
    left_join(
      b2 %>% select(-Outcome),
      by = "Y"
    ) %>%
    left_join(int_p, by = c("Y", "Outcome" = "Y_label")) %>%
    select(
      Outcome,
      starts_with("L1"),
      starts_with("L2"),
      `P-value for interaction` = int_p
    )
  
  ## Build Table S12-style flextable
  ft <- flextable(out)
  
  hdr2 <- c(
    "Outcome",
    "N","Mean","SD","Unadjusted difference: Intervention vs. Control (95% CI)","P-value",
    "N","Mean","SD","Unadjusted difference: Intervention vs. Control (95% CI)","P-value",
    ""
  )
  
  hdr1 <- c(
    "",
    lev1_lab, rep("", 4),
    lev2_lab, rep("", 4),
    "P-value for interactiona"
  )
  
  ft %>%
    add_header_row(values = hdr2, colwidths = rep(1, length(hdr2))) %>%
    add_header_row(values = hdr1, colwidths = rep(1, length(hdr1))) %>%
    merge_h(part = "header") %>%
    bold(part = "header") %>%
    align(align = "center", part = "header") %>%
    align(align = "left", part = "body") %>%
    valign(valign = "center", part = "all") %>%
    theme_booktabs() %>%
    font(fontname = "Times New Roman", part = "all") %>%
    fontsize(size = 7, part = "body") %>%     # body text
    fontsize(size = 8, part = "header") %>%   # header text
    padding(padding = 1, part = "all") %>%    # tighter rows
    set_table_properties(
      layout = "autofit",
      width = 1
    ) %>%
    autofit()
}


## ---------------------------------------------------------------
## 8. Export all pathogen tables to Word (landscape), S12-like
## ---------------------------------------------------------------
doc <- read_docx() %>%
  body_set_default_section(
    prop_section(
      page_size = page_size(orient = "landscape", width = 11, height = 8.5)
    )
  ) %>%
  body_add_par(
    "Supplementary tables: Effect modification by child pathogen status on cytokine ratios at age 28 months (T3)",
    style = "heading 1"
  ) %>%
  body_add_par(
    "Notes: Confidence intervals were adjusted for clustered observations using robust standard errors. P-values shown are Bonferroni corrected to control for familywise error rates.",
    style = "Normal"
  )

# for (i in seq_len(nrow(pathogen_lut))) {
#   
#   V_code <- pathogen_lut$V[i]
#   pathogen_name <- pathogen_lut$Pathogen[i]
#   
#   df_path <- dat %>% filter(V == V_code)
#   
#   if (nrow(df_path) == 0) {
#     warning(paste0("No rows found for V = ", V_code, " (", pathogen_name, "). Skipping."))
#     next
#   }
#   
#   res <- build_pathogen_table_S12_style(df_path, pathogen_name)
#   
#   doc <- doc %>%
#     body_add_par(
#       paste0("Table S12-", i, ". Effect modification by ", pathogen_name,
#              " status on cytokine ratios at age 28 months."),
#       style = "heading 2"
#     ) %>%
#     body_add_flextable(res$ft) %>%
#     body_add_par("aP-values shown are Bonferroni corrected to control for familywise error rates.", style = "Normal") %>%
#     body_add_par("", style = "Normal")
# }


ft_giardia <- build_pathogen_table_S12_style(dat %>% filter(V == "ch_pos_giardia"))


doc <- read_docx() %>%
  body_add_par(
    "Table S12-X. Effect modification by Giardia status on cytokine ratios at age 28 months.",
    style = "heading 2"
  ) %>%
  body_add_flextable(ft_giardia)

print(doc, target = here("tables/Table_S12_pathogen_example.docx"))

print(doc, target = here("tables/bangladesh-immune-pathogen-EMM-supp-tables-S12style.docx"))
