
# Can we create the following EMM supplementary tables using all cytokine ratios at 28 months 
#(T3 results to ensure pathogen to outcome time ordering, as you mentioned) (listed in Figure 2) by the following subgroups:
#   
#   1. symptoms (diarrhea, fever)
# 2. parasites (giardia, eh, crypto, trichuris, ascaris)
# 
# If the tables are too big, we could break it down further by cytokine groupings in Figure 2: All Th1/Th2 cytokine ratios (blue color in Figure 2), All Pro/IL-10 ratios (orange color in Fig 2), and Th1/Th17 ratios (green color in Figure 2).
# 


## ---------------------------------------------------------------
## 0. Libraries
## ---------------------------------------------------------------
library(here)
library(dplyr)
library(stringr)
library(tidyr)
library(flextable)
library(officer)

## ---------------------------------------------------------------
## 1. Read data
## ---------------------------------------------------------------
full_res <- read.csv(
  here("results/bangladesh-immune-posthoc-subgroup-results-t3.csv"),
  stringsAsFactors = FALSE
)

## ---------------------------------------------------------------
## 2. Cytokine label lookup (for Y)
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

## ---------------------------------------------------------------
## 3. Keep only T3 cytokine ratios (Y begins with t3_)
## ---------------------------------------------------------------
full_res <- full_res %>%
  filter(str_detect(Y, "^t3_"))

## ---------------------------------------------------------------
## 4. Assign cytokine grouping (Figure 2 color groupings)
## ---------------------------------------------------------------
full_res <- full_res %>%
  mutate(
    cytokine_group = case_when(
      Y %in% paste0("t3_", c(
        "ratio_th1_th2", "ratio_il12_il4", "ratio_ifn_il4",
        "ratio_il12_il5", "ratio_ifn_il5",
        "ratio_il12_il13", "ratio_ifn_il13"
      )) ~ "Group 1: Th1/Th2-related ratios",
      
      Y %in% paste0("t3_", c(
        "ratio_pro_il10", "ratio_il1_il10", "ratio_il6_il10",
        "ratio_tnf_il10", "ratio_il2_il10", "ratio_th1_il10",
        "ratio_th2_il10", "ratio_il12_il10", "ratio_ifn_il10",
        "ratio_il4_il10", "ratio_il5_il10", "ratio_il13_il10",
        "ratio_th17_il10", "ratio_il17_il10", "ratio_il21_il10",
        "ratio_gmc_il10"
      )) ~ "Group 2: Pro/IL-10-related ratios",
      
      Y %in% paste0("t3_", c(
        "ratio_th1_th17", "ratio_il12_il17", "ratio_ifn_il17",
        "ratio_il12_il21", "ratio_ifn_il21"
      )) ~ "Group 3: Th1/Th17-related ratios",
      
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(cytokine_group))

## ---------------------------------------------------------------
## 5. Clean Y labels (publication-ready)
## ---------------------------------------------------------------
clean_Y_label <- function(Y) {
  
  ratio <- Y %>%
    str_remove("^t[0-9]_") %>%
    str_remove("^ratio_")
  
  parts <- str_split(ratio, "_", simplify = TRUE)
  
  num_lab <- cytokine_lut$label[match(parts[, 1], cytokine_lut$raw)]
  den_lab <- cytokine_lut$label[match(parts[, 2], cytokine_lut$raw)]
  
  tibble(Y_label = paste0(num_lab, "/", den_lab))
}

Y_clean_df <- full_res %>%
  distinct(Y) %>%
  rowwise() %>%
  bind_cols(clean_Y_label(.$Y)) %>%
  ungroup()

full_res <- full_res %>%
  left_join(Y_clean_df, by = "Y")

## Optional: order ratios to match Figure 2 listing
ratio_order <- paste0("t3_", c(
  "ratio_th1_th2", "ratio_il12_il4", "ratio_ifn_il4", "ratio_il12_il5", "ratio_ifn_il5",
  "ratio_il12_il13", "ratio_ifn_il13",
  "ratio_pro_il10", "ratio_il1_il10", "ratio_il6_il10", "ratio_tnf_il10", "ratio_il2_il10",
  "ratio_th1_il10", "ratio_th2_il10", "ratio_il12_il10", "ratio_ifn_il10", "ratio_il4_il10",
  "ratio_il5_il10", "ratio_il13_il10", "ratio_th17_il10", "ratio_il17_il10", "ratio_il21_il10",
  "ratio_gmc_il10",
  "ratio_th1_th17", "ratio_il12_il17", "ratio_ifn_il17", "ratio_il12_il21", "ratio_ifn_il21"
))

full_res <- full_res %>%
  mutate(Y = factor(Y, levels = ratio_order)) %>%
  arrange(cytokine_group, Y, V)

## ---------------------------------------------------------------
## 6. Outcome subgroup definitions and labels (V)
## ---------------------------------------------------------------

V_outcome_lut <- tibble::tribble(
  ~V,                     ~Outcome,                                     ~Outcome_group,
  "diar7d_t2",            "Diarrhea (past 7 days, T2)",                  "Symptoms",
  "diar7d_t3",            "Diarrhea (past 7 days, T3)",                  "Symptoms",
  "fever7d_t2",           "Fever (past 7 days, T2)",                     "Symptoms",
  "fever7d_t3",           "Fever (past 7 days, T3)",                     "Symptoms",
  "ch_pos_giardia",       "Giardia positive",                            "Parasites",
  "ch_pos_entamoeba",     "Entamoeba positive",                          "Parasites",
  "ch_pos_crypto",        "Cryptosporidium positive",                    "Parasites",
  "ch_qpcr_pos_trichuris","Trichuris positive (qPCR)",                   "Parasites",
  "ch_qpcr_pos_ascaris",  "Ascaris positive (qPCR)",                     "Parasites"
)

full_res <- full_res %>%
  left_join(V_outcome_lut, by = "V")


V_modifier_lut <- tibble::tribble(
  ~V,                     ~Modifier_label,                     ~Time_label,
  "diar7d_t2",            "Diarrhea (past 7 days)",             "14 days",
  "diar7d_t3",            "Diarrhea (past 7 days)",             "28 days",
  "fever7d_t2",           "Fever (past 7 days)",                "14 days",
  "fever7d_t3",           "Fever (past 7 days)",                "28 days",
  "ch_pos_giardia",       "Giardia positive",                   NA,
  "ch_pos_entamoeba",     "Entamoeba positive",                 NA,
  "ch_pos_crypto",        "Cryptosporidium positive",           NA,
  "ch_qpcr_pos_trichuris","Trichuris positive (qPCR)",          NA,
  "ch_qpcr_pos_ascaris",  "Ascaris positive (qPCR)",            NA
)

full_res <- full_res %>%
  left_join(V_modifier_lut, by = "V") %>%
  mutate(
    Modifier = ifelse(
      Outcome_group == "Symptoms",
      paste0(Modifier_label, ", ", Time_label),
      Modifier_label
    )
  )


modifier_order <- c(
  "Diarrhea (past 7 days), 14 days",
  "Diarrhea (past 7 days), 28 days",
  "Fever (past 7 days), 14 days",
  "Fever (past 7 days), 28 days",
  "Giardia positive",
  "Entamoeba positive",
  "Cryptosporidium positive",
  "Trichuris positive (qPCR)",
  "Ascaris positive (qPCR)"
)

full_res <- full_res %>%
  mutate(Modifier = factor(Modifier, levels = modifier_order))


## ---------------------------------------------------------------
## 7. Standardized formatting for RD, CI, p
## ---------------------------------------------------------------
fmt_p <- function(p) format.pval(p, digits = 2, eps = 0.001)

full_res <- full_res %>%
  mutate(
    RD_num   = RD,
    ci_lb_num = ci.lb,
    ci_ub_num = ci.ub,
    RD_CI = sprintf("%.3f (%.3f, %.3f)", RD_num, ci_lb_num, ci_ub_num),
    p_fmt = fmt_p(`P.value`),
    int_p_fmt = fmt_p(`int_Pval`)
  )

## ---------------------------------------------------------------
## 8. Create wide tables (ratios as rows; outcomes as columns)
##    Each outcome gets two columns: RD (CI) and p
## ---------------------------------------------------------------
make_emm_table <- function(df, outcome_group_label, cytokine_group_label) {
  
  df %>%
    filter(
      Outcome_group == outcome_group_label,
      cytokine_group == cytokine_group_label
    ) %>%
    mutate(
      RD_CI = sprintf("%.3f (%.3f, %.3f)", RD, ci.lb, ci.ub),
      p_fmt = format.pval(`P.value`, digits = 2, eps = 0.001),
      int_p_fmt = format.pval(`int_Pval`, digits = 2, eps = 0.001),
      int_p_fmt = gsub("NA","",as.character(int_p_fmt))
    ) %>%
    select(
      `Cytokine ratio` = Y_label,
      Modifier,
      `RD (95% CI)` = RD_CI,
      `p-value` = p_fmt,
      `Int. p-value` = int_p_fmt
    ) %>%
    arrange(`Cytokine ratio`, Modifier)
}

## Table titles (6 total)
cyto_groups <- c(
  "Group 1: Th1/Th2-related ratios",
  "Group 2: Pro/IL-10-related ratios",
  "Group 3: Th1/Th17-related ratios"
)

## Build 6 data.frames
sym_g1 <- make_emm_table(full_res, "Symptoms",  cyto_groups[1])
sym_g2 <- make_emm_table(full_res, "Symptoms",  cyto_groups[2])
sym_g3 <- make_emm_table(full_res, "Symptoms",  cyto_groups[3])

par_g1 <- make_emm_table(full_res, "Parasites", cyto_groups[1])
par_g2 <- make_emm_table(full_res, "Parasites", cyto_groups[2])
par_g3 <- make_emm_table(full_res, "Parasites", cyto_groups[3])

## Convert to flextables
as_ft <- function(df) {
  flextable(df) %>%
    theme_booktabs() %>%
    fontsize(size = 8, part = "all") %>%     # ??? smaller text
    fontsize(size = 9, part = "header") %>%  # slightly larger header
    bold(part = "header") %>%
    align(align = "left", part = "all") %>%
    valign(valign = "center", part = "all") %>%
    padding(padding = 2, part = "all") %>%   # tighter cells
    autofit()
}

ft_sym_g1 <- as_ft(sym_g1)
ft_sym_g2 <- as_ft(sym_g2)
ft_sym_g3 <- as_ft(sym_g3)

ft_par_g1 <- as_ft(par_g1)
ft_par_g2 <- as_ft(par_g2)
ft_par_g3 <- as_ft(par_g3)


## ---------------------------------------------------------------
## 9. Export to Word (single document with 6 tables)
## ---------------------------------------------------------------
doc <- read_docx() %>%
  body_add_par("Supplementary EMM tables (T3 cytokine ratios)", style = "heading 1") %>%
  
  body_add_par("Symptoms outcomes", style = "heading 1") %>%
  body_add_par("Table S1. Group 1: Th1/Th2-related ratios", style = "heading 2") %>%
  body_add_flextable(ft_sym_g1) %>%
  body_add_par("", style = "Normal") %>%
  
  body_add_par("Table S2. Group 2: Pro/IL-10-related ratios", style = "heading 2") %>%
  body_add_flextable(ft_sym_g2) %>%
  body_add_par("", style = "Normal") %>%
  
  body_add_par("Table S3. Group 3: Th1/Th17-related ratios", style = "heading 2") %>%
  body_add_flextable(ft_sym_g3) %>%
  body_add_par("", style = "Normal") %>%
  
  body_add_par("Parasite outcomes", style = "heading 1") %>%
  body_add_par("Table S4. Group 1: Th1/Th2-related ratios", style = "heading 2") %>%
  body_add_flextable(ft_par_g1) %>%
  body_add_par("", style = "Normal") %>%
  
  body_add_par("Table S5. Group 2: Pro/IL-10-related ratios", style = "heading 2") %>%
  body_add_flextable(ft_par_g2) %>%
  body_add_par("", style = "Normal") %>%
  
  body_add_par("Table S6. Group 3: Th1/Th17-related ratios", style = "heading 2") %>%
  body_add_flextable(ft_par_g3)

print(doc, target = here("tables/bangladesh-immune-posthoc-EMM-supp-tables.docx"))
