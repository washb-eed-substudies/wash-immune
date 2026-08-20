## ===============================================================
## EMM supplementary tables - one standalone docx per timepoint.
##
## All table-building logic lives in _emm_helpers.R so that this script
## and manuscript/supplement.qmd cannot drift apart. This file only
## chooses the page setup and writes the files.
##
## For the full rationale on which modifier belongs to which timepoint,
## see the header of _emm_helpers.R and of src/09.
##
## Output: tables/bangladesh-immune-posthoc-EMM-supp-tables-{t2,t3,delta}.docx
## The number of tables per timepoint is driven by which modifiers
## survived the temporality assignment and the sparse-cell screen in
## src/09 - it is no longer a fixed 9.
## ===============================================================

rm(list = ls())

source(here::here("table scripts/_emm_helpers.R"))

dat <- load_emm_results()

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
    body_add_par(emm_intro_text(tp), style = "Normal")

  table_idx <- 0
  for (i in seq_len(nrow(modifier_lut))) {
    V_code   <- modifier_lut$V[i]
    mod_name <- modifier_lut$Modifier_name[i]
    df_mod   <- df_tp %>% filter(V == V_code)
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
      body_add_flextable(ft)

    mod_note <- modifier_lut$Modifier_note[i]
    if (!is.na(mod_note) && nzchar(mod_note)) {
      doc <- doc %>% body_add_par(mod_note, style = "Normal")
    }
    doc <- doc %>% body_add_par("", style = "Normal")
  }

  out_path <- here(sprintf(
    "tables/bangladesh-immune-posthoc-EMM-supp-tables-%s.docx",
    tolower(tp)
  ))
  print(doc, target = out_path)
  message("Wrote ", out_path)
}
