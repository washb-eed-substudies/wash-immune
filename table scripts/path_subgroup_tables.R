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
    )

  plan <- emm_table_plan(df_tp, tp)

  table_idx <- 0
  for (k in seq_len(nrow(plan))) {
    V_code   <- plan$V[k]
    mod_name <- modifier_heading(V_code, tp)
    df_mod   <- df_tp %>% filter(V == V_code)
    if (nrow(df_mod) == 0) next

    ft <- build_modifier_table(df_mod)
    if (is.null(ft)) next

    if (plan$numbered[k]) {
      table_idx <- table_idx + 1
      tbl_label <- paste0(tp, "-", table_idx)
    } else {
      tbl_label <- "SX"
    }

    doc <- doc %>%
      body_add_par(
        paste0("Table ", tbl_label,
               ". Effect modification by ", mod_name, "."),
        style = "heading 2"
      ) %>%
      body_add_flextable(ft)

    ## Footnotes are paragraphs under the table, at 8pt. The specimen
    ## note is now part of that block rather than a separate paragraph.
    for (p in note_fpars(attr(ft, "table_notes")))
      doc <- doc %>% body_add_fpar(p)
    doc <- doc %>% body_add_par("", style = "Normal")
  }

  out_path <- here(sprintf(
    "tables/bangladesh-immune-posthoc-EMM-supp-tables-%s.docx",
    tolower(tp)
  ))
  print(doc, target = out_path)
  message("Wrote ", out_path)
}
