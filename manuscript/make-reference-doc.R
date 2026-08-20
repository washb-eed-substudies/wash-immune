## ===============================================================
## Build the Word reference document used by manuscript/supplement.qmd.
##
## The supplementary tables are wide (up to 14 columns), so the whole
## supplement is set landscape with narrow margins. Quarto/pandoc take
## page setup and paragraph styles from a reference .docx; we start from
## pandoc's own default so every style pandoc emits ("Body Text",
## "Table Caption", "Heading 1", ...) exists, then override the section.
##
## Only needs re-running if the page setup changes. Committing the
## generated file means a plain `quarto render` works with no setup.
##
## Usage:  Rscript manuscript/make-reference-doc.R
## ===============================================================

library(here)
library(officer)
library(magrittr)

base_doc <- here("manuscript/_pandoc-reference-base.docx")
out_doc  <- here("manuscript/reference-landscape.docx")

if (!file.exists(base_doc)) {
  stop("Missing ", base_doc, "\n",
       "Regenerate it with:\n",
       "  pandoc --print-default-data-file reference.docx > ",
       "manuscript/_pandoc-reference-base.docx")
}

read_docx(path = base_doc) %>%
  body_set_default_section(prop_section(
    page_size    = page_size(orient = "landscape", width = 11, height = 8.5),
    page_margins = page_mar(top = 0.6, bottom = 0.6, left = 0.6, right = 0.6,
                            header = 0.3, footer = 0.3, gutter = 0)
  )) %>%
  print(target = out_doc)

message("Wrote ", out_doc)
