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

## officer has no API for defining a paragraph style, so "Table Note" is
## injected into styles.xml directly. The supplement applies it to the
## footnote paragraphs under each table. Keeping it a real Word style,
## rather than hard-coded run formatting, means a coauthor can restyle
## every note in the document at once.
TABLE_NOTE_STYLE <- paste0(
  '<w:style w:type="paragraph" w:customStyle="1" w:styleId="TableNote">',
  '<w:name w:val="Table Note"/>',
  '<w:basedOn w:val="Normal"/>',
  '<w:qFormat/>',
  '<w:pPr><w:spacing w:before="0" w:after="40"/>',
  '<w:contextualSpacing w:val="true"/></w:pPr>',
  '<w:rPr><w:rFonts w:ascii="Times New Roman" w:hAnsi="Times New Roman"',
  ' w:cs="Times New Roman"/><w:sz w:val="16"/><w:szCs w:val="16"/></w:rPr>',
  '</w:style>')

add_table_note_style <- function(path) {
  tmp <- file.path(tempdir(), "refdoc")
  unlink(tmp, recursive = TRUE)
  dir.create(tmp, recursive = TRUE)
  utils::unzip(path, exdir = tmp)

  sty <- file.path(tmp, "word", "styles.xml")
  x <- paste(readLines(sty, warn = FALSE, encoding = "UTF-8"), collapse = "")
  if (!grepl('w:styleId="TableNote"', x, fixed = TRUE)) {
    x <- sub("</w:styles>", paste0(TABLE_NOTE_STYLE, "</w:styles>"), x,
             fixed = TRUE)
    con <- file(sty, open = "wb")
    writeBin(charToRaw(x), con)
    close(con)
  }
  zip::zipr(path, list.files(tmp, all.files = TRUE, no.. = TRUE,
                             full.names = TRUE))
  invisible(path)
}

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

add_table_note_style(out_doc)

message("Wrote ", out_doc)
