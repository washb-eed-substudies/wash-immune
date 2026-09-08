rm(list=ls())
library("xtable")
source(here::here("0-config.R"))
setwd(paste0(dropboxDir,"Data/Cleaned/Audrie/"))

load(here('results/immune_N_tr_means.RData'))
load(here('results/immune_unadj_glm.RData'))
load(here('results/immune_adj_sex_age_glm.RData'))
load(here('results/immune_adj_glm.RData'))
load(here('results/immune_ipcw.RData'))

## Display rule matching bonpval() in "table scripts/immune-main-tables.R",
## so the sum score table renders p-values like the rest of the paper and
## never prints a small p-value as "0".
fmt_p <- function(p){
  if (p < 0.001)
    return ("<0.001")
  if (p < 0.01)
    return ("<0.01")
  as.character(round(p, 2))
}

## Bonferroni, project convention: double the raw p-value, cap at 1.
bonpval <- function(var) fmt_p(min(1, var[6] * 2))

## The uncorrected p-value. The pre-analysis plan calls for both.
rawpval <- function(var) fmt_p(var[6])

makecival<-function(var){
  rounded<-round(var, 2)
  paste(rounded[1], " (", rounded[2], ", ", rounded[3], ")", sep="")
}

outcome <- c("Outcome, Arm", "Sum score at Year 1", "Control", "Nutrition + WSH",
             "Sum score at Year 2", "Control", "Nutrition + WSH")

N <- c("N", " ", sumscore_t2_N_tr$N[1], sumscore_t2_N_tr$N[2],
       " ", sumscore_t3_N_tr$N[1], sumscore_t3_N_tr$N[2])

mean <- c("Mean", " ", round(sumscore_t2_N_tr$Mean[1], 2), round(sumscore_t2_N_tr$Mean[2], 2),
          " ", round(sumscore_t3_N_tr$Mean[1], 2), round(sumscore_t3_N_tr$Mean[2], 2))

sd <- c("SD", " ", round(sumscore_t2_N_tr$SD[1], 2), round(sumscore_t2_N_tr$SD[2], 2),
        " ", round(sumscore_t3_N_tr$SD[1], 2), round(sumscore_t3_N_tr$SD[2], 2))

unadj <- c("95% CI", " ", " ", makecival(sumscore_t2_unadj_L), " ", " ", makecival(sumscore_t3_unadj_L))

unadjrawpval <- c("P-value", " ", " ", rawpval(sumscore_t2_unadj_L), " "," ", rawpval(sumscore_t3_unadj_L))

unadjpval <- c("Bonferroni P-value", " ", " ", bonpval(sumscore_t2_unadj_L), " "," ", bonpval(sumscore_t3_unadj_L))

agesex <- c("95% CI", " ", " ", makecival(sumscore_t2_adj_sex_age_L), " ", " ", makecival(sumscore_t3_adj_sex_age_L))

agesexrawpval <- c("P-value", " ", " ", rawpval(sumscore_t2_adj_sex_age_L), " ", " ", rawpval(sumscore_t3_adj_sex_age_L))

agesexpval <- c("Bonferroni P-value", " ", " ", bonpval(sumscore_t2_adj_sex_age_L), " ", " ", bonpval(sumscore_t3_adj_sex_age_L))

adj <- c("95% CI", " ", " ", makecival(sumscore_t2_adj_L), " ", " ", makecival(sumscore_t3_adj_L))

adjrawpval <- c("P-value", " ", " ", rawpval(sumscore_t2_adj_L), " "," ", rawpval(sumscore_t3_adj_L))

adjpval <- c("Bonferroni P-value", " ", " ", bonpval(sumscore_t2_adj_L), " "," ", bonpval(sumscore_t3_adj_L))

makeipcwcival<-function(vec){
  rounded<-round(vec, 2)
  paste(rounded[1], " (", rounded[3], ", ", rounded[4], ")", sep="")
}

ipcw <- c("95% CI", " ", " ", makeipcwcival(sumscore_t2_adj_ipcw_L$`unlist(sumscore_t2_adj_ipcw$estimates$ATE`), " ", " ", makeipcwcival(sumscore_t3_adj_ipcw_L$`unlist(sumscore_t3_adj_ipcw$estimates$ATE)`))

ipcwbonpval <- function(pval) fmt_p(min(1, pval * 2))
ipcwrawp <- function(pval) fmt_p(pval)

ipcwrawpval <- c("P-value", " ", " ", ipcwrawp(sumscore_t2_adj_ipcw_L$`unlist(sumscore_t2_adj_ipcw$estimates$ATE)`[5]), " "," ", ipcwrawp(sumscore_t3_adj_ipcw_L$`unlist(sumscore_t3_adj_ipcw$estimates$ATE)`[5]))

ipcwpval <- c("Bonferroni P-value", " ", " ", ipcwbonpval(sumscore_t2_adj_ipcw_L$`unlist(sumscore_t2_adj_ipcw$estimates$ATE)`[5]), " "," ", ipcwbonpval(sumscore_t3_adj_ipcw_L$`unlist(sumscore_t3_adj_ipcw$estimates$ATE)`[5]))


tbl <- data.table(
  " " = outcome,
  " " = N,
  " " = mean, 
  " " = sd,
  "Unadjusted difference: Intervention vs. Control" = unadj,
  " " = unadjrawpval,
  " " = unadjpval,
  "Age-sex-adjusted difference: Intervention vs. Control" = agesex,
  " " = agesexrawpval,
  " " = agesexpval,
  "Fully adjusted difference: Intervention vs. Control" = adj,
  " " = adjrawpval,
  " " = adjpval,
  "IPCW adjusted difference: Intervention vs. Control" = ipcw,
  " " = ipcwrawpval,
  " " = ipcwpval
)

write.csv(tbl, here('tables/supplementary/sumscore_table.csv'))


## Machine-readable copy for the Quarto supplement (manuscript/supplement.qmd).
saveRDS(list(S8 = tbl), here("results/supplement_tables_8.rds"))
