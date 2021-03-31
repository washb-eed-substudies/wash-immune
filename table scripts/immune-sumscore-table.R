rm(list=ls())
library("xtable")
source(here::here("0-config.R"))
setwd(paste0(dropboxDir,"Data/Cleaned/Audrie/"))

load(here('results/immune_N_tr_means.RData'))
load(here('results/immune_unadj_glm.RData'))
load(here('results/immune_adj_sex_age_glm.RData'))
load(here('results/immune_adj_glm.RData'))
load(here('results/immune_ipcw.RData'))

bonpval <- function(var){
  bon = round(var[6] * 2, 2)
  if (bon > 1)
    bon = 1
  as.character(bon) 
}

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

unadjpval <- c("Bonferroni P-value", " ", " ", bonpval(sumscore_t2_unadj_L), " "," ", bonpval(sumscore_t3_unadj_L))

agesex <- c("95% CI", " ", " ", makecival(sumscore_t2_adj_sex_age_L), " ", " ", makecival(sumscore_t3_adj_sex_age_L))

agesexpval <- c("Bonferroni P-value", " ", " ", bonpval(sumscore_t2_adj_sex_age_L), " ", " ", bonpval(sumscore_t3_adj_sex_age_L))

adj <- c("95% CI", " ", " ", makecival(sumscore_t2_adj_L), " ", " ", makecival(sumscore_t3_adj_L))

adjpval <- c("Bonferroni P-value", " ", " ", bonpval(sumscore_t2_adj_L), " "," ", bonpval(sumscore_t3_adj_L))

makeipcwcival<-function(vec){
  rounded<-round(vec, 2)
  paste(rounded[1], " (", rounded[3], ", ", rounded[4], ")", sep="")
}

ipcw <- c("95% CI", " ", " ", makeipcwcival(sumscore_t2_adj_ipcw_L$`unlist(sumscore_t2_adj_ipcw$estimates$ATE`), " ", " ", makeipcwcival(sumscore_t3_adj_ipcw_L$`unlist(sumscore_t3_adj_ipcw$estimates$ATE)`))

ipcwbonpval <- function(pval){
  bon = round(pval * 2, 2)
  if (bon > 1)
    bon = 1
  as.character(bon) 
}

ipcwpval <- c("Bonferroni P-value", " ", " ", ipcwbonpval(sumscore_t2_adj_ipcw_L$`unlist(sumscore_t2_adj_ipcw$estimates$ATE)`[5]), " "," ", ipcwbonpval(sumscore_t3_adj_ipcw_L$`unlist(sumscore_t3_adj_ipcw$estimates$ATE)`[5]))


tbl <- data.table(
  " " = outcome,
  " " = N,
  " " = mean, 
  " " = sd,
  "Unadjusted difference: Intervention vs. Control" = unadj,
  " " = unadjpval,
  "Age-sex-adjusted difference: Intervention vs. Control" = agesex,
  " " = agesexpval,
  "Fully adjusted difference: Intervention vs. Control" = adj,
  " " = adjpval,
  "IPCW adjusted difference: Intervention vs. Control" = ipcw,
  " " = ipcwpval
)

write.csv(tbl, here('tables/supplementary/sumscore_table.csv'))
