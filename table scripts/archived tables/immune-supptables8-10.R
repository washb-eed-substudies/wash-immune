rm(list=ls())
library("xtable")
source(here::here("0-config.R"))

load(here("results/immune_unadj_glm.RData"))
load(here("results/immune_adj_glm.RData"))

bonpval <- function(pval){
  bon = round(pval * 2, 2)
  if (pval >= .5)
    bon = 1
  bon 
}

#### TABLE S8 ####
outcomes10<-c("Outcome", "Ln IL-1\u03b2 (pg/ml)", "Ln Il-6 (pg/ml)", "Ln TNF-\u03b1 (pg/ml)", "Ln CRP (mg/L)", "Ln IL-12 (pg/ml)",
             "Ln IFN-\u03b3 (pg/ml)", "Ln IL-4 (pg/ml)", "Ln IL-5 (pg/ml)", "Ln IL-13 (pg/ml)", "Ln IL-17A (pg/ml)", 
             "Ln IL-21 (pg/ml)", "Ln IL-10 (pg/ml)", "Ln IL-2 (pg/ml)", "Ln GM-CSF (pg/ml)", "Ln AGP (g/L)", "Ln IGF-1 (\u03bcg/L)",
             "Ln IL-1\u03b2/IL-10", "Ln IL-6/IL-10", "Ln TNF-\u03b1/IL-10", "Ln IL-12/IL-10", "Ln IFN-\u03b3/IL-10",
             "Ln IL-4/IL-10", "Ln IL-5/IL-10", "Ln IL-13/IL-10", "Ln IL-17A/IL-10", "Ln IL-21/IL-10",
             "Ln IL-2/IL-10", "Ln GM-CSF/IL-10", "Ln IL-12/IL-4", "Ln IFN-\u03b3/IL-4", "Ln IL-12/IL-5", "Ln IFN-\u03b3/IL-5",
             "Ln IL-12/IL-13", "Ln IFN-\u03b3/IL-13", "Ln IL-12/IL-17A", "Ln IFN-\u03b3/IL-17A", "Ln IL-12/IL-21", "Ln IFN-\u03b3/IL-21",
             "Ln Pro-inflammatory cytokines/IL-10", "Ln Th1/IL-10", "Ln Th2/IL-10", "Ln Th17/IL-10",
             "Ln Th1/Th2", "Ln Th1/Th17")

unadjdiffs10 <- c("Unadjusted Difference", 
                  as.character(round(t2_il1_unadj_L[1], 2)), as.character(round(t2_il6_unadj_L[1], 2)), as.character(round(t2_tnf_unadj_L[1], 2)),
                  as.character(round(t2_crp_unadj_L[1], 2)), as.character(round(t2_il12_unadj_L[1], 2)), as.character(round(t2_ifn_unadj_L[1], 2)),
                  as.character(round(t2_il4_unadj_L[1], 2)), as.character(round(t2_il5_unadj_L[1], 2)), as.character(round(t2_il13_unadj_L[1], 2)),
                  as.character(round(t2_il17_unadj_L[1], 2)), as.character(round(t2_il21_unadj_L[1], 2)), as.character(round(t2_il10_unadj_L[1], 2)),
                  as.character(round(t2_il2_unadj_L[1], 2)), as.character(round(t2_gmc_unadj_L[1], 2)), as.character(round(t2_agp_unadj_L[1], 2)),
                  as.character(round(t2_igf_unadj_L[1], 2)), as.character(round(t2_ratio_il1_il10_unadj_L[1], 2)),
                  as.character(round(t2_ratio_il6_il10_unadj_L[1], 2)), as.character(round(t2_ratio_tnf_il10_unadj_L[1], 2)), as.character(round(t2_ratio_il12_il10_unadj_L[1], 2)),
                  as.character(round(t2_ratio_ifn_il10_unadj_L[1], 2)), as.character(round(t2_ratio_il4_il10_unadj_L[1], 2)), as.character(round(t2_ratio_il5_il10_unadj_L[1], 2)),
                  as.character(round(t2_ratio_il13_il10_unadj_L[1], 2)), as.character(round(t2_ratio_il17_il10_unadj_L[1], 2)), as.character(round(t2_ratio_il21_il10_unadj_L[1], 2)),
                  as.character(round(t2_ratio_il2_il10_unadj_L[1], 2)), as.character(round(t2_ratio_gmc_il10_unadj_L[1], 2)), as.character(round(t2_ratio_il12_il4_unadj_L[1], 2)),
                  as.character(round(t2_ratio_ifn_il4_unadj_L[1], 2)), as.character(round(t2_ratio_il12_il5_unadj_L[1], 2)), as.character(round(t2_ratio_ifn_il5_unadj_L[1], 2)),
                  as.character(round(t2_ratio_il12_il13_unadj_L[1], 2)), as.character(round(t2_ratio_ifn_il13_unadj_L[1], 2)), as.character(round(t2_ratio_il12_il17_unadj_L[1], 2)),
                  as.character(round(t2_ratio_ifn_il17_unadj_L[1], 2)), as.character(round(t2_ratio_il12_il21_unadj_L[1], 2)), as.character(round(t2_ratio_ifn_il21_unadj_L[1], 2)),
                  as.character(round(t2_ratio_pro_il10_unadj_L[1], 2)), as.character(round(t2_ratio_th1_il10_unadj_L[1], 2)), as.character(round(t2_ratio_th2_il10_unadj_L[1], 2)),
                  as.character(round(t2_ratio_th17_il10_unadj_L[1], 2)), as.character(round(t2_ratio_th1_th2_unadj_L[1], 2)), as.character(round(t2_ratio_th1_th17_unadj_L[1], 2))
                  )

unadjpvals10 <- c("Unadjusted P-value", 
                  as.character(round(t2_il1_unadj_L[6], 2)), as.character(round(t2_il6_unadj_L[6], 2)), as.character(round(t2_tnf_unadj_L[6], 2)),
                  as.character(round(t2_crp_unadj_L[6], 2)), as.character(round(t2_il12_unadj_L[6], 2)), as.character(round(t2_ifn_unadj_L[6], 2)),
                  as.character(round(t2_il4_unadj_L[6], 2)), as.character(round(t2_il5_unadj_L[6], 2)), as.character(round(t2_il13_unadj_L[6], 2)),
                  as.character(round(t2_il17_unadj_L[6], 2)), as.character(round(t2_il21_unadj_L[6], 2)), as.character(round(t2_il10_unadj_L[6], 2)),
                  as.character(round(t2_il2_unadj_L[6], 2)), as.character(round(t2_gmc_unadj_L[6], 2)), as.character(round(t2_agp_unadj_L[6], 2)),
                  as.character(round(t2_igf_unadj_L[6], 2)), as.character(round(t2_ratio_il1_il10_unadj_L[6], 2)),
                  as.character(round(t2_ratio_il6_il10_unadj_L[6], 2)), as.character(round(t2_ratio_tnf_il10_unadj_L[6], 2)), as.character(round(t2_ratio_il12_il10_unadj_L[6], 2)),
                  as.character(round(t2_ratio_ifn_il10_unadj_L[6], 2)), as.character(round(t2_ratio_il4_il10_unadj_L[6], 2)), as.character(round(t2_ratio_il5_il10_unadj_L[6], 2)),
                  as.character(round(t2_ratio_il13_il10_unadj_L[6], 2)), as.character(round(t2_ratio_il17_il10_unadj_L[6], 2)), as.character(round(t2_ratio_il21_il10_unadj_L[6], 2)),
                  as.character(round(t2_ratio_il2_il10_unadj_L[6], 2)), as.character(round(t2_ratio_gmc_il10_unadj_L[6], 2)), as.character(round(t2_ratio_il12_il4_unadj_L[6], 2)),
                  as.character(round(t2_ratio_ifn_il4_unadj_L[6], 2)), as.character(round(t2_ratio_il12_il5_unadj_L[6], 2)), as.character(round(t2_ratio_ifn_il5_unadj_L[6], 2)),
                  as.character(round(t2_ratio_il12_il13_unadj_L[6], 2)), as.character(round(t2_ratio_ifn_il13_unadj_L[6], 2)), as.character(round(t2_ratio_il12_il17_unadj_L[6], 2)),
                  as.character(round(t2_ratio_ifn_il17_unadj_L[6], 2)), as.character(round(t2_ratio_il12_il21_unadj_L[6], 2)), as.character(round(t2_ratio_ifn_il21_unadj_L[6], 2)),
                  as.character(round(t2_ratio_pro_il10_unadj_L[6], 2)), as.character(round(t2_ratio_th1_il10_unadj_L[6], 2)), as.character(round(t2_ratio_th2_il10_unadj_L[6], 2)),
                  as.character(round(t2_ratio_th17_il10_unadj_L[6], 2)), as.character(round(t2_ratio_th1_th2_unadj_L[6], 2)), as.character(round(t2_ratio_th1_th17_unadj_L[6], 2)))

unadjbonpvals10 <- c("Bonferroni P-value", 
                     as.character(bonpval(t2_il1_unadj_L[6])), as.character(bonpval(t2_il6_unadj_L[6])), as.character(bonpval(t2_tnf_unadj_L[6])),
                     as.character(bonpval(t2_crp_unadj_L[6])), as.character(bonpval(t2_il12_unadj_L[6])), as.character(bonpval(t2_ifn_unadj_L[6])),
                     as.character(bonpval(t2_il4_unadj_L[6])), as.character(bonpval(t2_il5_unadj_L[6])), as.character(bonpval(t2_il13_unadj_L[6])),
                     as.character(bonpval(t2_il17_unadj_L[6])), as.character(bonpval(t2_il21_unadj_L[6])), as.character(bonpval(t2_il10_unadj_L[6])),
                     as.character(bonpval(t2_il2_unadj_L[6])), as.character(bonpval(t2_gmc_unadj_L[6])), as.character(bonpval(t2_agp_unadj_L[6])),
                     as.character(bonpval(t2_igf_unadj_L[6])), as.character(bonpval(t2_ratio_il1_il10_unadj_L[6])),
                     as.character(bonpval(t2_ratio_il6_il10_unadj_L[6])), as.character(bonpval(t2_ratio_tnf_il10_unadj_L[6])), as.character(bonpval(t2_ratio_il12_il10_unadj_L[6])),
                     as.character(bonpval(t2_ratio_ifn_il10_unadj_L[6])), as.character(bonpval(t2_ratio_il4_il10_unadj_L[6])), as.character(bonpval(t2_ratio_il5_il10_unadj_L[6])),
                     as.character(bonpval(t2_ratio_il13_il10_unadj_L[6])), as.character(bonpval(t2_ratio_il17_il10_unadj_L[6])), as.character(bonpval(t2_ratio_il21_il10_unadj_L[6])),
                     as.character(bonpval(t2_ratio_il2_il10_unadj_L[6])), as.character(bonpval(t2_ratio_gmc_il10_unadj_L[6])), as.character(bonpval(t2_ratio_il12_il4_unadj_L[6])),
                     as.character(bonpval(t2_ratio_ifn_il4_unadj_L[6])), as.character(bonpval(t2_ratio_il12_il5_unadj_L[6])), as.character(bonpval(t2_ratio_ifn_il5_unadj_L[6])),
                     as.character(bonpval(t2_ratio_il12_il13_unadj_L[6])), as.character(bonpval(t2_ratio_ifn_il13_unadj_L[6])), as.character(bonpval(t2_ratio_il12_il17_unadj_L[6])),
                     as.character(bonpval(t2_ratio_ifn_il17_unadj_L[6])), as.character(bonpval(t2_ratio_il12_il21_unadj_L[6])), as.character(bonpval(t2_ratio_ifn_il21_unadj_L[6])),
                     as.character(bonpval(t2_ratio_pro_il10_unadj_L[6])), as.character(bonpval(t2_ratio_th1_il10_unadj_L[6])), as.character(bonpval(t2_ratio_th2_il10_unadj_L[6])),
                     as.character(bonpval(t2_ratio_th17_il10_unadj_L[6])), as.character(bonpval(t2_ratio_th1_th2_unadj_L[6])), as.character(bonpval(t2_ratio_th1_th17_unadj_L[6]))
                     )

adjdiffs10 <- c("Adjusted Difference",  
                as.character(round(t2_il1_adj_L[1], 2)), as.character(round(t2_il6_adj_L[1], 2)), as.character(round(t2_tnf_adj_L[1], 2)),
                as.character(round(t2_crp_adj_L[1], 2)), as.character(round(t2_il12_adj_L[1], 2)), as.character(round(t2_ifn_adj_L[1], 2)),
                as.character(round(t2_il4_adj_L[1], 2)), as.character(round(t2_il5_adj_L[1], 2)), as.character(round(t2_il13_adj_L[1], 2)),
                as.character(round(t2_il17_adj_L[1], 2)), as.character(round(t2_il21_adj_L[1], 2)), as.character(round(t2_il10_adj_L[1], 2)),
                as.character(round(t2_il2_adj_L[1], 2)), as.character(round(t2_gmc_adj_L[1], 2)), as.character(round(t2_agp_adj_L[1], 2)),
                as.character(round(t2_igf_adj_L[1], 2)), as.character(round(t2_ratio_il1_il10_adj_L[1], 2)),
                as.character(round(t2_ratio_il6_il10_adj_L[1], 2)), as.character(round(t2_ratio_tnf_il10_adj_L[1], 2)), as.character(round(t2_ratio_il12_il10_adj_L[1], 2)),
                as.character(round(t2_ratio_ifn_il10_adj_L[1], 2)), as.character(round(t2_ratio_il4_il10_adj_L[1], 2)), as.character(round(t2_ratio_il5_il10_adj_L[1], 2)),
                as.character(round(t2_ratio_il13_il10_adj_L[1], 2)), as.character(round(t2_ratio_il17_il10_adj_L[1], 2)), as.character(round(t2_ratio_il21_il10_adj_L[1], 2)),
                as.character(round(t2_ratio_il2_il10_adj_L[1], 2)), as.character(round(t2_ratio_gmc_il10_adj_L[1], 2)), as.character(round(t2_ratio_il12_il4_adj_L[1], 2)),
                as.character(round(t2_ratio_ifn_il4_adj_L[1], 2)), as.character(round(t2_ratio_il12_il5_adj_L[1], 2)), as.character(round(t2_ratio_ifn_il5_adj_L[1], 2)),
                as.character(round(t2_ratio_il12_il13_adj_L[1], 2)), as.character(round(t2_ratio_ifn_il13_adj_L[1], 2)), as.character(round(t2_ratio_il12_il17_adj_L[1], 2)),
                as.character(round(t2_ratio_ifn_il17_adj_L[1], 2)), as.character(round(t2_ratio_il12_il21_adj_L[1], 2)), as.character(round(t2_ratio_ifn_il21_adj_L[1], 2)),
                as.character(round(t2_ratio_pro_il10_adj_L[1], 2)), as.character(round(t2_ratio_th1_il10_adj_L[1], 2)), as.character(round(t2_ratio_th2_il10_adj_L[1], 2)),
                as.character(round(t2_ratio_th17_il10_adj_L[1], 2)), as.character(round(t2_ratio_th1_th2_adj_L[1], 2)), as.character(round(t2_ratio_th1_th17_adj_L[1], 2)))

adjpvals10 <- c("Unadjusted P-value",
                as.character(round(t2_il1_adj_L[6], 2)), as.character(round(t2_il6_adj_L[6], 2)), as.character(round(t2_tnf_adj_L[6], 2)),
                as.character(round(t2_crp_adj_L[6], 2)), as.character(round(t2_il12_adj_L[6], 2)), as.character(round(t2_ifn_adj_L[6], 2)),
                as.character(round(t2_il4_adj_L[6], 2)), as.character(round(t2_il5_adj_L[6], 2)), as.character(round(t2_il13_adj_L[6], 2)),
                as.character(round(t2_il17_adj_L[6], 2)), as.character(round(t2_il21_adj_L[6], 2)), as.character(round(t2_il10_adj_L[6], 2)),
                as.character(round(t2_il2_adj_L[6], 2)), as.character(round(t2_gmc_adj_L[6], 2)), as.character(round(t2_agp_adj_L[6], 2)),
                as.character(round(t2_igf_adj_L[6], 2)), as.character(round(t2_ratio_il1_il10_adj_L[6], 2)),
                as.character(round(t2_ratio_il6_il10_adj_L[6], 2)), as.character(round(t2_ratio_tnf_il10_adj_L[6], 2)), as.character(round(t2_ratio_il12_il10_adj_L[6], 2)),
                as.character(round(t2_ratio_ifn_il10_adj_L[6], 2)), as.character(round(t2_ratio_il4_il10_adj_L[6], 2)), as.character(round(t2_ratio_il5_il10_adj_L[6], 2)),
                as.character(round(t2_ratio_il13_il10_adj_L[6], 2)), as.character(round(t2_ratio_il17_il10_adj_L[6], 2)), as.character(round(t2_ratio_il21_il10_adj_L[6], 2)),
                as.character(round(t2_ratio_il2_il10_adj_L[6], 2)), as.character(round(t2_ratio_gmc_il10_adj_L[6], 2)), as.character(round(t2_ratio_il12_il4_adj_L[6], 2)),
                as.character(round(t2_ratio_ifn_il4_adj_L[6], 2)), as.character(round(t2_ratio_il12_il5_adj_L[6], 2)), as.character(round(t2_ratio_ifn_il5_adj_L[6], 2)),
                as.character(round(t2_ratio_il12_il13_adj_L[6], 2)), as.character(round(t2_ratio_ifn_il13_adj_L[6], 2)), as.character(round(t2_ratio_il12_il17_adj_L[6], 2)),
                as.character(round(t2_ratio_ifn_il17_adj_L[6], 2)), as.character(round(t2_ratio_il12_il21_adj_L[6], 2)), as.character(round(t2_ratio_ifn_il21_adj_L[6], 2)),
                as.character(round(t2_ratio_pro_il10_adj_L[6], 2)), as.character(round(t2_ratio_th1_il10_adj_L[6], 2)), as.character(round(t2_ratio_th2_il10_adj_L[6], 2)),
                as.character(round(t2_ratio_th17_il10_adj_L[6], 2)), as.character(round(t2_ratio_th1_th2_adj_L[6], 2)), as.character(round(t2_ratio_th1_th17_adj_L[6], 2)))

adjbonpvals10 <- c("Bonferroni P-value", 
                   as.character(bonpval(t2_il1_adj_L[6])), as.character(bonpval(t2_il6_adj_L[6])), as.character(bonpval(t2_tnf_adj_L[6])),
                   as.character(bonpval(t2_crp_adj_L[6])), as.character(bonpval(t2_il12_adj_L[6])), as.character(bonpval(t2_ifn_adj_L[6])),
                   as.character(bonpval(t2_il4_adj_L[6])), as.character(bonpval(t2_il5_adj_L[6])), as.character(bonpval(t2_il13_adj_L[6])),
                   as.character(bonpval(t2_il17_adj_L[6])), as.character(bonpval(t2_il21_adj_L[6])), as.character(bonpval(t2_il10_adj_L[6])),
                   as.character(bonpval(t2_il2_adj_L[6])), as.character(bonpval(t2_gmc_adj_L[6])), as.character(bonpval(t2_agp_adj_L[6])),
                   as.character(bonpval(t2_igf_adj_L[6])), as.character(bonpval(t2_ratio_il1_il10_adj_L[6])),
                   as.character(bonpval(t2_ratio_il6_il10_adj_L[6])), as.character(bonpval(t2_ratio_tnf_il10_adj_L[6])), as.character(bonpval(t2_ratio_il12_il10_adj_L[6])),
                   as.character(bonpval(t2_ratio_ifn_il10_adj_L[6])), as.character(bonpval(t2_ratio_il4_il10_adj_L[6])), as.character(bonpval(t2_ratio_il5_il10_adj_L[6])),
                   as.character(bonpval(t2_ratio_il13_il10_adj_L[6])), as.character(bonpval(t2_ratio_il17_il10_adj_L[6])), as.character(bonpval(t2_ratio_il21_il10_adj_L[6])),
                   as.character(bonpval(t2_ratio_il2_il10_adj_L[6])), as.character(bonpval(t2_ratio_gmc_il10_adj_L[6])), as.character(bonpval(t2_ratio_il12_il4_adj_L[6])),
                   as.character(bonpval(t2_ratio_ifn_il4_adj_L[6])), as.character(bonpval(t2_ratio_il12_il5_adj_L[6])), as.character(bonpval(t2_ratio_ifn_il5_adj_L[6])),
                   as.character(bonpval(t2_ratio_il12_il13_adj_L[6])), as.character(bonpval(t2_ratio_ifn_il13_adj_L[6])), as.character(bonpval(t2_ratio_il12_il17_adj_L[6])),
                   as.character(bonpval(t2_ratio_ifn_il17_adj_L[6])), as.character(bonpval(t2_ratio_il12_il21_adj_L[6])), as.character(bonpval(t2_ratio_ifn_il21_adj_L[6])),
                   as.character(bonpval(t2_ratio_pro_il10_adj_L[6])), as.character(bonpval(t2_ratio_th1_il10_adj_L[6])), as.character(bonpval(t2_ratio_th2_il10_adj_L[6])),
                   as.character(bonpval(t2_ratio_th17_il10_adj_L[6])), as.character(bonpval(t2_ratio_th1_th2_adj_L[6])), as.character(bonpval(t2_ratio_th1_th17_adj_L[6])))

# Table S8: P-values of treatment estimates, unadjusted, and adjusted for multiple testing 
# (by controlling family-wise error rate using the Bonferroni correction) at age 14 months. 
tbls10 <- data.table(
  " " = outcomes10,
  "Unadjusted Analysis" = unadjdiffs10, 
  " " = unadjpvals10,
  " " = unadjbonpvals10, 
  "Adjusted Analysis" = adjdiffs10,
  " " = adjpvals10,
  " " = adjbonpvals10
  )

write.csv(tbls10, file=here('tables/supplementary/immune_supptable8.csv'))
print(xtable(tbls10), type="html", file=here("tables/supplementary/immune_supptable8.html"))



### TABLE S9 ####
outcomes11<-c("Outcome", "Ln IL-1\u03b2 (pg/ml)", "Ln Il-6 (pg/ml)", "Ln TNF-\u03b1 (pg/ml)", "Ln IL-12 (pg/ml)",
              "Ln IFN-\u03b3 (pg/ml)", "Ln IL-4 (pg/ml)", "Ln IL-5 (pg/ml)", "Ln IL-13 (pg/ml)", "Ln IL-17A (pg/ml)", 
              "Ln IL-21 (pg/ml)", "Ln IL-10 (pg/ml)", "Ln IL-2 (pg/ml)", "Ln GM-CSF (pg/ml)", "Ln IGF-1 (\u03bcg/L)",
              "Ln IL-1\u03b2/IL-10", "Ln IL-6/IL-10", "Ln TNF-\u03b1/IL-10", "Ln IL-12/IL-10", "Ln IFN-\u03b3/IL-10",
              "Ln IL-4/IL-10", "Ln IL-5/IL-10", "Ln IL-13/IL-10", "Ln IL-17A/IL-10", "Ln IL-21/IL-10",
              "Ln IL-2/IL-10", "Ln GM-CSF/IL-10", "Ln IL-12/IL-4", "Ln IFN-\u03b3/IL-4", "Ln IL-12/IL-5", "Ln IFN-\u03b3/IL-5",
              "Ln IL-12/IL-13", "Ln IFN-\u03b3/IL-13", "Ln IL-12/IL-17A", "Ln IFN-\u03b3/IL-17A", "Ln IL-12/IL-21", "Ln IFN-\u03b3/IL-21",
              "Ln Pro-inflammatory cytokines/IL-10", "Ln Th1/IL-10", "Ln Th2/IL-10", "Ln Th17/IL-10",
              "Ln Th1/Th2", "Ln Th1/Th17")

unadjdiffs11 <- c("Unadjusted Difference", 
                  as.character(round(t3_il1_unadj_L[1], 2)), as.character(round(t3_il6_unadj_L[1], 2)), as.character(round(t3_tnf_unadj_L[1], 2)),
                  as.character(round(t3_il12_unadj_L[1], 2)), as.character(round(t3_ifn_unadj_L[1], 2)),
                  as.character(round(t3_il4_unadj_L[1], 2)), as.character(round(t3_il5_unadj_L[1], 2)), as.character(round(t3_il13_unadj_L[1], 2)),
                  as.character(round(t3_il17_unadj_L[1], 2)), as.character(round(t3_il21_unadj_L[1], 2)), as.character(round(t3_il10_unadj_L[1], 2)),
                  as.character(round(t3_il2_unadj_L[1], 2)), as.character(round(t3_gmc_unadj_L[1], 2)),
                  as.character(round(t3_igf_unadj_L[1], 2)), as.character(round(t3_ratio_il1_il10_unadj_L[1], 2)),
                  as.character(round(t3_ratio_il6_il10_unadj_L[1], 2)), as.character(round(t3_ratio_tnf_il10_unadj_L[1], 2)), as.character(round(t3_ratio_il12_il10_unadj_L[1], 2)),
                  as.character(round(t3_ratio_ifn_il10_unadj_L[1], 2)), as.character(round(t3_ratio_il4_il10_unadj_L[1], 2)), as.character(round(t3_ratio_il5_il10_unadj_L[1], 2)),
                  as.character(round(t3_ratio_il13_il10_unadj_L[1], 2)), as.character(round(t3_ratio_il17_il10_unadj_L[1], 2)), as.character(round(t3_ratio_il21_il10_unadj_L[1], 2)),
                  as.character(round(t3_ratio_il2_il10_unadj_L[1], 2)), as.character(round(t3_ratio_gmc_il10_unadj_L[1], 2)), as.character(round(t3_ratio_il12_il4_unadj_L[1], 2)),
                  as.character(round(t3_ratio_ifn_il4_unadj_L[1], 2)), as.character(round(t3_ratio_il12_il5_unadj_L[1], 2)), as.character(round(t3_ratio_ifn_il5_unadj_L[1], 2)),
                  as.character(round(t3_ratio_il12_il13_unadj_L[1], 2)), as.character(round(t3_ratio_ifn_il13_unadj_L[1], 2)), as.character(round(t3_ratio_il12_il17_unadj_L[1], 2)),
                  as.character(round(t3_ratio_ifn_il17_unadj_L[1], 2)), as.character(round(t3_ratio_il12_il21_unadj_L[1], 2)), as.character(round(t3_ratio_ifn_il21_unadj_L[1], 2)),
                  as.character(round(t3_ratio_pro_il10_unadj_L[1], 2)), as.character(round(t3_ratio_th1_il10_unadj_L[1], 2)), as.character(round(t3_ratio_th2_il10_unadj_L[1], 2)),
                  as.character(round(t3_ratio_th17_il10_unadj_L[1], 2)), as.character(round(t3_ratio_th1_th2_unadj_L[1], 2)), as.character(round(t3_ratio_th1_th17_unadj_L[1], 2))
)

unadjpvals11 <- c("Unadjusted P-value", 
                  as.character(round(t3_il1_unadj_L[6], 2)), as.character(round(t3_il6_unadj_L[6], 2)), as.character(round(t3_tnf_unadj_L[6], 2)),
                  as.character(round(t3_il12_unadj_L[6], 2)), as.character(round(t3_ifn_unadj_L[6], 2)),
                  as.character(round(t3_il4_unadj_L[6], 2)), as.character(round(t3_il5_unadj_L[6], 2)), as.character(round(t3_il13_unadj_L[6], 2)),
                  as.character(round(t3_il17_unadj_L[6], 2)), as.character(round(t3_il21_unadj_L[6], 2)), as.character(round(t3_il10_unadj_L[6], 2)),
                  as.character(round(t3_il2_unadj_L[6], 2)), as.character(round(t3_gmc_unadj_L[6], 2)),
                  as.character(round(t3_igf_unadj_L[6], 2)), as.character(round(t3_ratio_il1_il10_unadj_L[6], 2)),
                  as.character(round(t3_ratio_il6_il10_unadj_L[6], 2)), as.character(round(t3_ratio_tnf_il10_unadj_L[6], 2)), as.character(round(t3_ratio_il12_il10_unadj_L[6], 2)),
                  as.character(round(t3_ratio_ifn_il10_unadj_L[6], 2)), as.character(round(t3_ratio_il4_il10_unadj_L[6], 2)), as.character(round(t3_ratio_il5_il10_unadj_L[6], 2)),
                  as.character(round(t3_ratio_il13_il10_unadj_L[6], 2)), as.character(round(t3_ratio_il17_il10_unadj_L[6], 2)), as.character(round(t3_ratio_il21_il10_unadj_L[6], 2)),
                  as.character(round(t3_ratio_il2_il10_unadj_L[6], 2)), as.character(round(t3_ratio_gmc_il10_unadj_L[6], 2)), as.character(round(t3_ratio_il12_il4_unadj_L[6], 2)),
                  as.character(round(t3_ratio_ifn_il4_unadj_L[6], 2)), as.character(round(t3_ratio_il12_il5_unadj_L[6], 2)), as.character(round(t3_ratio_ifn_il5_unadj_L[6], 2)),
                  as.character(round(t3_ratio_il12_il13_unadj_L[6], 2)), as.character(round(t3_ratio_ifn_il13_unadj_L[6], 2)), as.character(round(t3_ratio_il12_il17_unadj_L[6], 2)),
                  as.character(round(t3_ratio_ifn_il17_unadj_L[6], 2)), as.character(round(t3_ratio_il12_il21_unadj_L[6], 2)), as.character(round(t3_ratio_ifn_il21_unadj_L[6], 2)),
                  as.character(round(t3_ratio_pro_il10_unadj_L[6], 2)), as.character(round(t3_ratio_th1_il10_unadj_L[6], 2)), as.character(round(t3_ratio_th2_il10_unadj_L[6], 2)),
                  as.character(round(t3_ratio_th17_il10_unadj_L[6], 2)), as.character(round(t3_ratio_th1_th2_unadj_L[6], 2)), as.character(round(t3_ratio_th1_th17_unadj_L[6], 2)))

unadjbonpvals11 <- c("Bonferroni P-value",
                     as.character(bonpval(t3_il1_unadj_L[6])), as.character(bonpval(t3_il6_unadj_L[6])), as.character(bonpval(t3_tnf_unadj_L[6])),
                     as.character(bonpval(t3_il12_unadj_L[6])), as.character(bonpval(t3_ifn_unadj_L[6])),
                     as.character(bonpval(t3_il4_unadj_L[6])), as.character(bonpval(t3_il5_unadj_L[6])), as.character(bonpval(t3_il13_unadj_L[6])),
                     as.character(bonpval(t3_il17_unadj_L[6])), as.character(bonpval(t3_il21_unadj_L[6])), as.character(bonpval(t3_il10_unadj_L[6])),
                     as.character(bonpval(t3_il2_unadj_L[6])), as.character(bonpval(t3_gmc_unadj_L[6])),
                     as.character(bonpval(t3_igf_unadj_L[6])), as.character(bonpval(t3_ratio_il1_il10_unadj_L[6])),
                     as.character(bonpval(t3_ratio_il6_il10_unadj_L[6])), as.character(bonpval(t3_ratio_tnf_il10_unadj_L[6])), as.character(bonpval(t3_ratio_il12_il10_unadj_L[6])),
                     as.character(bonpval(t3_ratio_ifn_il10_unadj_L[6])), as.character(bonpval(t3_ratio_il4_il10_unadj_L[6])), as.character(bonpval(t3_ratio_il5_il10_unadj_L[6])),
                     as.character(bonpval(t3_ratio_il13_il10_unadj_L[6])), as.character(bonpval(t3_ratio_il17_il10_unadj_L[6])), as.character(bonpval(t3_ratio_il21_il10_unadj_L[6])),
                     as.character(bonpval(t3_ratio_il2_il10_unadj_L[6])), as.character(bonpval(t3_ratio_gmc_il10_unadj_L[6])), as.character(bonpval(t3_ratio_il12_il4_unadj_L[6])),
                     as.character(bonpval(t3_ratio_ifn_il4_unadj_L[6])), as.character(bonpval(t3_ratio_il12_il5_unadj_L[6])), as.character(bonpval(t3_ratio_ifn_il5_unadj_L[6])),
                     as.character(bonpval(t3_ratio_il12_il13_unadj_L[6])), as.character(bonpval(t3_ratio_ifn_il13_unadj_L[6])), as.character(bonpval(t3_ratio_il12_il17_unadj_L[6])),
                     as.character(bonpval(t3_ratio_ifn_il17_unadj_L[6])), as.character(bonpval(t3_ratio_il12_il21_unadj_L[6])), as.character(bonpval(t3_ratio_ifn_il21_unadj_L[6])),
                     as.character(bonpval(t3_ratio_pro_il10_unadj_L[6])), as.character(bonpval(t3_ratio_th1_il10_unadj_L[6])), as.character(bonpval(t3_ratio_th2_il10_unadj_L[6])),
                     as.character(bonpval(t3_ratio_th17_il10_unadj_L[6])), as.character(bonpval(t3_ratio_th1_th2_unadj_L[6])), as.character(bonpval(t3_ratio_th1_th17_unadj_L[6])))

adjdiffs11 <- c("Adjusted Difference",  
                as.character(round(t3_il1_adj_L[1], 2)), as.character(round(t3_il6_adj_L[1], 2)), as.character(round(t3_tnf_adj_L[1], 2)),
                as.character(round(t3_il12_adj_L[1], 2)), as.character(round(t3_ifn_adj_L[1], 2)),
                as.character(round(t3_il4_adj_L[1], 2)), as.character(round(t3_il5_adj_L[1], 2)), as.character(round(t3_il13_adj_L[1], 2)),
                as.character(round(t3_il17_adj_L[1], 2)), as.character(round(t3_il21_adj_L[1], 2)), as.character(round(t3_il10_adj_L[1], 2)),
                as.character(round(t3_il2_adj_L[1], 2)), as.character(round(t3_gmc_adj_L[1], 2)), 
                as.character(round(t3_igf_adj_L[1], 2)), as.character(round(t3_ratio_il1_il10_adj_L[1], 2)),
                as.character(round(t3_ratio_il6_il10_adj_L[1], 2)), as.character(round(t3_ratio_tnf_il10_adj_L[1], 2)), as.character(round(t3_ratio_il12_il10_adj_L[1], 2)),
                as.character(round(t3_ratio_ifn_il10_adj_L[1], 2)), as.character(round(t3_ratio_il4_il10_adj_L[1], 2)), as.character(round(t3_ratio_il5_il10_adj_L[1], 2)),
                as.character(round(t3_ratio_il13_il10_adj_L[1], 2)), as.character(round(t3_ratio_il17_il10_adj_L[1], 2)), as.character(round(t3_ratio_il21_il10_adj_L[1], 2)),
                as.character(round(t3_ratio_il2_il10_adj_L[1], 2)), as.character(round(t3_ratio_gmc_il10_adj_L[1], 2)), as.character(round(t3_ratio_il12_il4_adj_L[1], 2)),
                as.character(round(t3_ratio_ifn_il4_adj_L[1], 2)), as.character(round(t3_ratio_il12_il5_adj_L[1], 2)), as.character(round(t3_ratio_ifn_il5_adj_L[1], 2)),
                as.character(round(t3_ratio_il12_il13_adj_L[1], 2)), as.character(round(t3_ratio_ifn_il13_adj_L[1], 2)), as.character(round(t3_ratio_il12_il17_adj_L[1], 2)),
                as.character(round(t3_ratio_ifn_il17_adj_L[1], 2)), as.character(round(t3_ratio_il12_il21_adj_L[1], 2)), as.character(round(t3_ratio_ifn_il21_adj_L[1], 2)),
                as.character(round(t3_ratio_pro_il10_adj_L[1], 2)), as.character(round(t3_ratio_th1_il10_adj_L[1], 2)), as.character(round(t3_ratio_th2_il10_adj_L[1], 2)),
                as.character(round(t3_ratio_th17_il10_adj_L[1], 2)), as.character(round(t3_ratio_th1_th2_adj_L[1], 2)), as.character(round(t3_ratio_th1_th17_adj_L[1], 2)))

adjpvals11 <- c("Unadjusted P-value",
                as.character(round(t3_il1_adj_L[6], 2)), as.character(round(t3_il6_adj_L[6], 2)), as.character(round(t3_tnf_adj_L[6], 2)),
                as.character(round(t3_il12_adj_L[6], 2)), as.character(round(t3_ifn_adj_L[6], 2)),
                as.character(round(t3_il4_adj_L[6], 2)), as.character(round(t3_il5_adj_L[6], 2)), as.character(round(t3_il13_adj_L[6], 2)),
                as.character(round(t3_il17_adj_L[6], 2)), as.character(round(t3_il21_adj_L[6], 2)), as.character(round(t3_il10_adj_L[6], 2)),
                as.character(round(t3_il2_adj_L[6], 2)), as.character(round(t3_gmc_adj_L[6], 2)),
                as.character(round(t3_igf_adj_L[6], 2)), as.character(round(t3_ratio_il1_il10_adj_L[6], 2)),
                as.character(round(t3_ratio_il6_il10_adj_L[6], 2)), as.character(round(t3_ratio_tnf_il10_adj_L[6], 2)), as.character(round(t3_ratio_il12_il10_adj_L[6], 2)),
                as.character(round(t3_ratio_ifn_il10_adj_L[6], 2)), as.character(round(t3_ratio_il4_il10_adj_L[6], 2)), as.character(round(t3_ratio_il5_il10_adj_L[6], 2)),
                as.character(round(t3_ratio_il13_il10_adj_L[6], 2)), as.character(round(t3_ratio_il17_il10_adj_L[6], 2)), as.character(round(t3_ratio_il21_il10_adj_L[6], 2)),
                as.character(round(t3_ratio_il2_il10_adj_L[6], 2)), as.character(round(t3_ratio_gmc_il10_adj_L[6], 2)), as.character(round(t3_ratio_il12_il4_adj_L[6], 2)),
                as.character(round(t3_ratio_ifn_il4_adj_L[6], 2)), as.character(round(t3_ratio_il12_il5_adj_L[6], 2)), as.character(round(t3_ratio_ifn_il5_adj_L[6], 2)),
                as.character(round(t3_ratio_il12_il13_adj_L[6], 2)), as.character(round(t3_ratio_ifn_il13_adj_L[6], 2)), as.character(round(t3_ratio_il12_il17_adj_L[6], 2)),
                as.character(round(t3_ratio_ifn_il17_adj_L[6], 2)), as.character(round(t3_ratio_il12_il21_adj_L[6], 2)), as.character(round(t3_ratio_ifn_il21_adj_L[6], 2)),
                as.character(round(t3_ratio_pro_il10_adj_L[6], 2)), as.character(round(t3_ratio_th1_il10_adj_L[6], 2)), as.character(round(t3_ratio_th2_il10_adj_L[6], 2)),
                as.character(round(t3_ratio_th17_il10_adj_L[6], 2)), as.character(round(t3_ratio_th1_th2_adj_L[6], 2)), as.character(round(t3_ratio_th1_th17_adj_L[6], 2)))

adjbonpvals11 <- c("Bonferroni P-value",
                   as.character(bonpval(t3_il1_adj_L[6])), as.character(bonpval(t3_il6_adj_L[6])), as.character(bonpval(t3_tnf_adj_L[6])),
                   as.character(bonpval(t3_il12_adj_L[6])), as.character(bonpval(t3_ifn_adj_L[6])),
                   as.character(bonpval(t3_il4_adj_L[6])), as.character(bonpval(t3_il5_adj_L[6])), as.character(bonpval(t3_il13_adj_L[6])),
                   as.character(bonpval(t3_il17_adj_L[6])), as.character(bonpval(t3_il21_adj_L[6])), as.character(bonpval(t3_il10_adj_L[6])),
                   as.character(bonpval(t3_il2_adj_L[6])), as.character(bonpval(t3_gmc_adj_L[6])),
                   as.character(bonpval(t3_igf_adj_L[6])), as.character(bonpval(t3_ratio_il1_il10_adj_L[6])),
                   as.character(bonpval(t3_ratio_il6_il10_adj_L[6])), as.character(bonpval(t3_ratio_tnf_il10_adj_L[6])), as.character(bonpval(t3_ratio_il12_il10_adj_L[6])),
                   as.character(bonpval(t3_ratio_ifn_il10_adj_L[6])), as.character(bonpval(t3_ratio_il4_il10_adj_L[6])), as.character(bonpval(t3_ratio_il5_il10_adj_L[6])),
                   as.character(bonpval(t3_ratio_il13_il10_adj_L[6])), as.character(bonpval(t3_ratio_il17_il10_adj_L[6])), as.character(bonpval(t3_ratio_il21_il10_adj_L[6])),
                   as.character(bonpval(t3_ratio_il2_il10_adj_L[6])), as.character(bonpval(t3_ratio_gmc_il10_adj_L[6])), as.character(bonpval(t3_ratio_il12_il4_adj_L[6])),
                   as.character(bonpval(t3_ratio_ifn_il4_adj_L[6])), as.character(bonpval(t3_ratio_il12_il5_adj_L[6])), as.character(bonpval(t3_ratio_ifn_il5_adj_L[6])),
                   as.character(bonpval(t3_ratio_il12_il13_adj_L[6])), as.character(bonpval(t3_ratio_ifn_il13_adj_L[6])), as.character(bonpval(t3_ratio_il12_il17_adj_L[6])),
                   as.character(bonpval(t3_ratio_ifn_il17_adj_L[6])), as.character(bonpval(t3_ratio_il12_il21_adj_L[6])), as.character(bonpval(t3_ratio_ifn_il21_adj_L[6])),
                   as.character(bonpval(t3_ratio_pro_il10_adj_L[6])), as.character(bonpval(t3_ratio_th1_il10_adj_L[6])), as.character(bonpval(t3_ratio_th2_il10_adj_L[6])),
                   as.character(bonpval(t3_ratio_th17_il10_adj_L[6])), as.character(bonpval(t3_ratio_th1_th2_adj_L[6])), as.character(bonpval(t3_ratio_th1_th17_adj_L[6])))

# Table S9: P-values of treatment estimates, unadjusted, and adjusted for multiple testing 
# (by controlling family-wise error rate using the Bonferroni correction) at age 18 months. 
tbls11 <- data.table(
  " " = outcomes11,
  "Unadjusted Analysis" = unadjdiffs11, 
  " " = unadjpvals11,
  " " = unadjbonpvals11, 
  "Adjusted Analysis" = adjdiffs11,
  " " = adjpvals11,
  " " = adjbonpvals11
)

write.csv(tbls11, file=here('tables/supplementary/immune_supptable9.csv'))
print(xtable(tbls11), type="html", file=here("tables/supplementary/immune_supptable9.html"))


#### TABLE 10 ####
outcomes12<-c("Outcome", "Ln \u0394IL-1\u03b2 (pg/ml)", "Ln \u0394Il-6 (pg/ml)", "Ln \u0394TNF-\u03b1 (pg/ml)", "Ln \u0394IL-12 (pg/ml)",
              "Ln \u0394IFN-\u03b3 (pg/ml)", "Ln \u0394IL-4 (pg/ml)", "Ln \u0394IL-5 (pg/ml)", "Ln \u0394IL-13 (pg/ml)", "Ln \u0394IL-17A (pg/ml)", 
              "Ln \u0394IL-21 (pg/ml)", "Ln \u0394IL-10 (pg/ml)", "Ln \u0394IL-2 (pg/ml)", "Ln \u0394GM-CSF (pg/ml)", "Ln \u0394IGF-1 (\u03bcg/L)",
              "Ln \u0394IL-1\u03b2/IL-10", "Ln \u0394IL-6/IL-10", "Ln \u0394TNF-\u03b1/IL-10", "Ln \u0394IL-12/IL-10", "Ln \u0394IFN-\u03b3/IL-10",
              "Ln \u0394IL-4/IL-10", "Ln \u0394IL-5/IL-10", "Ln \u0394IL-13/IL-10", "Ln \u0394IL-17A/IL-10", "Ln \u0394IL-21/IL-10",
              "Ln \u0394IL-2/IL-10", "Ln \u0394GM-CSF/IL-10", "Ln \u0394IL-12/IL-4", "Ln \u0394IFN-\u03b3/IL-4", "Ln \u0394IL-12/IL-5", "Ln \u0394IFN-\u03b3/IL-5",
              "Ln \u0394IL-12/IL-13", "Ln \u0394IFN-\u03b3/IL-13", "Ln \u0394IL-12/IL-17A", "Ln \u0394IFN-\u03b3/IL-17A", "Ln \u0394IL-12/IL-21", "Ln \u0394IFN-\u03b3/IL-21",
              "Ln \u0394Pro-inflammatory cytokines/IL-10", "Ln \u0394Th1/IL-10", "Ln \u0394Th2/IL-10", "Ln \u0394Th17/IL-10",
              "Ln \u0394Th1/Th2", "Ln \u0394Th1/Th17")

unadjdiffs12 <- c("Unadjusted Difference", 
                  as.character(round(d23_ln_il1_unadj_L[1], 2)), as.character(round(d23_ln_il6_unadj_L[1], 2)), as.character(round(d23_ln_tnf_unadj_L[1], 2)),
                  as.character(round(d23_ln_il12_unadj_L[1], 2)), as.character(round(d23_ln_ifn_unadj_L[1], 2)),
                  as.character(round(d23_ln_il4_unadj_L[1], 2)), as.character(round(d23_ln_il5_unadj_L[1], 2)), as.character(round(d23_ln_il13_unadj_L[1], 2)),
                  as.character(round(d23_ln_il17_unadj_L[1], 2)), as.character(round(d23_ln_il21_unadj_L[1], 2)), as.character(round(d23_ln_il10_unadj_L[1], 2)),
                  as.character(round(d23_ln_il2_unadj_L[1], 2)), as.character(round(d23_ln_gmc_unadj_L[1], 2)),
                  as.character(round(d23_ln_igf_unadj_L[1], 2)), as.character(round(d23_ratio_il1_il10_unadj_L[1], 2)),
                  as.character(round(d23_ratio_il6_il10_unadj_L[1], 2)), as.character(round(d23_ratio_tnf_il10_unadj_L[1], 2)), as.character(round(d23_ratio_il12_il10_unadj_L[1], 2)),
                  as.character(round(d23_ratio_ifn_il10_unadj_L[1], 2)), as.character(round(d23_ratio_il4_il10_unadj_L[1], 2)), as.character(round(d23_ratio_il5_il10_unadj_L[1], 2)),
                  as.character(round(d23_ratio_il13_il10_unadj_L[1], 2)), as.character(round(d23_ratio_il17_il10_unadj_L[1], 2)), as.character(round(d23_ratio_il21_il10_unadj_L[1], 2)),
                  as.character(round(d23_ratio_il2_il10_unadj_L[1], 2)), as.character(round(d23_ratio_gmc_il10_unadj_L[1], 2)), as.character(round(d23_ratio_il12_il4_unadj_L[1], 2)),
                  as.character(round(d23_ratio_ifn_il4_unadj_L[1], 2)), as.character(round(d23_ratio_il12_il5_unadj_L[1], 2)), as.character(round(d23_ratio_ifn_il5_unadj_L[1], 2)),
                  as.character(round(d23_ratio_il12_il13_unadj_L[1], 2)), as.character(round(d23_ratio_ifn_il13_unadj_L[1], 2)), as.character(round(d23_ratio_il12_il17_unadj_L[1], 2)),
                  as.character(round(d23_ratio_ifn_il17_unadj_L[1], 2)), as.character(round(d23_ratio_il12_il21_unadj_L[1], 2)), as.character(round(d23_ratio_ifn_il21_unadj_L[1], 2)),
                  as.character(round(d23_ratio_pro_il10_unadj_L[1], 2)), as.character(round(d23_ratio_th1_il10_unadj_L[1], 2)), as.character(round(d23_ratio_th2_il10_unadj_L[1], 2)),
                  as.character(round(d23_ratio_th17_il10_unadj_L[1], 2)), as.character(round(d23_ratio_th1_th2_unadj_L[1], 2)), as.character(round(d23_ratio_th1_th17_unadj_L[1], 2))
)

unadjpvals12 <- c("Unadjusted P-value", 
                  as.character(round(d23_ln_il1_unadj_L[6], 2)), as.character(round(d23_ln_il6_unadj_L[6], 2)), as.character(round(d23_ln_tnf_unadj_L[6], 2)),
                  as.character(round(d23_ln_il12_unadj_L[6], 2)), as.character(round(d23_ln_ifn_unadj_L[6], 2)),
                  as.character(round(d23_ln_il4_unadj_L[6], 2)), as.character(round(d23_ln_il5_unadj_L[6], 2)), as.character(round(d23_ln_il13_unadj_L[6], 2)),
                  as.character(round(d23_ln_il17_unadj_L[6], 2)), as.character(round(d23_ln_il21_unadj_L[6], 2)), as.character(round(d23_ln_il10_unadj_L[6], 2)),
                  as.character(round(d23_ln_il2_unadj_L[6], 2)), as.character(round(d23_ln_gmc_unadj_L[6], 2)),
                  as.character(round(d23_ln_igf_unadj_L[6], 2)), as.character(round(d23_ratio_il1_il10_unadj_L[6], 2)),
                  as.character(round(d23_ratio_il6_il10_unadj_L[6], 2)), as.character(round(d23_ratio_tnf_il10_unadj_L[6], 2)), as.character(round(d23_ratio_il12_il10_unadj_L[6], 2)),
                  as.character(round(d23_ratio_ifn_il10_unadj_L[6], 2)), as.character(round(d23_ratio_il4_il10_unadj_L[6], 2)), as.character(round(d23_ratio_il5_il10_unadj_L[6], 2)),
                  as.character(round(d23_ratio_il13_il10_unadj_L[6], 2)), as.character(round(d23_ratio_il17_il10_unadj_L[6], 2)), as.character(round(d23_ratio_il21_il10_unadj_L[6], 2)),
                  as.character(round(d23_ratio_il2_il10_unadj_L[6], 2)), as.character(round(d23_ratio_gmc_il10_unadj_L[6], 2)), as.character(round(d23_ratio_il12_il4_unadj_L[6], 2)),
                  as.character(round(d23_ratio_ifn_il4_unadj_L[6], 2)), as.character(round(d23_ratio_il12_il5_unadj_L[6], 2)), as.character(round(d23_ratio_ifn_il5_unadj_L[6], 2)),
                  as.character(round(d23_ratio_il12_il13_unadj_L[6], 2)), as.character(round(d23_ratio_ifn_il13_unadj_L[6], 2)), as.character(round(d23_ratio_il12_il17_unadj_L[6], 2)),
                  as.character(round(d23_ratio_ifn_il17_unadj_L[6], 2)), as.character(round(d23_ratio_il12_il21_unadj_L[6], 2)), as.character(round(d23_ratio_ifn_il21_unadj_L[6], 2)),
                  as.character(round(d23_ratio_pro_il10_unadj_L[6], 2)), as.character(round(d23_ratio_th1_il10_unadj_L[6], 2)), as.character(round(d23_ratio_th2_il10_unadj_L[6], 2)),
                  as.character(round(d23_ratio_th17_il10_unadj_L[6], 2)), as.character(round(d23_ratio_th1_th2_unadj_L[6], 2)), as.character(round(d23_ratio_th1_th17_unadj_L[6], 2)))

unadjbonpvals12 <- c("Bonferroni P-value",
                     as.character(bonpval(d23_ln_il1_unadj_L[6])), as.character(bonpval(d23_ln_il6_unadj_L[6])), as.character(bonpval(d23_ln_tnf_unadj_L[6])),
                     as.character(bonpval(d23_ln_il12_unadj_L[6])), as.character(bonpval(d23_ln_ifn_unadj_L[6])),
                     as.character(bonpval(d23_ln_il4_unadj_L[6])), as.character(bonpval(d23_ln_il5_unadj_L[6])), as.character(bonpval(d23_ln_il13_unadj_L[6])),
                     as.character(bonpval(d23_ln_il17_unadj_L[6])), as.character(bonpval(d23_ln_il21_unadj_L[6])), as.character(bonpval(d23_ln_il10_unadj_L[6])),
                     as.character(bonpval(d23_ln_il2_unadj_L[6])), as.character(bonpval(d23_ln_gmc_unadj_L[6])),
                     as.character(bonpval(d23_ln_igf_unadj_L[6])), as.character(bonpval(d23_ratio_il1_il10_unadj_L[6])),
                     as.character(bonpval(d23_ratio_il6_il10_unadj_L[6])), as.character(bonpval(d23_ratio_tnf_il10_unadj_L[6])), as.character(bonpval(d23_ratio_il12_il10_unadj_L[6])),
                     as.character(bonpval(d23_ratio_ifn_il10_unadj_L[6])), as.character(bonpval(d23_ratio_il4_il10_unadj_L[6])), as.character(bonpval(d23_ratio_il5_il10_unadj_L[6])),
                     as.character(bonpval(d23_ratio_il13_il10_unadj_L[6])), as.character(bonpval(d23_ratio_il17_il10_unadj_L[6])), as.character(bonpval(d23_ratio_il21_il10_unadj_L[6])),
                     as.character(bonpval(d23_ratio_il2_il10_unadj_L[6])), as.character(bonpval(d23_ratio_gmc_il10_unadj_L[6])), as.character(bonpval(d23_ratio_il12_il4_unadj_L[6])),
                     as.character(bonpval(d23_ratio_ifn_il4_unadj_L[6])), as.character(bonpval(d23_ratio_il12_il5_unadj_L[6])), as.character(bonpval(d23_ratio_ifn_il5_unadj_L[6])),
                     as.character(bonpval(d23_ratio_il12_il13_unadj_L[6])), as.character(bonpval(d23_ratio_ifn_il13_unadj_L[6])), as.character(bonpval(d23_ratio_il12_il17_unadj_L[6])),
                     as.character(bonpval(d23_ratio_ifn_il17_unadj_L[6])), as.character(bonpval(d23_ratio_il12_il21_unadj_L[6])), as.character(bonpval(d23_ratio_ifn_il21_unadj_L[6])),
                     as.character(bonpval(d23_ratio_pro_il10_unadj_L[6])), as.character(bonpval(d23_ratio_th1_il10_unadj_L[6])), as.character(bonpval(d23_ratio_th2_il10_unadj_L[6])),
                     as.character(bonpval(d23_ratio_th17_il10_unadj_L[6])), as.character(bonpval(d23_ratio_th1_th2_unadj_L[6])), as.character(bonpval(d23_ratio_th1_th17_unadj_L[6])))

adjdiffs12 <- c("Adjusted Difference",  
                as.character(round(d23_ln_il1_adj_L[1], 2)), as.character(round(d23_ln_il6_adj_L[1], 2)), as.character(round(d23_ln_tnf_adj_L[1], 2)),
                as.character(round(d23_ln_il12_adj_L[1], 2)), as.character(round(d23_ln_ifn_adj_L[1], 2)),
                as.character(round(d23_ln_il4_adj_L[1], 2)), as.character(round(d23_ln_il5_adj_L[1], 2)), as.character(round(d23_ln_il13_adj_L[1], 2)),
                as.character(round(d23_ln_il17_adj_L[1], 2)), as.character(round(d23_ln_il21_adj_L[1], 2)), as.character(round(d23_ln_il10_adj_L[1], 2)),
                as.character(round(d23_ln_il2_adj_L[1], 2)), as.character(round(d23_ln_gmc_adj_L[1], 2)),
                as.character(round(d23_ln_igf_adj_L[1], 2)), as.character(round(d23_ratio_il1_il10_adj_L[1], 2)),
                as.character(round(d23_ratio_il6_il10_adj_L[1], 2)), as.character(round(d23_ratio_tnf_il10_adj_L[1], 2)), as.character(round(d23_ratio_il12_il10_adj_L[1], 2)),
                as.character(round(d23_ratio_ifn_il10_adj_L[1], 2)), as.character(round(d23_ratio_il4_il10_adj_L[1], 2)), as.character(round(d23_ratio_il5_il10_adj_L[1], 2)),
                as.character(round(d23_ratio_il13_il10_adj_L[1], 2)), as.character(round(d23_ratio_il17_il10_adj_L[1], 2)), as.character(round(d23_ratio_il21_il10_adj_L[1], 2)),
                as.character(round(d23_ratio_il2_il10_adj_L[1], 2)), as.character(round(d23_ratio_gmc_il10_adj_L[1], 2)), as.character(round(d23_ratio_il12_il4_adj_L[1], 2)),
                as.character(round(d23_ratio_ifn_il4_adj_L[1], 2)), as.character(round(d23_ratio_il12_il5_adj_L[1], 2)), as.character(round(d23_ratio_ifn_il5_adj_L[1], 2)),
                as.character(round(d23_ratio_il12_il13_adj_L[1], 2)), as.character(round(d23_ratio_ifn_il13_adj_L[1], 2)), as.character(round(d23_ratio_il12_il17_adj_L[1], 2)),
                as.character(round(d23_ratio_ifn_il17_adj_L[1], 2)), as.character(round(d23_ratio_il12_il21_adj_L[1], 2)), as.character(round(d23_ratio_ifn_il21_adj_L[1], 2)),
                as.character(round(d23_ratio_pro_il10_adj_L[1], 2)), as.character(round(d23_ratio_th1_il10_adj_L[1], 2)), as.character(round(d23_ratio_th2_il10_adj_L[1], 2)),
                as.character(round(d23_ratio_th17_il10_adj_L[1], 2)), as.character(round(d23_ratio_th1_th2_adj_L[1], 2)), as.character(round(d23_ratio_th1_th17_adj_L[1], 2)))

adjpvals12 <- c("Unadjusted P-value",
                as.character(round(d23_ln_il1_adj_L[6], 2)), as.character(round(d23_ln_il6_adj_L[6], 2)), as.character(round(d23_ln_tnf_adj_L[6], 2)),
                as.character(round(d23_ln_il12_adj_L[6], 2)), as.character(round(d23_ln_ifn_adj_L[6], 2)),
                as.character(round(d23_ln_il4_adj_L[6], 2)), as.character(round(d23_ln_il5_adj_L[6], 2)), as.character(round(d23_ln_il13_adj_L[6], 2)),
                as.character(round(d23_ln_il17_adj_L[6], 2)), as.character(round(d23_ln_il21_adj_L[6], 2)), as.character(round(d23_ln_il10_adj_L[6], 2)),
                as.character(round(d23_ln_il2_adj_L[6], 2)), as.character(round(d23_ln_gmc_adj_L[6], 2)),
                as.character(round(d23_ln_igf_adj_L[6], 2)), as.character(round(d23_ratio_il1_il10_adj_L[6], 2)),
                as.character(round(d23_ratio_il6_il10_adj_L[6], 2)), as.character(round(d23_ratio_tnf_il10_adj_L[6], 2)), as.character(round(d23_ratio_il12_il10_adj_L[6], 2)),
                as.character(round(d23_ratio_ifn_il10_adj_L[6], 2)), as.character(round(d23_ratio_il4_il10_adj_L[6], 2)), as.character(round(d23_ratio_il5_il10_adj_L[6], 2)),
                as.character(round(d23_ratio_il13_il10_adj_L[6], 2)), as.character(round(d23_ratio_il17_il10_adj_L[6], 2)), as.character(round(d23_ratio_il21_il10_adj_L[6], 2)),
                as.character(round(d23_ratio_il2_il10_adj_L[6], 2)), as.character(round(d23_ratio_gmc_il10_adj_L[6], 2)), as.character(round(d23_ratio_il12_il4_adj_L[6], 2)),
                as.character(round(d23_ratio_ifn_il4_adj_L[6], 2)), as.character(round(d23_ratio_il12_il5_adj_L[6], 2)), as.character(round(d23_ratio_ifn_il5_adj_L[6], 2)),
                as.character(round(d23_ratio_il12_il13_adj_L[6], 2)), as.character(round(d23_ratio_ifn_il13_adj_L[6], 2)), as.character(round(d23_ratio_il12_il17_adj_L[6], 2)),
                as.character(round(d23_ratio_ifn_il17_adj_L[6], 2)), as.character(round(d23_ratio_il12_il21_adj_L[6], 2)), as.character(round(d23_ratio_ifn_il21_adj_L[6], 2)),
                as.character(round(d23_ratio_pro_il10_adj_L[6], 2)), as.character(round(d23_ratio_th1_il10_adj_L[6], 2)), as.character(round(d23_ratio_th2_il10_adj_L[6], 2)),
                as.character(round(d23_ratio_th17_il10_adj_L[6], 2)), as.character(round(d23_ratio_th1_th2_adj_L[6], 2)), as.character(round(d23_ratio_th1_th17_adj_L[6], 2)))

adjbonpvals12 <- c("Bonferroni P-value",
                   as.character(bonpval(d23_ln_il1_adj_L[6])), as.character(bonpval(d23_ln_il6_adj_L[6])), as.character(bonpval(d23_ln_tnf_adj_L[6])),
                   as.character(bonpval(d23_ln_il12_adj_L[6])), as.character(bonpval(d23_ln_ifn_adj_L[6])),
                   as.character(bonpval(d23_ln_il4_adj_L[6])), as.character(bonpval(d23_ln_il5_adj_L[6])), as.character(bonpval(d23_ln_il13_adj_L[6])),
                   as.character(bonpval(d23_ln_il17_adj_L[6])), as.character(bonpval(d23_ln_il21_adj_L[6])), as.character(bonpval(d23_ln_il10_adj_L[6])),
                   as.character(bonpval(d23_ln_il2_adj_L[6])), as.character(bonpval(d23_ln_gmc_adj_L[6])),
                   as.character(bonpval(d23_ln_igf_adj_L[6])), as.character(bonpval(d23_ratio_il1_il10_adj_L[6])),
                   as.character(bonpval(d23_ratio_il6_il10_adj_L[6])), as.character(bonpval(d23_ratio_tnf_il10_adj_L[6])), as.character(bonpval(d23_ratio_il12_il10_adj_L[6])),
                   as.character(bonpval(d23_ratio_ifn_il10_adj_L[6])), as.character(bonpval(d23_ratio_il4_il10_adj_L[6])), as.character(bonpval(d23_ratio_il5_il10_adj_L[6])),
                   as.character(bonpval(d23_ratio_il13_il10_adj_L[6])), as.character(bonpval(d23_ratio_il17_il10_adj_L[6])), as.character(bonpval(d23_ratio_il21_il10_adj_L[6])),
                   as.character(bonpval(d23_ratio_il2_il10_adj_L[6])), as.character(bonpval(d23_ratio_gmc_il10_adj_L[6])), as.character(bonpval(d23_ratio_il12_il4_adj_L[6])),
                   as.character(bonpval(d23_ratio_ifn_il4_adj_L[6])), as.character(bonpval(d23_ratio_il12_il5_adj_L[6])), as.character(bonpval(d23_ratio_ifn_il5_adj_L[6])),
                   as.character(bonpval(d23_ratio_il12_il13_adj_L[6])), as.character(bonpval(d23_ratio_ifn_il13_adj_L[6])), as.character(bonpval(d23_ratio_il12_il17_adj_L[6])),
                   as.character(bonpval(d23_ratio_ifn_il17_adj_L[6])), as.character(bonpval(d23_ratio_il12_il21_adj_L[6])), as.character(bonpval(d23_ratio_ifn_il21_adj_L[6])),
                   as.character(bonpval(d23_ratio_pro_il10_adj_L[6])), as.character(bonpval(d23_ratio_th1_il10_adj_L[6])), as.character(bonpval(d23_ratio_th2_il10_adj_L[6])),
                   as.character(bonpval(d23_ratio_th17_il10_adj_L[6])), as.character(bonpval(d23_ratio_th1_th2_adj_L[6])), as.character(bonpval(d23_ratio_th1_th17_adj_L[6])))

# Table 10: P-values of change in treatment estimates, unadjusted, and adjusted for multiple testing 
# (by controlling family-wise error rate using the Bonferroni correction). 
tbls12 <- data.table(
  " " = outcomes12,
  "Unadjusted Analysis" = unadjdiffs12, 
  " " = unadjpvals12,
  " " = unadjbonpvals12, 
  "Adjusted Analysis" = adjdiffs12,
  " " = adjpvals12,
  " " = adjbonpvals12
)

write.csv(tbls12, file=here('tables/supplementary/immune_supptable10.csv'))
print(xtable(tbls12), type="html", file=here("tables/supplementary/immune_supptable10.html"))




