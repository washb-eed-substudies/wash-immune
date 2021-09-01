rm(list=ls())
library("xtable")
source(here::here("0-config.R"))

load(here("results/immune_subgroup.RData"))
load(here("results/immune_N_sex_tr_means.RData"))

objects_ind_t2 <- list(il1_t2_N_sex_tr, il6_t2_N_sex_tr, tnf_t2_N_sex_tr, crp_t2_N_sex_tr, il12_t2_N_sex_tr, 
                    ifn_t2_N_sex_tr, il4_t2_N_sex_tr, il5_t2_N_sex_tr, il13_t2_N_sex_tr, il17_t2_N_sex_tr,
                    il21_t2_N_sex_tr, il10_t2_N_sex_tr, il2_t2_N_sex_tr, gmc_t2_N_sex_tr, agp_t2_N_sex_tr, igf_t2_N_sex_tr)
objects_ind_res_t2 <- list(t2_il1_subgroup_L, t2_il6_subgroup_L, t2_tnf_subgroup_L, t2_crp_subgroup_L, t2_il12_subgroup_L, 
                    t2_ifn_subgroup_L, t2_il4_subgroup_L, t2_il5_subgroup_L, t2_il13_subgroup_L, t2_il17_subgroup_L, 
                    t2_il21_subgroup_L, t2_il10_subgroup_L, t2_il2_subgroup_L, t2_gmc_subgroup_L, t2_agp_subgroup_L, t2_igf_subgroup_L)

objects_ind_t3 <- list(il1_t3_N_sex_tr, il6_t3_N_sex_tr, tnf_t3_N_sex_tr, crp_t3_N_sex_tr, il12_t3_N_sex_tr, 
                    ifn_t3_N_sex_tr, il4_t3_N_sex_tr, il5_t3_N_sex_tr, il13_t3_N_sex_tr, il17_t3_N_sex_tr,
                    il21_t3_N_sex_tr, il10_t3_N_sex_tr, il2_t3_N_sex_tr, gmc_t3_N_sex_tr, agp_t3_N_sex_tr, igf_t3_N_sex_tr)
objects_ind_res_t3 <- list(t3_il1_subgroup_L, t3_il6_subgroup_L, t3_tnf_subgroup_L, t3_crp_subgroup_L, t3_il12_subgroup_L, 
                        t3_ifn_subgroup_L, t3_il4_subgroup_L, t3_il5_subgroup_L, t3_il13_subgroup_L, t3_il17_subgroup_L, 
                        t3_il21_subgroup_L, t3_il10_subgroup_L, t3_il2_subgroup_L, t3_gmc_subgroup_L, t3_agp_subgroup_L, t3_igf_subgroup_L)

objects_rat_t2 <- list(t2_ratio_il1_il10_N_sex_tr, t2_ratio_il6_il10_N_sex_tr, t2_ratio_tnf_il10_N_sex_tr, t2_ratio_il12_il10_N_sex_tr, 
                    t2_ratio_ifn_il10_N_sex_tr, t2_ratio_il4_il10_N_sex_tr, t2_ratio_il5_il10_N_sex_tr, t2_ratio_il13_il10_N_sex_tr,
                    t2_ratio_il17_il10_N_sex_tr, t2_ratio_il21_il10_N_sex_tr, t2_ratio_il2_il10_N_sex_tr, t2_ratio_gmc_il10_N_sex_tr,
                    t2_ratio_il12_il4_N_sex_tr, t2_ratio_ifn_il4_N_sex_tr, t2_ratio_il12_il5_N_sex_tr, t2_ratio_ifn_il5_N_sex_tr,
                    t2_ratio_il12_il13_N_sex_tr, t2_ratio_ifn_il13_N_sex_tr, t2_ratio_il12_il17_N_sex_tr, t2_ratio_ifn_il17_N_sex_tr,
                    t2_ratio_il12_il21_N_sex_tr, t2_ratio_ifn_il21_N_sex_tr, t2_ratio_pro_il10_N_sex_tr, t2_ratio_th1_il10_N_sex_tr,
                    t2_ratio_th2_il10_N_sex_tr, t2_ratio_th17_il10_N_sex_tr, t2_ratio_th1_th2_N_sex_tr, t2_ratio_th1_th17_N_sex_tr)
objects_rat_res_t2 <- list(t2_ratio_il1_il10_subgroup_L, t2_ratio_il6_il10_subgroup_L, t2_ratio_tnf_il10_subgroup_L, t2_ratio_il12_il10_subgroup_L,
                           t2_ratio_ifn_il10_subgroup_L, t2_ratio_il4_il10_subgroup_L, t2_ratio_il5_il10_subgroup_L, t2_ratio_il13_il10_subgroup_L,
                           t2_ratio_il17_il10_subgroup_L, t2_ratio_il21_il10_subgroup_L, t2_ratio_il2_il10_subgroup_L, t2_ratio_gmc_il10_subgroup_L,
                           t2_ratio_il12_il4_subgroup_L, t2_ratio_ifn_il4_subgroup_L, t2_ratio_il12_il5_subgroup_L, t2_ratio_ifn_il5_subgroup_L,
                           t2_ratio_il12_il13_subgroup_L, t2_ratio_ifn_il13_subgroup_L, t2_ratio_il12_il17_subgroup_L, t2_ratio_ifn_il17_subgroup_L,
                           t2_ratio_il12_il21_subgroup_L, t2_ratio_ifn_il21_subgroup_L, t2_ratio_pro_il10_subgroup_L, t2_ratio_th1_il10_subgroup_L,
                           t2_ratio_th2_il10_subgroup_L, t2_ratio_th17_il10_subgroup_L, t2_ratio_th1_th2_subgroup_L, t2_ratio_th1_th17_subgroup_L)

objects_rat_t3 <- list(t3_ratio_il1_il10_N_sex_tr, t3_ratio_il6_il10_N_sex_tr, t3_ratio_tnf_il10_N_sex_tr, t3_ratio_il12_il10_N_sex_tr, 
                       t3_ratio_ifn_il10_N_sex_tr, t3_ratio_il4_il10_N_sex_tr, t3_ratio_il5_il10_N_sex_tr, t3_ratio_il13_il10_N_sex_tr,
                       t3_ratio_il17_il10_N_sex_tr, t3_ratio_il21_il10_N_sex_tr, t3_ratio_il2_il10_N_sex_tr, t3_ratio_gmc_il10_N_sex_tr,
                       t3_ratio_il12_il4_N_sex_tr, t3_ratio_ifn_il4_N_sex_tr, t3_ratio_il12_il5_N_sex_tr, t3_ratio_ifn_il5_N_sex_tr,
                       t3_ratio_il12_il13_N_sex_tr, t3_ratio_ifn_il13_N_sex_tr, t3_ratio_il12_il17_N_sex_tr, t3_ratio_ifn_il17_N_sex_tr,
                       t3_ratio_il12_il21_N_sex_tr, t3_ratio_ifn_il21_N_sex_tr, t3_ratio_pro_il10_N_sex_tr, t3_ratio_th1_il10_N_sex_tr,
                       t3_ratio_th2_il10_N_sex_tr, t3_ratio_th17_il10_N_sex_tr, t3_ratio_th1_th2_N_sex_tr, t3_ratio_th1_th17_N_sex_tr)
objects_rat_res_t3 <- list(t3_ratio_il1_il10_subgroup_L, t3_ratio_il6_il10_subgroup_L, t3_ratio_tnf_il10_subgroup_L, t3_ratio_il12_il10_subgroup_L,
                           t3_ratio_ifn_il10_subgroup_L, t3_ratio_il4_il10_subgroup_L, t3_ratio_il5_il10_subgroup_L, t3_ratio_il13_il10_subgroup_L,
                           t3_ratio_il17_il10_subgroup_L, t3_ratio_il21_il10_subgroup_L, t3_ratio_il2_il10_subgroup_L, t3_ratio_gmc_il10_subgroup_L,
                           t3_ratio_il12_il4_subgroup_L, t3_ratio_ifn_il4_subgroup_L, t3_ratio_il12_il5_subgroup_L, t3_ratio_ifn_il5_subgroup_L,
                           t3_ratio_il12_il13_subgroup_L, t3_ratio_ifn_il13_subgroup_L, t3_ratio_il12_il17_subgroup_L, t3_ratio_ifn_il17_subgroup_L,
                           t3_ratio_il12_il21_subgroup_L, t3_ratio_ifn_il21_subgroup_L, t3_ratio_pro_il10_subgroup_L, t3_ratio_th1_il10_subgroup_L,
                           t3_ratio_th2_il10_subgroup_L, t3_ratio_th17_il10_subgroup_L, t3_ratio_th1_th2_subgroup_L, t3_ratio_th1_th17_subgroup_L)


outcome_ind <- c("Outcome, Arm", "Ln IL-1b (pg/ml)", "Control", "Nutrition + WSH", 
             "Ln IL-6 (pg/ml)", "Control", "Nutrition + WSH", 
             "Ln TNF-a (pg/ml)", "Control", "Nutrition + WSH",
             "Ln CRP (mg/L)", "Control", "Nutrition + WSH", 
             "Ln IL-12 (pg/ml)", "Control", "Nutrition + WSH", 
             "Ln IFN-y (pg/ml)", "Control", "Nutrition + WSH", 
             "Ln IL-4 (pg/ml)", "Control", "Nutrition + WSH", 
             "Ln IL-5 (pg/ml)", "Control", "Nutrition + WSH", 
             "Ln IL-13 (pg/ml)", "Control", "Nutrition + WSH", 
             "Ln IL-17A (pg/ml)", "Control", "Nutrition + WSH", 
             "Ln IL-21 (pg/ml)", "Control", "Nutrition + WSH", 
             "Ln IL-10 (pg/ml)", "Control", "Nutrition + WSH", 
             "Ln IL-2 (pg/ml)", "Control", "Nutrition + WSH", 
             "Ln GM-CSF (pg/ml)", "Control", "Nutrition + WSH", 
             "Ln AGP (g/L)", "Control", "Nutrition + WSH",
             "Ln IGF-1 (ug/L)", "Control", "Nutrition + WSH")


tr <- c("Control", "Nutrition + WSH")
outcome_rat <- c( "Outcome", "Ln IL-1b/IL-10", tr,
                 "Ln IL-6/IL-10",  tr,
                 "Ln TNF-a/IL-10", tr, 
                 "Ln IL-12/IL-10", tr,  
                 "Ln IFN-y/IL-10", tr, 
                 "Ln IL-4/IL-10", tr, 
                 "Ln IL-5/IL-10", tr,
                 "Ln IL-13/IL-10", tr,
                 "Ln IL-17A/IL-10", tr,
                 "Ln IL-21/IL-10", tr,
                 "Ln IL-2/IL-10", tr,
                 "Ln GM-CSF/IL-10", tr,
                 "Ln IL-12/IL-4", tr,
                 "Ln IFN-y/IL-4", tr,
                 "Ln IL-12/IL-5", tr,
                 "Ln IFN-y/IL-5", tr,
                 "Ln IL-12/IL-13", tr,
                 "Ln IFN-y/IL-13", tr, 
                 "Ln IL-12/IL-17A", tr,
                 "Ln IFN-y/IL-17A", tr, 
                 "Ln IL-12/IL-21", tr,
                 "Ln IFN-y/IL-21", tr, 
                 "Ln Pro-inflammatory cytokines/IL-10", tr,
                 "Ln Th1/IL-10", tr,
                 "Ln Th2/IL-10", tr,
                 "Ln Th17/IL-10", tr,
                 "Ln Th1/Th2", tr,
                 "Ln Th1/Th17", tr)

tbl <- function(outcomes, object_N_mean, object_res){
  nf <- c("N")
  meanf <- c("Mean")
  sdf <- c("SD")
  nm <- c("N")
  meanm <- c("Mean") 
  sdm <- c("SD")
  resf <- c("Unadjusted difference: Intervention vs.Control (95% CI)")
  resm <- c("Unadjusted difference:Intervention vs.Control (95% CI)")
  pvalf <- c("P-value")
  pvalm <- c("P-value")
  pvalint <- c("")
  
  for (i in 1:length(object_N_mean)){
    sex_tr <- object_N_mean[[i]] %>% filter(tr %in% c("Control", "Nutrition + WSH"))
    sex_tr[,4:5] <- round(sex_tr[,4:5], 2)
    sex_tr_f <- sex_tr[1:2,]
    sex_tr_m <- sex_tr[3:4,]
    res <- object_res[[i]]
    res[,2:8] <- round(res[,2:8], 2)

    nf <- c(nf, "", as.character(sex_tr_f[1, 3]), as.character(sex_tr_f[2, 3]))
    nm <- c(nm, "", as.character(sex_tr_m[1, 3]), as.character(sex_tr_m[2, 3]))
    meanf <- c(meanf, "", as.character(sex_tr_f[1, 4]), as.character(sex_tr_f[2, 4]))
    meanm <- c(meanm, "", as.character(sex_tr_m[1, 4]), as.character(sex_tr_m[2, 4]))
    sdf <- c(sdf, "", as.character(sex_tr_f[1, 5]), as.character(sex_tr_f[2, 5]))
    sdm <- c(sdm, "", as.character(sex_tr_m[1, 5]), as.character(sex_tr_m[2, 5]))
  
    resf <- c(resf, "","",paste0(res[1, 2], " (", res[1, 4], ", ", res[1, 5], ")", sep=""))
    resm <- c(resm, "","",paste0(res[2, 2], " (", res[2, 4], ", ", res[2, 5], ")", sep=""))
  
    pvalf <- c(pvalf,"","",as.character(res[1, 7]))
    pvalm <- c(pvalm,"","",as.character(res[2, 7]))
    
    pvalint <- c(pvalint,"","",as.character(res[1, 8]))
  }
  
  data.table( 
    " " = outcomes,
    "Female" = nf,
    " " = meanf,
    " " = sdf,
    " " =  resf, 
    " " = pvalf,
    "Male" = nm,
    " " = meanm,
    " " = sdm,
    " " = resm,
    " " = pvalm,
    "P-value for interaction" = pvalint
  )
}

tbls9 <- tbl(outcome_ind, objects_ind_t2, objects_ind_res_t2)
tbls10 <- tbl(outcome_rat, objects_rat_t2, objects_rat_res_t2)
tbls11 <- tbl(outcome_ind, objects_ind_t3, objects_ind_res_t3)
tbls12 <- tbl(outcome_rat, objects_rat_t3, objects_rat_res_t3)


write.csv(tbls9, file=here("tables/supplementary/immune_supptable9.csv"))
write.csv(tbls10, file=here("tables/supplementary/immune_supptable10.csv"))
write.csv(tbls11, file=here("tables/supplementary/immune_supptable11.csv"))
write.csv(tbls12, file=here("tables/supplementary/immune_supptable12.csv"))


