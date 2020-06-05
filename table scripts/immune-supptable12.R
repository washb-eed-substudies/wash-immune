rm(list=ls())
library("xtable")
source(here::here("0-config.R"))

load(here("results/immune_subgroup.RData"))
load(here("results/immune_N_sex.RData"))


####Table S12 ####
#Supplemental Table 12: Risk difference in inflammation biomarkers of children in the WASH + Nutrition treatment arm relative to control, stratified by sex

#2 uncertainties to flag:

####I entered N female as [1, 2] and N male as [2, 2] (eg female = as.character(t2_ratio_il1_il10_N_sex[1, 2])). 

#Headings

outcometbl12 <- c( "Outcome", paste("Ln IL-1", "β", "/IL-10", sep=""),  
                 "Ln IL-6/IL-10",  
                 paste("Ln TNF-", "α", "/IL-10", sep=""), 
                 "Ln IL-12/IL-10", 
                 paste("Ln IFN-", "γ", "/IL-10", sep=""), 
                 "Ln IL-4/IL-10", 
                 "Ln IL-5/IL-10", 
                 "Ln IL-13/IL-10", 
                 "Ln IL-17A/IL-10", 
                 "Ln IL-21/IL-10", 
                 "Ln IL-2/IL-10", 
                 "Ln GM-CSF/IL-10", 
                 "Ln IL-12/IL-4", 
                 paste("Ln IFN-", "γ", "/IL-4", sep=""), 
                 "Ln IL-12/IL-5", 
                 paste("Ln IFN-", "γ", "/IL-5", sep=""), 
                 "Ln IL-12/IL-13", 
                 paste("Ln IFN-", "γ", "/IL-13", sep=""), 
                 "Ln IL-12/IL-17A", 
                 paste("Ln IFN-", "γ", "/IL-17A", sep=""), 
                 "Ln IL-12/IL-21",
                 paste("Ln IFN-", "γ", "/IL-21", sep=""), 
                 "Ln Pro-inflammatory cytokines/IL-10", 
                 "Ln Th1/IL-10", 
                 "Ln Th2/IL-10",  
                 "Ln Th17/IL-10", 
                 "Ln Th1/Th2",  
                 "Ln Th1/Th17")

# #Calculate N
# 
#Female

####error here, can't calculate n for each group
as.character(t2_ratio_il1_il10_N_sex[1, 2])
ntbl12f <- c("N", 
             as.character(t2_ratio_il1_il10_N_sex[1, 2]), 
              as.character(t2_ratio_il6_il10_N_sex[1, 2]), 
              as.character(t2_ratio_tnf_il10_N_sex[1, 2]), 
              as.character(t2_ratio_il12_il10_N_sex[1, 2]), 
              as.character(t2_ratio_ifn_il10_N_sex[1, 2]), 
              as.character(t2_ratio_il4_il10_N_sex[1, 2]), 
              as.character(t2_ratio_il5_il10_N_sex[1, 2]), 
              as.character(t2_ratio_il13_il10_N_sex[1, 2]),     
              as.character(t2_ratio_il17_il10_N_sex[1, 2]), 
              as.character(t2_ratio_il21_il10_N_sex[1, 2]), 
              as.character(t2_ratio_il2_il10_N_sex[1, 2]),  
              as.character(t2_ratio_gmc_il10_N_sex[1, 2]), 
              as.character(t2_ratio_il12_il4_N_sex[1, 2]), 
              as.character(t2_ratio_ifn_il4_N_sex[1, 2]), 
              as.character(t2_ratio_il12_il5_N_sex[1, 2]), 
              as.character(t2_ratio_ifn_il5_N_sex[1, 2]), 
              as.character(t2_ratio_il12_il13_N_sex[1, 2]), 
              as.character(t2_ratio_ifn_il13_N_sex[1, 2]), 
              as.character(t2_ratio_il12_il17_N_sex[1, 2]), 
              as.character(t2_ratio_ifn_il17_N_sex[1, 2]), 
              as.character(t2_ratio_il12_il21_N_sex[1, 2]), 
              as.character(t2_ratio_ifn_il21_N_sex[1, 2]), 
              as.character(t2_ratio_pro_il10_N_sex[1, 2]), 
              as.character(t2_ratio_th1_il10_N_sex[1, 2]), 
              as.character(t2_ratio_th2_il10_N_sex[1, 2]), 
              as.character(t2_ratio_th17_il10_N_sex[1, 2]), 
              as.character(t2_ratio_th1_th2_N_sex[1, 2]), 
              as.character(t2_ratio_th1_th17_N_sex[1, 2]))

#male

ntbl12m <- c("N", 
             as.character(t2_ratio_il1_il10_N_sex[2, 2]), 
             as.character(t2_ratio_il6_il10_N_sex[2, 2]), 
             as.character(t2_ratio_tnf_il10_N_sex[2, 2]), 
             as.character(t2_ratio_il12_il10_N_sex[2, 2]), 
             as.character(t2_ratio_ifn_il10_N_sex[2, 2]), 
             as.character(t2_ratio_il4_il10_N_sex[2, 2]), 
             as.character(t2_ratio_il5_il10_N_sex[2, 2]), 
             as.character(t2_ratio_il13_il10_N_sex[2, 2]),     
             as.character(t2_ratio_il17_il10_N_sex[2, 2]), 
             as.character(t2_ratio_il21_il10_N_sex[2, 2]), 
             as.character(t2_ratio_il2_il10_N_sex[2, 2]),  
             as.character(t2_ratio_gmc_il10_N_sex[2, 2]), 
             as.character(t2_ratio_il12_il4_N_sex[2, 2]), 
             as.character(t2_ratio_ifn_il4_N_sex[2, 2]), 
             as.character(t2_ratio_il12_il5_N_sex[2, 2]), 
             as.character(t2_ratio_ifn_il5_N_sex[2, 2]), 
             as.character(t2_ratio_il12_il13_N_sex[2, 2]), 
             as.character(t2_ratio_ifn_il13_N_sex[2, 2]), 
             as.character(t2_ratio_il12_il17_N_sex[2, 2]), 
             as.character(t2_ratio_ifn_il17_N_sex[2, 2]), 
             as.character(t2_ratio_il12_il21_N_sex[2, 2]), 
             as.character(t2_ratio_ifn_il21_N_sex[2, 2]), 
             as.character(t2_ratio_pro_il10_N_sex[2, 2]), 
             as.character(t2_ratio_th1_il10_N_sex[2, 2]), 
             as.character(t2_ratio_th2_il10_N_sex[2, 2]), 
             as.character(t2_ratio_th17_il10_N_sex[2, 2]), 
             as.character(t2_ratio_th1_th2_N_sex[2, 2]), 
             as.character(t2_ratio_th1_th17_N_sex[2, 2]))

#Mean by subgroup

#Female

meantbl12f <- c("Mean", 
             as.character( round (t2_ratio_il1_il10_N_sex$mean[1], 2)), 
             as.character( round (t2_ratio_il6_il10_N_sex$mean[1], 2)), 
             as.character( round (t2_ratio_tnf_il10_N_sex$mean[1], 2)), 
             as.character( round (t2_ratio_il12_il10_N_sex$mean[1], 2)), 
             as.character( round (t2_ratio_ifn_il10_N_sex$mean[1], 2)), 
             as.character( round (t2_ratio_il4_il10_N_sex$mean[1], 2)), 
             as.character( round (t2_ratio_il5_il10_N_sex$mean[1], 2)), 
             as.character( round (t2_ratio_il13_il10_N_sex$mean[1], 2)),     
             as.character( round (t2_ratio_il17_il10_N_sex$mean[1], 2)), 
             as.character( round (t2_ratio_il21_il10_N_sex$mean[1], 2)), 
             as.character( round (t2_ratio_il2_il10_N_sex$mean[1], 2)),  
             as.character( round (t2_ratio_gmc_il10_N_sex$mean[1], 2)), 
             as.character( round (t2_ratio_il12_il4_N_sex$mean[1], 2)), 
             as.character( round (t2_ratio_ifn_il4_N_sex$mean[1], 2)), 
             as.character( round (t2_ratio_il12_il5_N_sex$mean[1], 2)), 
             as.character( round (t2_ratio_ifn_il5_N_sex$mean[1], 2)), 
             as.character( round (t2_ratio_il12_il13_N_sex$mean[1], 2)), 
             as.character( round (t2_ratio_ifn_il13_N_sex$mean[1], 2)), 
             as.character( round (t2_ratio_il12_il17_N_sex$mean[1], 2)), 
             as.character( round (t2_ratio_ifn_il17_N_sex$mean[1], 2)), 
             as.character( round (t2_ratio_il12_il21_N_sex$mean[1], 2)), 
             as.character( round (t2_ratio_ifn_il21_N_sex$mean[1], 2)), 
             as.character( round (t2_ratio_pro_il10_N_sex$mean[1], 2)), 
             as.character( round (t2_ratio_th1_il10_N_sex$mean[1], 2)), 
             as.character( round (t2_ratio_th2_il10_N_sex$mean[1], 2)), 
             as.character( round (t2_ratio_th17_il10_N_sex$mean[1], 2)), 
             as.character( round (t2_ratio_th1_th2_N_sex$mean[1], 2)), 
             as.character( round (t2_ratio_th1_th17_N_sex$mean[1], 2)))

#male

meantbl12m <- c("Mean", 
                as.character( round (t2_ratio_il1_il10_N_sex$mean[2], 2)), 
                as.character( round (t2_ratio_il6_il10_N_sex$mean[2], 2)), 
                as.character( round (t2_ratio_tnf_il10_N_sex$mean[2], 2)), 
                as.character( round (t2_ratio_il12_il10_N_sex$mean[2], 2)), 
                as.character( round (t2_ratio_ifn_il10_N_sex$mean[2], 2)), 
                as.character( round (t2_ratio_il4_il10_N_sex$mean[2], 2)), 
                as.character( round (t2_ratio_il5_il10_N_sex$mean[2], 2)), 
                as.character( round (t2_ratio_il13_il10_N_sex$mean[2], 2)),     
                as.character( round (t2_ratio_il17_il10_N_sex$mean[2], 2)), 
                as.character( round (t2_ratio_il21_il10_N_sex$mean[2], 2)), 
                as.character( round (t2_ratio_il2_il10_N_sex$mean[2], 2)),  
                as.character( round (t2_ratio_gmc_il10_N_sex$mean[2], 2)), 
                as.character( round (t2_ratio_il12_il4_N_sex$mean[2], 2)), 
                as.character( round (t2_ratio_ifn_il4_N_sex$mean[2], 2)), 
                as.character( round (t2_ratio_il12_il5_N_sex$mean[2], 2)), 
                as.character( round (t2_ratio_ifn_il5_N_sex$mean[2], 2)), 
                as.character( round (t2_ratio_il12_il13_N_sex$mean[2], 2)), 
                as.character( round (t2_ratio_ifn_il13_N_sex$mean[2], 2)), 
                as.character( round (t2_ratio_il12_il17_N_sex$mean[2], 2)), 
                as.character( round (t2_ratio_ifn_il17_N_sex$mean[2], 2)), 
                as.character( round (t2_ratio_il12_il21_N_sex$mean[2], 2)), 
                as.character( round (t2_ratio_ifn_il21_N_sex$mean[2], 2)), 
                as.character( round (t2_ratio_pro_il10_N_sex$mean[2], 2)), 
                as.character( round (t2_ratio_th1_il10_N_sex$mean[2], 2)), 
                as.character( round (t2_ratio_th2_il10_N_sex$mean[2], 2)), 
                as.character( round (t2_ratio_th17_il10_N_sex$mean[2], 2)), 
                as.character( round (t2_ratio_th1_th2_N_sex$mean[2], 2)), 
                as.character( round (t2_ratio_th1_th17_N_sex$mean[2], 2)))

#Standard deviation

#Female


sdtbl12f <- c("SD", 
                as.character( round (t2_ratio_il1_il10_N_sex$sd[1], 2)), 
                as.character( round (t2_ratio_il6_il10_N_sex$sd[1], 2)), 
                as.character( round (t2_ratio_tnf_il10_N_sex$sd[1], 2)), 
                as.character( round (t2_ratio_il12_il10_N_sex$sd[1], 2)), 
                as.character( round (t2_ratio_ifn_il10_N_sex$sd[1], 2)), 
                as.character( round (t2_ratio_il4_il10_N_sex$sd[1], 2)), 
                as.character( round (t2_ratio_il5_il10_N_sex$sd[1], 2)), 
                as.character( round (t2_ratio_il13_il10_N_sex$sd[1], 2)),     
                as.character( round (t2_ratio_il17_il10_N_sex$sd[1], 2)), 
                as.character( round (t2_ratio_il21_il10_N_sex$sd[1], 2)), 
                as.character( round (t2_ratio_il2_il10_N_sex$sd[1], 2)),  
                as.character( round (t2_ratio_gmc_il10_N_sex$sd[1], 2)), 
                as.character( round (t2_ratio_il12_il4_N_sex$sd[1], 2)), 
                as.character( round (t2_ratio_ifn_il4_N_sex$sd[1], 2)), 
                as.character( round (t2_ratio_il12_il5_N_sex$sd[1], 2)), 
                as.character( round (t2_ratio_ifn_il5_N_sex$sd[1], 2)), 
                as.character( round (t2_ratio_il12_il13_N_sex$sd[1], 2)), 
                as.character( round (t2_ratio_ifn_il13_N_sex$sd[1], 2)), 
                as.character( round (t2_ratio_il12_il17_N_sex$sd[1], 2)), 
                as.character( round (t2_ratio_ifn_il17_N_sex$sd[1], 2)), 
                as.character( round (t2_ratio_il12_il21_N_sex$sd[1], 2)), 
                as.character( round (t2_ratio_ifn_il21_N_sex$sd[1], 2)), 
                as.character( round (t2_ratio_pro_il10_N_sex$sd[1], 2)), 
                as.character( round (t2_ratio_th1_il10_N_sex$sd[1], 2)), 
                as.character( round (t2_ratio_th2_il10_N_sex$sd[1], 2)), 
                as.character( round (t2_ratio_th17_il10_N_sex$sd[1], 2)), 
                as.character( round (t2_ratio_th1_th2_N_sex$sd[1], 2)), 
                as.character( round (t2_ratio_th1_th17_N_sex$sd[1], 2)))

#male

sdtbl12m <- c("SD", 
                as.character( round (t2_ratio_il1_il10_N_sex$sd[2], 2)), 
                as.character( round (t2_ratio_il6_il10_N_sex$sd[2], 2)), 
                as.character( round (t2_ratio_tnf_il10_N_sex$sd[2], 2)), 
                as.character( round (t2_ratio_il12_il10_N_sex$sd[2], 2)), 
                as.character( round (t2_ratio_ifn_il10_N_sex$sd[2], 2)), 
                as.character( round (t2_ratio_il4_il10_N_sex$sd[2], 2)), 
                as.character( round (t2_ratio_il5_il10_N_sex$sd[2], 2)), 
                as.character( round (t2_ratio_il13_il10_N_sex$sd[2], 2)),     
                as.character( round (t2_ratio_il17_il10_N_sex$sd[2], 2)), 
                as.character( round (t2_ratio_il21_il10_N_sex$sd[2], 2)), 
                as.character( round (t2_ratio_il2_il10_N_sex$sd[2], 2)),  
                as.character( round (t2_ratio_gmc_il10_N_sex$sd[2], 2)), 
                as.character( round (t2_ratio_il12_il4_N_sex$sd[2], 2)), 
                as.character( round (t2_ratio_ifn_il4_N_sex$sd[2], 2)), 
                as.character( round (t2_ratio_il12_il5_N_sex$sd[2], 2)), 
                as.character( round (t2_ratio_ifn_il5_N_sex$sd[2], 2)), 
                as.character( round (t2_ratio_il12_il13_N_sex$sd[2], 2)), 
                as.character( round (t2_ratio_ifn_il13_N_sex$sd[2], 2)), 
                as.character( round (t2_ratio_il12_il17_N_sex$sd[2], 2)), 
                as.character( round (t2_ratio_ifn_il17_N_sex$sd[2], 2)), 
                as.character( round (t2_ratio_il12_il21_N_sex$sd[2], 2)), 
                as.character( round (t2_ratio_ifn_il21_N_sex$sd[2], 2)), 
                as.character( round (t2_ratio_pro_il10_N_sex$sd[2], 2)), 
                as.character( round (t2_ratio_th1_il10_N_sex$sd[2], 2)), 
                as.character( round (t2_ratio_th2_il10_N_sex$sd[2], 2)), 
                as.character( round (t2_ratio_th17_il10_N_sex$sd[2], 2)), 
                as.character( round (t2_ratio_th1_th2_N_sex$sd[2], 2)), 
                as.character( round (t2_ratio_th1_th17_N_sex$sd[2], 2)))

#Risk difference

#female

rdtbl12f <- c(as.character(round(t2_ratio_il1_il10_subgroup_L$RD[1], 2)), 
              as.character(round(t2_ratio_il6_il10_subgroup_L$RD[1], 2)), 
              as.character(round(t2_ratio_tnf_il10_subgroup_L$RD[1], 2)), 
              as.character(round(t2_ratio_il12_il10_subgroup_L$RD[1], 2)), 
              as.character(round(t2_ratio_ifn_il10_subgroup_L$RD[1], 2)), 
              as.character(round(t2_ratio_il4_il10_subgroup_L$RD[1], 2)), 
              as.character(round(t2_ratio_il5_il10_subgroup_L$RD[1], 2)), 
              as.character(round(t2_ratio_il13_il10_subgroup_L$RD[1], 2)),     
              as.character(round(t2_ratio_il17_il10_subgroup_L$RD[1], 2)), 
              as.character(round(t2_ratio_il21_il10_subgroup_L$RD[1], 2)), 
              as.character(round(t2_ratio_il2_il10_subgroup_L$RD[1], 2)),  
              as.character(round(t2_ratio_gmc_il10_subgroup_L$RD[1], 2)), 
              as.character(round(t2_ratio_il12_il4_subgroup_L$RD[1], 2)), 
              as.character(round(t2_ratio_ifn_il4_subgroup_L$RD[1], 2)), 
              as.character(round(t2_ratio_il12_il5_subgroup_L$RD[1], 2)), 
              as.character(round(t2_ratio_ifn_il5_subgroup_L$RD[1], 2)), 
              as.character(round(t2_ratio_il12_il13_subgroup_L$RD[1], 2)), 
              as.character(round(t2_ratio_ifn_il13_subgroup_L$RD[1], 2)), 
              as.character(round(t2_ratio_il12_il17_subgroup_L$RD[1], 2)), 
              as.character(round(t2_ratio_ifn_il17_subgroup_L$RD[1], 2)), 
              as.character(round(t2_ratio_il12_il21_subgroup_L$RD[1], 2)), 
              as.character(round(t2_ratio_ifn_il21_subgroup_L$RD[1], 2)), 
              as.character(round(t2_ratio_pro_il10_subgroup_L$RD[1], 2)), 
              as.character(round(t2_ratio_th1_il10_subgroup_L$RD[1], 2)), 
              as.character(round(t2_ratio_th2_il10_subgroup_L$RD[1], 2)), 
              as.character(round(t2_ratio_th17_il10_subgroup_L$RD[1], 2)), 
              as.character(round(t2_ratio_th1_th2_subgroup_L$RD[1], 2)), 
              as.character(round(t2_ratio_th1_th17_subgroup_L$RD[1], 2)))

#male

rdtbl12m <- c(as.character(round(t2_ratio_il1_il10_subgroup_L$RD[2], 2)), 
               as.character(round(t2_ratio_il6_il10_subgroup_L$RD[2], 2)),  
               as.character(round(t2_ratio_tnf_il10_subgroup_L$RD[2], 2)), 
             as.character(round(t2_ratio_il12_il10_subgroup_L$RD[2], 2)), 
              as.character(round(t2_ratio_ifn_il10_subgroup_L$RD[2], 2)), 
              as.character(round(t2_ratio_il4_il10_subgroup_L$RD[2], 2)), 
              as.character(round(t2_ratio_il5_il10_subgroup_L$RD[2], 2)),  
               as.character(round(t2_ratio_il13_il10_subgroup_L$RD[2], 2)),    
               as.character(round(t2_ratio_il17_il10_subgroup_L$RD[2], 2)),  
               as.character(round(t2_ratio_il21_il10_subgroup_L$RD[2], 2)),  
               as.character(round(t2_ratio_il2_il10_subgroup_L$RD[2], 2)),  
              as.character(round(t2_ratio_gmc_il10_subgroup_L$RD[2], 2)),  
               as.character(round(t2_ratio_il12_il4_subgroup_L$RD[2], 2)),  
              as.character(round(t2_ratio_ifn_il4_subgroup_L$RD[2], 2)),  
              as.character(round(t2_ratio_il12_il5_subgroup_L$RD[2], 2)),  
              as.character(round(t2_ratio_ifn_il5_subgroup_L$RD[2], 2)), 
               as.character(round(t2_ratio_il12_il13_subgroup_L$RD[2], 2)), 
             as.character(round(t2_ratio_ifn_il13_subgroup_L$RD[2], 2)), 
             as.character(round(t2_ratio_il12_il17_subgroup_L$RD[2], 2)),  
             as.character(round(t2_ratio_ifn_il17_subgroup_L$RD[2], 2)),  
             as.character(round(t2_ratio_il12_il21_subgroup_L$RD[2], 2)), 
              as.character(round(t2_ratio_ifn_il21_subgroup_L$RD[2], 2)),
               as.character(round(t2_ratio_pro_il10_subgroup_L$RD[2], 2)),  
               as.character(round(t2_ratio_th1_il10_subgroup_L$RD[2], 2)),  
              as.character(round(t2_ratio_th2_il10_subgroup_L$RD[2], 2)),    
               as.character(round(t2_ratio_th17_il10_subgroup_L$RD[2], 2)), 
             as.character(round(t2_ratio_th1_th2_subgroup_L$RD[2], 2)),  
               as.character(round(t2_ratio_th1_th17_subgroup_L$RD[2], 2)))

#Confidence Interval

#lower bound
#female

cilbtbl12f <- c( as.character(round(t2_ratio_il1_il10_subgroup_L$ci.lb[1], 2)),
                 as.character(round(t2_ratio_il6_il10_subgroup_L$ci.lb[1], 2)), 
                 as.character(round(t2_ratio_tnf_il10_subgroup_L$ci.lb[1], 2)), 
                 as.character(round(t2_ratio_il12_il10_subgroup_L$ci.lb[1], 2)), 
                 as.character(round(t2_ratio_ifn_il10_subgroup_L$ci.lb[1], 2)), 
                 as.character(round(t2_ratio_il4_il10_subgroup_L$ci.lb[1], 2)), 
                 as.character(round(t2_ratio_il5_il10_subgroup_L$ci.lb[1], 2)), 
                 as.character(round(t2_ratio_il13_il10_subgroup_L$ci.lb[1], 2)), 
                 as.character(round(t2_ratio_il17_il10_subgroup_L$ci.lb[1], 2)), 
                 as.character(round(t2_ratio_il21_il10_subgroup_L$ci.lb[1], 2)), 
                 as.character(round(t2_ratio_il2_il10_subgroup_L$ci.lb[1], 2)), 
                 as.character(round(t2_ratio_gmc_il10_subgroup_L$ci.lb[1], 2)), 
                 as.character(round(t2_ratio_il12_il4_subgroup_L$ci.lb[1], 2)), 
                 as.character(round(t2_ratio_ifn_il4_subgroup_L$ci.lb[1], 2)), 
                 as.character(round(t2_ratio_il12_il5_subgroup_L$ci.lb[1], 2)), 
                 as.character(round(t2_ratio_ifn_il5_subgroup_L$ci.lb[1], 2)),
                 as.character(round(t2_ratio_il12_il13_subgroup_L$ci.lb[1], 2)), 
                 as.character(round(t2_ratio_ifn_il13_subgroup_L$ci.lb[1], 2)), 
                 as.character(round(t2_ratio_il12_il17_subgroup_L$ci.lb[1], 2)),
                 as.character(round(t2_ratio_ifn_il17_subgroup_L$ci.lb[1], 2)), 
                 as.character(round(t2_ratio_il12_il21_subgroup_L$ci.lb[1], 2)), 
                 as.character(round(t2_ratio_ifn_il21_subgroup_L$ci.lb[1], 2)), 
                 as.character(round(t2_ratio_pro_il10_subgroup_L$ci.lb[1], 2)), 
                 as.character(round(t2_ratio_th1_il10_subgroup_L$ci.lb[1], 2)), 
                 as.character(round(t2_ratio_th2_il10_subgroup_L$ci.lb[1], 2)), 
                 as.character(round(t2_ratio_th17_il10_subgroup_L$ci.lb[1], 2)), 
                 as.character(round(t2_ratio_th1_th2_subgroup_L$ci.lb[1], 2)), 
                 as.character(round(t2_ratio_th1_th17_subgroup_L$ci.lb[1], 2)))

#male

cilbtbl12m <- c( as.character(round(t2_ratio_il1_il10_subgroup_L$ci.lb[2], 2)), 
                  as.character(round(t2_ratio_il6_il10_subgroup_L$ci.lb[2], 2)),  
                 as.character(round(t2_ratio_tnf_il10_subgroup_L$ci.lb[2], 2)), 
                  as.character(round(t2_ratio_il12_il10_subgroup_L$ci.lb[2], 2)), 
                  as.character(round(t2_ratio_ifn_il10_subgroup_L$ci.lb[2], 2)), 
                  as.character(round(t2_ratio_il4_il10_subgroup_L$ci.lb[2], 2)), 
                  as.character(round(t2_ratio_il5_il10_subgroup_L$ci.lb[2], 2)),  
                 as.character(round(t2_ratio_il13_il10_subgroup_L$ci.lb[2], 2)), 
                  as.character(round(t2_ratio_il17_il10_subgroup_L$ci.lb[2], 2)),  
                  as.character(round(t2_ratio_il21_il10_subgroup_L$ci.lb[2], 2)),  
                 as.character(round(t2_ratio_il2_il10_subgroup_L$ci.lb[2], 2)),  
                 as.character(round(t2_ratio_gmc_il10_subgroup_L$ci.lb[2], 2)),  
                  as.character(round(t2_ratio_il12_il4_subgroup_L$ci.lb[2], 2)),  
                  as.character(round(t2_ratio_ifn_il4_subgroup_L$ci.lb[2], 2)),  
                  as.character(round(t2_ratio_il12_il5_subgroup_L$ci.lb[2], 2)),  
                as.character(round(t2_ratio_ifn_il5_subgroup_L$ci.lb[2], 2)), 
                 as.character(round(t2_ratio_il12_il13_subgroup_L$ci.lb[2], 2)), 
                  as.character(round(t2_ratio_ifn_il13_subgroup_L$ci.lb[2], 2)), 
                  as.character(round(t2_ratio_il12_il17_subgroup_L$ci.lb[2], 2)),  
                  as.character(round(t2_ratio_ifn_il17_subgroup_L$ci.lb[2], 2)),  
                 as.character(round(t2_ratio_il12_il21_subgroup_L$ci.lb[2], 2)), 
                  as.character(round(t2_ratio_ifn_il21_subgroup_L$ci.lb[2], 2)),
                  as.character(round(t2_ratio_pro_il10_subgroup_L$ci.lb[2], 2)),  
                 as.character(round(t2_ratio_th1_il10_subgroup_L$ci.lb[2], 2)),  
                 as.character(round(t2_ratio_th2_il10_subgroup_L$ci.lb[2], 2)),    
                  as.character(round(t2_ratio_th17_il10_subgroup_L$ci.lb[2], 2)), 
                 as.character(round(t2_ratio_th1_th2_subgroup_L$ci.lb[2], 2)),  
                  as.character(round(t2_ratio_th1_th17_subgroup_L$ci.lb[2], 2)))

#Upper bound

#Female
ciubtbl12f <- c( as.character(round(t2_ratio_il1_il10_subgroup_L$ci.ub[1], 2)), 
                 as.character(round(t2_ratio_il6_il10_subgroup_L$ci.ub[1], 2)), 
                 as.character(round(t2_ratio_tnf_il10_subgroup_L$ci.ub[1], 2)), 
                 as.character(round(t2_ratio_il12_il10_subgroup_L$ci.ub[1], 2)), 
                 as.character(round(t2_ratio_ifn_il10_subgroup_L$ci.ub[1], 2)), 
                 as.character(round(t2_ratio_il4_il10_subgroup_L$ci.ub[1], 2)), 
                 as.character(round(t2_ratio_il5_il10_subgroup_L$ci.ub[1], 2)), 
                 as.character(round(t2_ratio_il13_il10_subgroup_L$ci.ub[1], 2)), 
                 as.character(round(t2_ratio_il17_il10_subgroup_L$ci.ub[1], 2)), 
                 as.character(round(t2_ratio_il21_il10_subgroup_L$ci.ub[1], 2)), 
                 as.character(round(t2_ratio_il2_il10_subgroup_L$ci.ub[1], 2)), 
                 as.character(round(t2_ratio_gmc_il10_subgroup_L$ci.ub[1], 2)), 
                 as.character(round(t2_ratio_il12_il4_subgroup_L$ci.ub[1], 2)), 
                 as.character(round(t2_ratio_ifn_il4_subgroup_L$ci.ub[1], 2)),
                 as.character(round(t2_ratio_il12_il5_subgroup_L$ci.ub[1], 2)), 
                 as.character(round(t2_ratio_ifn_il5_subgroup_L$ci.ub[1], 2)), 
                 as.character(round(t2_ratio_il12_il13_subgroup_L$ci.ub[1], 2)),
                 as.character(round(t2_ratio_ifn_il13_subgroup_L$ci.ub[1], 2)), 
                 as.character(round(t2_ratio_il12_il17_subgroup_L$ci.ub[1], 2)), 
                 as.character(round(t2_ratio_ifn_il17_subgroup_L$ci.ub[1], 2)), 
                 as.character(round(t2_ratio_il12_il21_subgroup_L$ci.ub[1], 2)), 
                 as.character(round(t2_ratio_ifn_il21_subgroup_L$ci.ub[1], 2)), 
                 as.character(round(t2_ratio_pro_il10_subgroup_L$ci.ub[1], 2)), 
                 as.character(round(t2_ratio_th1_il10_subgroup_L$ci.ub[1], 2)), 
                 as.character(round(t2_ratio_th2_il10_subgroup_L$ci.ub[1], 2)), 
                 as.character(round(t2_ratio_th17_il10_subgroup_L$ci.ub[1], 2)), 
                 as.character(round(t2_ratio_th1_th2_subgroup_L$ci.ub[1], 2)), 
                 as.character(round(t2_ratio_th1_th17_subgroup_L$ci.ub[1], 2)))

#Male
ciubtbl12m <- c( as.character(round(t2_ratio_il1_il10_subgroup_L$ci.ub[2], 2)), 
                 as.character(round(t2_ratio_il6_il10_subgroup_L$ci.ub[2], 2)),  
                 as.character(round(t2_ratio_tnf_il10_subgroup_L$ci.ub[2], 2)), 
                 as.character(round(t2_ratio_il12_il10_subgroup_L$ci.ub[2], 2)), 
                 as.character(round(t2_ratio_ifn_il10_subgroup_L$ci.ub[2], 2)), 
                 as.character(round(t2_ratio_il4_il10_subgroup_L$ci.ub[2], 2)), 
                 as.character(round(t2_ratio_il5_il10_subgroup_L$ci.ub[2], 2)),  
                 as.character(round(t2_ratio_il13_il10_subgroup_L$ci.ub[2], 2)), 
                 as.character(round(t2_ratio_il17_il10_subgroup_L$ci.ub[2], 2)),  
                 as.character(round(t2_ratio_il21_il10_subgroup_L$ci.ub[2], 2)),  
                as.character(round(t2_ratio_il2_il10_subgroup_L$ci.ub[2], 2)),  
                as.character(round(t2_ratio_gmc_il10_subgroup_L$ci.ub[2], 2)),  
               as.character(round(t2_ratio_il12_il4_subgroup_L$ci.ub[2], 2)),  
                as.character(round(t2_ratio_ifn_il4_subgroup_L$ci.ub[2], 2)),  
                as.character(round(t2_ratio_il12_il5_subgroup_L$ci.ub[2], 2)),  
                as.character(round(t2_ratio_ifn_il5_subgroup_L$ci.ub[2], 2)), 
                as.character(round(t2_ratio_il12_il13_subgroup_L$ci.ub[2], 2)), 
               as.character(round(t2_ratio_ifn_il13_subgroup_L$ci.ub[2], 2)), 
               as.character(round(t2_ratio_il12_il17_subgroup_L$ci.ub[2], 2)),  
               as.character(round(t2_ratio_ifn_il17_subgroup_L$ci.ub[2], 2)),  
               as.character(round(t2_ratio_il12_il21_subgroup_L$ci.ub[2], 2)), 
                as.character(round(t2_ratio_ifn_il21_subgroup_L$ci.ub[2], 2)),
                 as.character(round(t2_ratio_pro_il10_subgroup_L$ci.ub[2], 2)),  
                as.character(round(t2_ratio_th1_il10_subgroup_L$ci.ub[2], 2)),  
                as.character(round(t2_ratio_th2_il10_subgroup_L$ci.ub[2], 2)),    
                as.character(round(t2_ratio_th17_il10_subgroup_L$ci.ub[2], 2)), 
                as.character(round(t2_ratio_th1_th2_subgroup_L$ci.ub[2], 2)),  
                as.character(round(t2_ratio_th1_th17_subgroup_L$ci.ub[2], 2)))


#Intervention Main Effect P-value stratified by sex

#female

pvaltbl12f <- c( "P-value" ,  as.character(round(t2_ratio_il1_il10_subgroup_L$`P-value`[1], 3)), 
                 as.character(round(t2_ratio_il6_il10_subgroup_L$`P-value`[1], 3)),  
                 as.character(round(t2_ratio_tnf_il10_subgroup_L$`P-value`[1], 3)), 
                 as.character(round(t2_ratio_il12_il10_subgroup_L$`P-value`[1], 3)), 
                 as.character(round(t2_ratio_ifn_il10_subgroup_L$`P-value`[1], 3)), 
                 as.character(round(t2_ratio_il4_il10_subgroup_L$`P-value`[1], 3)), 
                 as.character(round(t2_ratio_il5_il10_subgroup_L$`P-value`[1], 3)),  
                 as.character(round(t2_ratio_il13_il10_subgroup_L$`P-value`[1], 3)), 
                  as.character(round(t2_ratio_il17_il10_subgroup_L$`P-value`[1], 3)),  
                 as.character(round(t2_ratio_il21_il10_subgroup_L$`P-value`[1], 3)),  
                  as.character(round(t2_ratio_il2_il10_subgroup_L$`P-value`[1], 3)),  
                  as.character(round(t2_ratio_gmc_il10_subgroup_L$`P-value`[1], 3)),  
                  as.character(round(t2_ratio_il12_il4_subgroup_L$`P-value`[1], 3)),  
                  as.character(round(t2_ratio_ifn_il4_subgroup_L$`P-value`[1], 3)),  
                  as.character(round(t2_ratio_il12_il5_subgroup_L$`P-value`[1], 3)),  
                  as.character(round(t2_ratio_ifn_il5_subgroup_L$`P-value`[1], 3)), 
                  as.character(round(t2_ratio_il12_il13_subgroup_L$`P-value`[1], 3)), 
                  as.character(round(t2_ratio_ifn_il13_subgroup_L$`P-value`[1], 3)), 
                  as.character(round(t2_ratio_il12_il17_subgroup_L$`P-value`[1], 3)),  
                  as.character(round(t2_ratio_ifn_il17_subgroup_L$`P-value`[1], 3)),  
                  as.character(round(t2_ratio_il12_il21_subgroup_L$`P-value`[1], 3)), 
                  as.character(round(t2_ratio_ifn_il21_subgroup_L$`P-value`[1], 3)),
                  as.character(round(t2_ratio_pro_il10_subgroup_L$`P-value`[1], 3)),  
                  as.character(round(t2_ratio_th1_il10_subgroup_L$`P-value`[1], 3)),  
                  as.character(round(t2_ratio_th2_il10_subgroup_L$`P-value`[1], 3)),    
                  as.character(round(t2_ratio_th17_il10_subgroup_L$`P-value`[1], 3)), 
                  as.character(round(t2_ratio_th1_th2_subgroup_L$`P-value`[1], 3)),  
                  as.character(round(t2_ratio_th1_th17_subgroup_L$`P-value`[1], 3)))

#male

pvaltbl12m <- c( "P-value" ,  as.character(round(t2_ratio_il1_il10_subgroup_L$`P-value`[2], 3)), 
                 as.character(round(t2_ratio_il6_il10_subgroup_L$`P-value`[2], 3)),  
                 as.character(round(t2_ratio_tnf_il10_subgroup_L$`P-value`[2], 3)), 
                 as.character(round(t2_ratio_il12_il10_subgroup_L$`P-value`[2], 3)), 
                 as.character(round(t2_ratio_ifn_il10_subgroup_L$`P-value`[2], 3)), 
                 as.character(round(t2_ratio_il4_il10_subgroup_L$`P-value`[2], 3)), 
                 as.character(round(t2_ratio_il5_il10_subgroup_L$`P-value`[2], 3)),  
                 as.character(round(t2_ratio_il13_il10_subgroup_L$`P-value`[2], 3)), 
                 as.character(round(t2_ratio_il17_il10_subgroup_L$`P-value`[2], 3)),  
                 as.character(round(t2_ratio_il21_il10_subgroup_L$`P-value`[2], 3)),  
                 as.character(round(t2_ratio_il2_il10_subgroup_L$`P-value`[2], 3)),  
                 as.character(round(t2_ratio_gmc_il10_subgroup_L$`P-value`[2], 3)),  
                 as.character(round(t2_ratio_il12_il4_subgroup_L$`P-value`[2], 3)),  
                 as.character(round(t2_ratio_ifn_il4_subgroup_L$`P-value`[2], 3)),  
                 as.character(round(t2_ratio_il12_il5_subgroup_L$`P-value`[2], 3)),  
                 as.character(round(t2_ratio_ifn_il5_subgroup_L$`P-value`[2], 3)), 
                 as.character(round(t2_ratio_il12_il13_subgroup_L$`P-value`[2], 3)), 
                 as.character(round(t2_ratio_ifn_il13_subgroup_L$`P-value`[2], 3)), 
                 as.character(round(t2_ratio_il12_il17_subgroup_L$`P-value`[2], 3)),  
                 as.character(round(t2_ratio_ifn_il17_subgroup_L$`P-value`[2], 3)),  
                 as.character(round(t2_ratio_il12_il21_subgroup_L$`P-value`[2], 3)), 
                 as.character(round(t2_ratio_ifn_il21_subgroup_L$`P-value`[2], 3)),
                 as.character(round(t2_ratio_pro_il10_subgroup_L$`P-value`[2], 3)),  
                 as.character(round(t2_ratio_th1_il10_subgroup_L$`P-value`[2], 3)),  
                 as.character(round(t2_ratio_th2_il10_subgroup_L$`P-value`[2], 3)),    
                 as.character(round(t2_ratio_th17_il10_subgroup_L$`P-value`[2], 3)), 
                 as.character(round(t2_ratio_th1_th2_subgroup_L$`P-value`[2], 3)),  
                 as.character(round(t2_ratio_th1_th17_subgroup_L$`P-value`[2], 3)))

#P-value for interaction

pinttbl12 <- c( " " ,  as.character(round(t2_ratio_il1_il10_subgroup_L$int_Pval[1], 3)), 
                 as.character(round(t2_ratio_il6_il10_subgroup_L$int_Pval[1], 3)),  
                 as.character(round(t2_ratio_tnf_il10_subgroup_L$int_Pval[1], 3)), 
                 as.character(round(t2_ratio_il12_il10_subgroup_L$int_Pval[1], 3)), 
                 as.character(round(t2_ratio_ifn_il10_subgroup_L$int_Pval[1], 3)), 
                 as.character(round(t2_ratio_il4_il10_subgroup_L$int_Pval[1], 3)), 
                 as.character(round(t2_ratio_il5_il10_subgroup_L$int_Pval[1], 3)),  
                 as.character(round(t2_ratio_il13_il10_subgroup_L$int_Pval[1], 3)), 
                 as.character(round(t2_ratio_il17_il10_subgroup_L$int_Pval[1], 3)),  
                 as.character(round(t2_ratio_il21_il10_subgroup_L$int_Pval[1], 3)),  
                 as.character(round(t2_ratio_il2_il10_subgroup_L$int_Pval[1], 3)),  
                 as.character(round(t2_ratio_gmc_il10_subgroup_L$int_Pval[1], 3)),  
                 as.character(round(t2_ratio_il12_il4_subgroup_L$int_Pval[1], 3)),  
                 as.character(round(t2_ratio_ifn_il4_subgroup_L$int_Pval[1], 3)),  
                 as.character(round(t2_ratio_il12_il5_subgroup_L$int_Pval[1], 3)),  
                 as.character(round(t2_ratio_ifn_il5_subgroup_L$int_Pval[1], 3)), 
                 as.character(round(t2_ratio_il12_il13_subgroup_L$int_Pval[1], 3)), 
                 as.character(round(t2_ratio_ifn_il13_subgroup_L$int_Pval[1], 3)), 
                 as.character(round(t2_ratio_il12_il17_subgroup_L$int_Pval[1], 3)),  
                 as.character(round(t2_ratio_ifn_il17_subgroup_L$int_Pval[1], 3)),  
                 as.character(round(t2_ratio_il12_il21_subgroup_L$int_Pval[1], 3)), 
                 as.character(round(t2_ratio_ifn_il21_subgroup_L$int_Pval[1], 3)),
                 as.character(round(t2_ratio_pro_il10_subgroup_L$int_Pval[1], 3)),  
                 as.character(round(t2_ratio_th1_il10_subgroup_L$int_Pval[1], 3)),  
                 as.character(round(t2_ratio_th2_il10_subgroup_L$int_Pval[1], 3)),    
                 as.character(round(t2_ratio_th17_il10_subgroup_L$int_Pval[1], 3)), 
                 as.character(round(t2_ratio_th1_th2_subgroup_L$int_Pval[1], 3)),  
                 as.character(round(t2_ratio_th1_th17_subgroup_L$int_Pval[1], 3)))

#Combine RD and CI vectors
RD_formatted_f <- paste0(rdtbl12f, "(", cilbtbl12f, ", ", ciubtbl12f,")")
RD_formatted_f <- c("Unadjusted difference: Intervention vs. Control (95% CI)", RD_formatted_f)
RD_formatted_m <- paste0(rdtbl12m, "(", cilbtbl12m, ", ", ciubtbl12m,")")
RD_formatted_m <- c("Unadjusted difference: Intervention vs. Control (95% CI)", RD_formatted_m)


#Supplemental Table 12: Risk difference in inflammation biomarkers of children in the WASH + Nutrition treatment arm relative to control, stratified by sex

supptbl12 <- data.table( 
  " " = outcometbl12,
  "Female" = ntbl12f,
  " " = meantbl12f,
  " " = sdtbl12f,
  " " =  RD_formatted_f, 
  " " = pvaltbl12f,
  "Male" = ntbl12m,
  " " = meantbl12m,
  " " = RD_formatted_m,
  " " = pvaltbl12m,
  "P-value for interaction" = pinttbl12
)



write.csv(supptbl12, file=here("tables/supplementary/immune_supptable12.csv"))
print(xtable(supptbl12), type="html", file=here("tables/supplementary/immune_supptable12.html"))

