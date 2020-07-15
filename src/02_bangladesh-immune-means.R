
#---------------------------------------
# bangladesh-immune-ages-unadj-analysis.R
#
# audrie lin (audrielin@berkeley.edu)
#
# EE in the WASH B package 
# functions and saving output for the replication 
# compare.R script
# 
# input: 
# bangladesh-dm-ee-anthro-diar-ee-med-plasma-blind-tr-enrol-covariates-lab.csv (from 3-bangladesh-dm-immun-plasma-immun-3.do)
#
# output: 
# immune_N_means.RData (N's and means of immune data)
# immune_unadj_glm.RData (unadj glm of immune data)
#---------------------------------------

#Clear out R environment (remove and loaded data)
rm(list=ls())


######################
###Load in packages
######################

source(here::here("0-config.R"))

######################
###Load in data
######################

#Set working directory to load in blinded treatment assignment and enrolment information
setwd(paste0(dropboxDir,"Data/Cleaned/Audrie/")) #Set working directory

#Load in enrollment data,blinded tr data, stool data for adjusted analysis. Use read.dta() to read the .dta files, or read.csv() to 
#read .csv files. Use stringAsFactors=TRUE so that any character-based variable will be read in as a factor.
lab<-readRDS(paste0(dropboxDir,"Data/Cleaned/Audrie/bangladesh-immune-analysis-dataset.rds"))

table(lab$tr) #crosstab of numbers in each treatment


# re-order the treatment factor for convenience, dropping the arms not included in immune
lab$tr <- factor(lab$tr,levels=c("Control","Nutrition + WSH"))

#calculate N's and mean of immune biomarkers t2 by arm


igf_t2_N_tr<-lab %>%
  subset(t2_ln_igf!="NA") %>%
  group_by (tr) %>%
  do(as.data.frame(washb_mean(Y=.$t2_ln_igf, id=.$block,  print=F)))

igf_t2_N_tr


crp_t2_N_tr<-lab %>%
  subset(t2_ln_crp!="NA") %>%
  group_by (tr) %>%
  do(as.data.frame(washb_mean(Y=.$t2_ln_crp,  id=.$block,  print=F)))

crp_t2_N_tr

agp_t2_N_tr<-lab %>%
  subset(t2_ln_agp!="NA") %>%
  group_by (tr) %>%
  do(as.data.frame(washb_mean(Y=.$t2_ln_agp,  id=.$block,  print=F)))

agp_t2_N_tr

gmc_t2_N_tr<-lab %>%
  subset(t2_ln_gmc!="NA") %>%
  group_by (tr) %>%
  do(as.data.frame(washb_mean(Y=.$t2_ln_gmc, id=.$block,  print=F)))

gmc_t2_N_tr

ifn_t2_N_tr<-lab %>%
  subset(t2_ln_ifn!="NA") %>%
  group_by (tr) %>%
  do(as.data.frame(washb_mean(Y=.$t2_ln_ifn,  id=.$block,  print=F)))

ifn_t2_N_tr

il10_t2_N_tr<-lab %>%
  subset(t2_ln_il10!="NA") %>%
  group_by (tr) %>%
  do(as.data.frame(washb_mean(Y=.$t2_ln_il10,  id=.$block,  print=F)))

il10_t2_N_tr

il12_t2_N_tr<-lab %>%
  subset(t2_ln_il12!="NA") %>%
  group_by (tr) %>%
  do(as.data.frame(washb_mean(Y=.$t2_ln_il12,  id=.$block,  print=F)))

il12_t2_N_tr

il13_t2_N_tr<-lab %>%
  subset(t2_ln_il13!="NA") %>%
  group_by (tr) %>%
  do(as.data.frame(washb_mean(Y=.$t2_ln_il13, id=.$block,  print=F)))

il13_t2_N_tr

il17_t2_N_tr<-lab %>%
  subset(t2_ln_il17!="NA") %>%
  group_by (tr) %>%
  do(as.data.frame(washb_mean(Y=.$t2_ln_il17,  id=.$block,  print=F)))

il17_t2_N_tr

il1_t2_N_tr<-lab %>%
  subset(t2_ln_il1!="NA") %>%
  group_by (tr) %>%
  do(as.data.frame(washb_mean(Y=.$t2_ln_il1, id=.$block,  print=F)))

il1_t2_N_tr

il2_t2_N_tr<-lab %>%
  subset(t2_ln_il2!="NA") %>%
  group_by (tr) %>%
  do(as.data.frame(washb_mean(Y=.$t2_ln_il2,  id=.$block,  print=F)))

il2_t2_N_tr

il21_t2_N_tr<-lab %>%
  subset(t2_ln_il21!="NA") %>%
  group_by (tr) %>%
  do(as.data.frame(washb_mean(Y=.$t2_ln_il21,  id=.$block,  print=F)))

il21_t2_N_tr

il4_t2_N_tr<-lab %>%
  subset(t2_ln_il4!="NA") %>%
  group_by (tr) %>%
  do(as.data.frame(washb_mean(Y=.$t2_ln_il4,  id=.$block,  print=F)))

il4_t2_N_tr

il5_t2_N_tr<-lab %>%
  subset(t2_ln_il5!="NA") %>%
  group_by (tr) %>%
  do(as.data.frame(washb_mean(Y=.$t2_ln_il5, id=.$block,  print=F)))

il5_t2_N_tr

il6_t2_N_tr<-lab %>%
  subset(t2_ln_il6!="NA") %>%
  group_by (tr) %>%
  do(as.data.frame(washb_mean(Y=.$t2_ln_il6, id=.$block,  print=F)))

il6_t2_N_tr

tnf_t2_N_tr<-lab %>%
  subset(t2_ln_tnf!="NA") %>%
  group_by (tr) %>%
  do(as.data.frame(washb_mean(Y=.$t2_ln_tnf,  id=.$block,  print=F)))

tnf_t2_N_tr

t2_ratio_il1_il10_N_tr<-lab %>%
  subset(t2_ratio_il1_il10!="NA") %>%
  group_by (tr) %>%
  do(as.data.frame(washb_mean(Y=.$t2_ratio_il1_il10, id=.$block,  print=F)))

t2_ratio_il1_il10_N_tr

t2_ratio_il6_il10_N_tr<-lab %>%
  subset(t2_ratio_il6_il10!="NA") %>%
  group_by (tr) %>%
  do(as.data.frame(washb_mean(Y=.$t2_ratio_il6_il10,  id=.$block,  print=F)))

t2_ratio_il6_il10_N_tr

t2_ratio_tnf_il10_N_tr<-lab %>%
  subset(t2_ratio_tnf_il10!="NA") %>%
  group_by (tr) %>%
  do(as.data.frame(washb_mean(Y=.$t2_ratio_tnf_il10,  id=.$block,  print=F)))

t2_ratio_tnf_il10_N_tr

t2_ratio_il12_il10_N_tr<-lab %>%
  subset(t2_ratio_il12_il10!="NA") %>%
  group_by (tr) %>%
  do(as.data.frame(washb_mean(Y=.$t2_ratio_il12_il10,  id=.$block,  print=F)))

t2_ratio_il12_il10_N_tr

t2_ratio_ifn_il10_N_tr<-lab %>%
  subset(t2_ratio_ifn_il10!="NA") %>%
  group_by (tr) %>%
  do(as.data.frame(washb_mean(Y=.$t2_ratio_ifn_il10, id=.$block,  print=F)))

t2_ratio_ifn_il10_N_tr

t2_ratio_il4_il10_N_tr<-lab %>%
  subset(t2_ratio_il4_il10!="NA") %>%
  group_by (tr) %>%
  do(as.data.frame(washb_mean(Y=.$t2_ratio_il4_il10,  id=.$block,  print=F)))

t2_ratio_il4_il10_N_tr

t2_ratio_il5_il10_N_tr<-lab %>%
  subset(t2_ratio_il5_il10!="NA") %>%
  group_by (tr) %>%
  do(as.data.frame(washb_mean(Y=.$t2_ratio_il5_il10, id=.$block,  print=F)))

t2_ratio_il5_il10_N_tr

t2_ratio_il13_il10_N_tr<-lab %>%
  subset(t2_ratio_il13_il10!="NA") %>%
  group_by (tr) %>%
  do(as.data.frame(washb_mean(Y=.$t2_ratio_il13_il10, id=.$block,  print=F)))

t2_ratio_il13_il10_N_tr

t2_ratio_il17_il10_N_tr<-lab %>%
  subset(t2_ratio_il17_il10!="NA") %>%
  group_by (tr) %>%
  do(as.data.frame(washb_mean(Y=.$t2_ratio_il17_il10,  id=.$block,  print=F)))

t2_ratio_il17_il10_N_tr

t2_ratio_il21_il10_N_tr<-lab %>%
  subset(t2_ratio_il21_il10!="NA") %>%
  group_by (tr) %>%
  do(as.data.frame(washb_mean(Y=.$t2_ratio_il21_il10, id=.$block,  print=F)))

t2_ratio_il21_il10_N_tr

t2_ratio_il2_il10_N_tr<-lab %>%
  subset(t2_ratio_il2_il10!="NA") %>%
  group_by (tr) %>%
  do(as.data.frame(washb_mean(Y=.$t2_ratio_il2_il10,  id=.$block,  print=F)))

t2_ratio_il2_il10_N_tr

t2_ratio_gmc_il10_N_tr<-lab %>%
  subset(t2_ratio_gmc_il10!="NA") %>%
  group_by (tr) %>%
  do(as.data.frame(washb_mean(Y=.$t2_ratio_gmc_il10,  id=.$block,  print=F)))

t2_ratio_gmc_il10_N_tr

t2_ratio_il12_il4_N_tr<-lab %>%
  subset(t2_ratio_il12_il4!="NA") %>%
  group_by (tr) %>%
  do(as.data.frame(washb_mean(Y=.$t2_ratio_il12_il4,  id=.$block,  print=F)))

t2_ratio_il12_il4_N_tr

t2_ratio_ifn_il4_N_tr<-lab %>%
  subset(t2_ratio_ifn_il4!="NA") %>%
  group_by (tr) %>%
  do(as.data.frame(washb_mean(Y=.$t2_ratio_ifn_il4,  id=.$block,  print=F)))

t2_ratio_ifn_il4_N_tr

t2_ratio_il12_il5_N_tr<-lab %>%
  subset(t2_ratio_il12_il5!="NA") %>%
  group_by (tr) %>%
  do(as.data.frame(washb_mean(Y=.$t2_ratio_il12_il5, id=.$block,  print=F)))

t2_ratio_il12_il5_N_tr

t2_ratio_ifn_il5_N_tr<-lab %>%
  subset(t2_ratio_ifn_il5!="NA") %>%
  group_by (tr) %>%
  do(as.data.frame(washb_mean(Y=.$t2_ratio_ifn_il5, id=.$block,  print=F)))

t2_ratio_ifn_il5_N_tr

t2_ratio_il12_il13_N_tr<-lab %>%
  subset(t2_ratio_il12_il13!="NA") %>%
  group_by (tr) %>%
  do(as.data.frame(washb_mean(Y=.$t2_ratio_il12_il13,  id=.$block,  print=F)))

t2_ratio_il12_il13_N_tr

t2_ratio_ifn_il13_N_tr<-lab %>%
  subset(t2_ratio_ifn_il13!="NA") %>%
  group_by (tr) %>%
  do(as.data.frame(washb_mean(Y=.$t2_ratio_ifn_il13, id=.$block,  print=F)))

t2_ratio_ifn_il13_N_tr

t2_ratio_il12_il17_N_tr<-lab %>%
  subset(t2_ratio_il12_il17!="NA") %>%
  group_by (tr) %>%
  do(as.data.frame(washb_mean(Y=.$t2_ratio_il12_il17,  id=.$block,  print=F)))

t2_ratio_il12_il17_N_tr

t2_ratio_ifn_il17_N_tr<-lab %>%
  subset(t2_ratio_ifn_il17!="NA") %>%
  group_by (tr) %>%
  do(as.data.frame(washb_mean(Y=.$t2_ratio_ifn_il17, id=.$block,  print=F)))

t2_ratio_ifn_il17_N_tr

t2_ratio_il12_il21_N_tr<-lab %>%
  subset(t2_ratio_il12_il21!="NA") %>%
  group_by (tr) %>%
  do(as.data.frame(washb_mean(Y=.$t2_ratio_il12_il21,  id=.$block,  print=F)))

t2_ratio_il12_il21_N_tr

t2_ratio_ifn_il21_N_tr<-lab %>%
  subset(t2_ratio_ifn_il21!="NA") %>%
  group_by (tr) %>%
  do(as.data.frame(washb_mean(Y=.$t2_ratio_ifn_il21,  id=.$block,  print=F)))

t2_ratio_ifn_il21_N_tr

t2_ratio_pro_il10_N_tr<-lab %>%
  subset(t2_ratio_pro_il10!="NA") %>%
  group_by (tr) %>%
  do(as.data.frame(washb_mean(Y=.$t2_ratio_pro_il10, id=.$block,  print=F)))

t2_ratio_pro_il10_N_tr

t2_ratio_th1_il10_N_tr<-lab %>%
  subset(t2_ratio_th1_il10!="NA") %>%
  group_by (tr) %>%
  do(as.data.frame(washb_mean(Y=.$t2_ratio_th1_il10,  id=.$block,  print=F)))

t2_ratio_th1_il10_N_tr

t2_ratio_th2_il10_N_tr<-lab %>%
  subset(t2_ratio_th2_il10!="NA") %>%
  group_by (tr) %>%
  do(as.data.frame(washb_mean(Y=.$t2_ratio_th2_il10,  id=.$block,  print=F)))

t2_ratio_th2_il10_N_tr

t2_ratio_th17_il10_N_tr<-lab %>%
  subset(t2_ratio_th17_il10!="NA") %>%
  group_by (tr) %>%
  do(as.data.frame(washb_mean(Y=.$t2_ratio_th17_il10,  id=.$block,  print=F)))

t2_ratio_th17_il10_N_tr

t2_ratio_th1_th2_N_tr<-lab %>%
  subset(t2_ratio_th1_th2!="NA") %>%
  group_by (tr) %>%
  do(as.data.frame(washb_mean(Y=.$t2_ratio_th1_th2, id=.$block,  print=F)))

t2_ratio_th1_th2_N_tr

t2_ratio_th1_th17_N_tr<-lab %>%
  subset(t2_ratio_th1_th17!="NA") %>%
  group_by (tr) %>%
  do(as.data.frame(washb_mean(Y=.$t2_ratio_th1_th17, id=.$block,  print=F)))

t2_ratio_th1_th17_N_tr



#calculate N's and mean of biomarkers at t3 by arm
igf_t3_N_tr<-lab %>%
  subset(t3_ln_igf!="NA") %>%
  group_by (tr) %>%
  do(as.data.frame(washb_mean(Y=.$t3_ln_igf, id=.$block,  print=F)))

igf_t3_N_tr

gmc_t3_N_tr<-lab %>%
  subset(t3_ln_gmc!="NA") %>%
  group_by (tr) %>%
  do(as.data.frame(washb_mean(Y=.$t3_ln_gmc,id=.$block,  print=F)))

gmc_t3_N_tr

ifn_t3_N_tr<-lab %>%
  subset(t3_ln_ifn!="NA") %>%
  group_by (tr) %>%
  do(as.data.frame(washb_mean(Y=.$t3_ln_ifn,  id=.$block,  print=F)))

ifn_t3_N_tr

il10_t3_N_tr<-lab %>%
  subset(t3_ln_il10!="NA") %>%
  group_by (tr) %>%
  do(as.data.frame(washb_mean(Y=.$t3_ln_il10,  id=.$block,  print=F)))

il10_t3_N_tr

il12_t3_N_tr<-lab %>%
  subset(t3_ln_il12!="NA") %>%
  group_by (tr) %>%
  do(as.data.frame(washb_mean(Y=.$t3_ln_il12, id=.$block,  print=F)))

il12_t3_N_tr

il13_t3_N_tr<-lab %>%
  subset(t3_ln_il13!="NA") %>%
  group_by (tr) %>%
  do(as.data.frame(washb_mean(Y=.$t3_ln_il13,  id=.$block,  print=F)))

il13_t3_N_tr

il17_t3_N_tr<-lab %>%
  subset(t3_ln_il17!="NA") %>%
  group_by (tr) %>%
  do(as.data.frame(washb_mean(Y=.$t3_ln_il17,  id=.$block,  print=F)))

il17_t3_N_tr

il1_t3_N_tr<-lab %>%
  subset(t3_ln_il1!="NA") %>%
  group_by (tr) %>%
  do(as.data.frame(washb_mean(Y=.$t3_ln_il1, id=.$block,  print=F)))

il1_t3_N_tr

il2_t3_N_tr<-lab %>%
  subset(t3_ln_il2!="NA") %>%
  group_by (tr) %>%
  do(as.data.frame(washb_mean(Y=.$t3_ln_il2, id=.$block,  print=F)))

il2_t3_N_tr

il21_t3_N_tr<-lab %>%
  subset(t3_ln_il21!="NA") %>%
  group_by (tr) %>%
  do(as.data.frame(washb_mean(Y=.$t3_ln_il21, id=.$block,  print=F)))

il21_t3_N_tr

il4_t3_N_tr<-lab %>%
  subset(t3_ln_il4!="NA") %>%
  group_by (tr) %>%
  do(as.data.frame(washb_mean(Y=.$t3_ln_il4, id=.$block,  print=F)))

il4_t3_N_tr

il5_t3_N_tr<-lab %>%
  subset(t3_ln_il5!="NA") %>%
  group_by (tr) %>%
  do(as.data.frame(washb_mean(Y=.$t3_ln_il5,  id=.$block,  print=F)))

il5_t3_N_tr

il6_t3_N_tr<-lab %>%
  subset(t3_ln_il6!="NA") %>%
  group_by (tr) %>%
  do(as.data.frame(washb_mean(Y=.$t3_ln_il6,  id=.$block,  print=F)))

il6_t3_N_tr

tnf_t3_N_tr<-lab %>%
  subset(t3_ln_tnf!="NA") %>%
  group_by (tr) %>%
  do(as.data.frame(washb_mean(Y=.$t3_ln_tnf, id=.$block,  print=F)))

tnf_t3_N_tr

t3_ratio_il1_il10_N_tr<-lab %>%
  subset(t3_ratio_il1_il10!="NA") %>%
  group_by (tr) %>%
  do(as.data.frame(washb_mean(Y=.$t3_ratio_il1_il10,  id=.$block,  print=F)))

t3_ratio_il1_il10_N_tr

t3_ratio_il6_il10_N_tr<-lab %>%
  subset(t3_ratio_il6_il10!="NA") %>%
  group_by (tr) %>%
  do(as.data.frame(washb_mean(Y=.$t3_ratio_il6_il10,  id=.$block,  print=F)))

t3_ratio_il6_il10_N_tr

t3_ratio_tnf_il10_N_tr<-lab %>%
  subset(t3_ratio_tnf_il10!="NA") %>%
  group_by (tr) %>%
  do(as.data.frame(washb_mean(Y=.$t3_ratio_tnf_il10,  id=.$block,  print=F)))

t3_ratio_tnf_il10_N_tr

t3_ratio_il12_il10_N_tr<-lab %>%
  subset(t3_ratio_il12_il10!="NA") %>%
  group_by (tr) %>%
  do(as.data.frame(washb_mean(Y=.$t3_ratio_il12_il10,  id=.$block,  print=F)))

t3_ratio_il12_il10_N_tr

t3_ratio_ifn_il10_N_tr<-lab %>%
  subset(t3_ratio_ifn_il10!="NA") %>%
  group_by (tr) %>%
  do(as.data.frame(washb_mean(Y=.$t3_ratio_ifn_il10,  id=.$block,  print=F)))

t3_ratio_ifn_il10_N_tr

t3_ratio_il4_il10_N_tr<-lab %>%
  subset(t3_ratio_il4_il10!="NA") %>%
  group_by (tr) %>%
  do(as.data.frame(washb_mean(Y=.$t3_ratio_il4_il10,  id=.$block,  print=F)))

t3_ratio_il4_il10_N_tr

t3_ratio_il5_il10_N_tr<-lab %>%
  subset(t3_ratio_il5_il10!="NA") %>%
  group_by (tr) %>%
  do(as.data.frame(washb_mean(Y=.$t3_ratio_il5_il10,  id=.$block,  print=F)))

t3_ratio_il5_il10_N_tr

t3_ratio_il13_il10_N_tr<-lab %>%
  subset(t3_ratio_il13_il10!="NA") %>%
  group_by (tr) %>%
  do(as.data.frame(washb_mean(Y=.$t3_ratio_il13_il10,  id=.$block,  print=F)))

t3_ratio_il13_il10_N_tr

t3_ratio_il17_il10_N_tr<-lab %>%
  subset(t3_ratio_il17_il10!="NA") %>%
  group_by (tr) %>%
  do(as.data.frame(washb_mean(Y=.$t3_ratio_il17_il10, id=.$block,  print=F)))

t3_ratio_il17_il10_N_tr

t3_ratio_il21_il10_N_tr<-lab %>%
  subset(t3_ratio_il21_il10!="NA") %>%
  group_by (tr) %>%
  do(as.data.frame(washb_mean(Y=.$t3_ratio_il21_il10, id=.$block,  print=F)))

t3_ratio_il21_il10_N_tr

t3_ratio_il2_il10_N_tr<-lab %>%
  subset(t3_ratio_il2_il10!="NA") %>%
  group_by (tr) %>%
  do(as.data.frame(washb_mean(Y=.$t3_ratio_il2_il10, id=.$block,  print=F)))

t3_ratio_il2_il10_N_tr

t3_ratio_gmc_il10_N_tr<-lab %>%
  subset(t3_ratio_gmc_il10!="NA") %>%
  group_by (tr) %>%
  do(as.data.frame(washb_mean(Y=.$t3_ratio_gmc_il10,  id=.$block,  print=F)))

t3_ratio_gmc_il10_N_tr

t3_ratio_il12_il4_N_tr<-lab %>%
  subset(t3_ratio_il12_il4!="NA") %>%
  group_by (tr) %>%
  do(as.data.frame(washb_mean(Y=.$t3_ratio_il12_il4,  id=.$block,  print=F)))

t3_ratio_il12_il4_N_tr

t3_ratio_ifn_il4_N_tr<-lab %>%
  subset(t3_ratio_ifn_il4!="NA") %>%
  group_by (tr) %>%
  do(as.data.frame(washb_mean(Y=.$t3_ratio_ifn_il4,  id=.$block,  print=F)))

t3_ratio_ifn_il4_N_tr

t3_ratio_il12_il5_N_tr<-lab %>%
  subset(t3_ratio_il12_il5!="NA") %>%
  group_by (tr) %>%
  do(as.data.frame(washb_mean(Y=.$t3_ratio_il12_il5,  id=.$block,  print=F)))

t3_ratio_il12_il5_N_tr

t3_ratio_ifn_il5_N_tr<-lab %>%
  subset(t3_ratio_ifn_il5!="NA") %>%
  group_by (tr) %>%
  do(as.data.frame(washb_mean(Y=.$t3_ratio_ifn_il5, id=.$block,  print=F)))

t3_ratio_ifn_il5_N_tr

t3_ratio_il12_il13_N_tr<-lab %>%
  subset(t3_ratio_il12_il13!="NA") %>%
  group_by (tr) %>%
  do(as.data.frame(washb_mean(Y=.$t3_ratio_il12_il13,  id=.$block,  print=F)))

t3_ratio_il12_il13_N_tr

t3_ratio_ifn_il13_N_tr<-lab %>%
  subset(t3_ratio_ifn_il13!="NA") %>%
  group_by (tr) %>%
  do(as.data.frame(washb_mean(Y=.$t3_ratio_ifn_il13, id=.$block,  print=F)))

t3_ratio_ifn_il13_N_tr

t3_ratio_il12_il17_N_tr<-lab %>%
  subset(t3_ratio_il12_il17!="NA") %>%
  group_by (tr) %>%
  do(as.data.frame(washb_mean(Y=.$t3_ratio_il12_il17,  id=.$block,  print=F)))

t3_ratio_il12_il17_N_tr

t3_ratio_ifn_il17_N_tr<-lab %>%
  subset(t3_ratio_ifn_il17!="NA") %>%
  group_by (tr) %>%
  do(as.data.frame(washb_mean(Y=.$t3_ratio_ifn_il17,  id=.$block,  print=F)))

t3_ratio_ifn_il17_N_tr

t3_ratio_il12_il21_N_tr<-lab %>%
  subset(t3_ratio_il12_il21!="NA") %>%
  group_by (tr) %>%
  do(as.data.frame(washb_mean(Y=.$t3_ratio_il12_il21,  id=.$block,  print=F)))

t3_ratio_il12_il21_N_tr

t3_ratio_ifn_il21_N_tr<-lab %>%
  subset(t3_ratio_ifn_il21!="NA") %>%
  group_by (tr) %>%
  do(as.data.frame(washb_mean(Y=.$t3_ratio_ifn_il21, id=.$block,  print=F)))

t3_ratio_ifn_il21_N_tr

t3_ratio_pro_il10_N_tr<-lab %>%
  subset(t3_ratio_pro_il10!="NA") %>%
  group_by (tr) %>%
  do(as.data.frame(washb_mean(Y=.$t3_ratio_pro_il10,  id=.$block,  print=F)))

t3_ratio_pro_il10_N_tr

t3_ratio_th1_il10_N_tr<-lab %>%
  subset(t3_ratio_th1_il10!="NA") %>%
  group_by (tr) %>%
  do(as.data.frame(washb_mean(Y=.$t3_ratio_th1_il10,  id=.$block,  print=F)))

t3_ratio_th1_il10_N_tr

t3_ratio_th2_il10_N_tr<-lab %>%
  subset(t3_ratio_th2_il10!="NA") %>%
  group_by (tr) %>%
  do(as.data.frame(washb_mean(Y=.$t3_ratio_th2_il10,  id=.$block,  print=F)))

t3_ratio_th2_il10_N_tr

t3_ratio_th17_il10_N_tr<-lab %>%
  subset(t3_ratio_th17_il10!="NA") %>%
  group_by (tr) %>%
  do(as.data.frame(washb_mean(Y=.$t3_ratio_th17_il10,  id=.$block,  print=F)))

t3_ratio_th17_il10_N_tr

t3_ratio_th1_th2_N_tr<-lab %>%
  subset(t3_ratio_th1_th2!="NA") %>%
  group_by (tr) %>%
  do(as.data.frame(washb_mean(Y=.$t3_ratio_th1_th2, id=.$block,  print=F)))

t3_ratio_th1_th2_N_tr

t3_ratio_th1_th17_N_tr<-lab %>%
  subset(t3_ratio_th1_th17!="NA") %>%
  group_by (tr) %>%
  do(as.data.frame(washb_mean(Y=.$t3_ratio_th1_th17,  id=.$block,  print=F)))

t3_ratio_th1_th17_N_tr  

# calculating N, mean, sd by arm between t2 and t3
d23_ln_il1_N_tr<-lab %>%
  subset(d23_ln_il1!="NA") %>%
  group_by (tr) %>%
  do(as.data.frame(washb_mean(Y=.$d23_ln_il1,  id=.$block,  print=F)))

d23_ln_il1_N_tr  

d23_ln_il6_N_tr<-lab %>%
  subset(d23_ln_il6!="NA") %>%
  group_by (tr) %>%
  do(as.data.frame(washb_mean(Y=.$d23_ln_il6,  id=.$block,  print=F)))

d23_ln_il6_N_tr  

d23_ln_tnf_N_tr<-lab %>%
  subset(d23_ln_tnf!="NA") %>%
  group_by (tr) %>%
  do(as.data.frame(washb_mean(Y=.$d23_ln_tnf,  id=.$block,  print=F)))

d23_ln_tnf_N_tr  

d23_ln_il12_N_tr<-lab %>%
  subset(d23_ln_il12!="NA") %>%
  group_by (tr) %>%
  do(as.data.frame(washb_mean(Y=.$d23_ln_il12,  id=.$block,  print=F)))

d23_ln_il12_N_tr  

d23_ln_ifn_N_tr<-lab %>%
  subset(d23_ln_ifn!="NA") %>%
  group_by (tr) %>%
  do(as.data.frame(washb_mean(Y=.$d23_ln_ifn,  id=.$block,  print=F)))

d23_ln_ifn_N_tr  

d23_ln_il4_N_tr<-lab %>%
  subset(d23_ln_il4!="NA") %>%
  group_by (tr) %>%
  do(as.data.frame(washb_mean(Y=.$d23_ln_il4,  id=.$block,  print=F)))

d23_ln_il4_N_tr  

d23_ln_il5_N_tr<-lab %>%
  subset(d23_ln_il5!="NA") %>%
  group_by (tr) %>%
  do(as.data.frame(washb_mean(Y=.$d23_ln_il5, id=.$block,  print=F)))

d23_ln_il5_N_tr  

d23_ln_il13_N_tr<-lab %>%
  subset(d23_ln_il13!="NA") %>%
  group_by (tr) %>%
  do(as.data.frame(washb_mean(Y=.$d23_ln_il13,  id=.$block,  print=F)))

d23_ln_il13_N_tr  

d23_ln_il17_N_tr<-lab %>%
  subset(d23_ln_il17!="NA") %>%
  group_by (tr) %>%
  do(as.data.frame(washb_mean(Y=.$d23_ln_il17,id=.$block,  print=F)))

d23_ln_il17_N_tr

d23_ln_il21_N_tr<-lab %>%
  subset(d23_ln_il21!="NA") %>%
  group_by (tr) %>%
  do(as.data.frame(washb_mean(Y=.$d23_ln_il21,  id=.$block,  print=F)))

d23_ln_il21_N_tr  

d23_ln_il10_N_tr<-lab %>%
  subset(d23_ln_il10!="NA") %>%
  group_by (tr) %>%
  do(as.data.frame(washb_mean(Y=.$d23_ln_il10,  id=.$block,  print=F)))

d23_ln_il10_N_tr  

d23_ln_il2_N_tr<-lab %>%
  subset(d23_ln_il2!="NA") %>%
  group_by (tr) %>%
  do(as.data.frame(washb_mean(Y=.$d23_ln_il2,  id=.$block,  print=F)))

d23_ln_il2_N_tr  

d23_ln_gmc_N_tr<-lab %>%
  subset(d23_ln_gmc!="NA") %>%
  group_by (tr) %>%
  do(as.data.frame(washb_mean(Y=.$d23_ln_gmc, id=.$block,  print=F)))

d23_ln_gmc_N_tr  

d23_ln_igf_N_tr<-lab %>%
  subset(d23_ln_igf!="NA") %>%
  group_by (tr) %>%
  do(as.data.frame(washb_mean(Y=.$d23_ln_igf, id=.$block,  print=F)))

d23_ln_igf_N_tr  

d23_ratio_il1_il10_N_tr<-lab %>%
  subset(d23_ratio_il1_il10!="NA") %>%
  group_by (tr) %>%
  do(as.data.frame(washb_mean(Y=.$d23_ratio_il1_il10,id=.$block,  print=F)))

d23_ratio_il1_il10_N_tr

d23_ratio_il6_il10_N_tr<-lab %>%
  subset(d23_ratio_il6_il10!="NA") %>%
  group_by (tr) %>%
  do(as.data.frame(washb_mean(Y=.$d23_ratio_il6_il10, id=.$block,  print=F)))

d23_ratio_il6_il10_N_tr

d23_ratio_tnf_il10_N_tr<-lab %>%
  subset(d23_ratio_tnf_il10!="NA") %>%
  group_by (tr) %>%
  do(as.data.frame(washb_mean(Y=.$d23_ratio_tnf_il10, id=.$block,  print=F)))

d23_ratio_tnf_il10_N_tr

d23_ratio_il12_il10_N_tr<-lab %>%
  subset(d23_ratio_il12_il10!="NA") %>%
  group_by (tr) %>%
  do(as.data.frame(washb_mean(Y=.$d23_ratio_il12_il10,  id=.$block,  print=F)))

d23_ratio_il12_il10_N_tr

d23_ratio_ifn_il10_N_tr<-lab %>%
  subset(d23_ratio_ifn_il10!="NA") %>%
  group_by (tr) %>%
  do(as.data.frame(washb_mean(Y=.$d23_ratio_ifn_il10, id=.$block,  print=F)))

d23_ratio_ifn_il10_N_tr

d23_ratio_il4_il10_N_tr<-lab %>%
  subset(d23_ratio_il4_il10!="NA") %>%
  group_by (tr) %>%
  do(as.data.frame(washb_mean(Y=.$d23_ratio_il4_il10,  id=.$block,  print=F)))

d23_ratio_il4_il10_N_tr

d23_ratio_il5_il10_N_tr<-lab %>%
  subset(d23_ratio_il5_il10!="NA") %>%
  group_by (tr) %>%
  do(as.data.frame(washb_mean(Y=.$d23_ratio_il5_il10, id=.$block,  print=F)))

d23_ratio_il5_il10_N_tr

d23_ratio_il13_il10_N_tr<-lab %>%
  subset(d23_ratio_il13_il10!="NA") %>%
  group_by (tr) %>%
  do(as.data.frame(washb_mean(Y=.$d23_ratio_il13_il10,  id=.$block,  print=F)))

d23_ratio_il13_il10_N_tr

d23_ratio_il17_il10_N_tr<-lab %>%
  subset(d23_ratio_il17_il10!="NA") %>%
  group_by (tr) %>%
  do(as.data.frame(washb_mean(Y=.$d23_ratio_il17_il10,  id=.$block,  print=F)))

d23_ratio_il17_il10_N_tr

d23_ratio_il21_il10_N_tr<-lab %>%
  subset(d23_ratio_il21_il10!="NA") %>%
  group_by (tr) %>%
  do(as.data.frame(washb_mean(Y=.$d23_ratio_il21_il10,  id=.$block,  print=F)))

d23_ratio_il21_il10_N_tr

d23_ratio_il2_il10_N_tr<-lab %>%
  subset(d23_ratio_il2_il10!="NA") %>%
  group_by (tr) %>%
  do(as.data.frame(washb_mean(Y=.$d23_ratio_il2_il10,  id=.$block,  print=F)))

d23_ratio_il2_il10_N_tr

d23_ratio_gmc_il10_N_tr<-lab %>%
  subset(d23_ratio_gmc_il10!="NA") %>%
  group_by (tr) %>%
  do(as.data.frame(washb_mean(Y=.$d23_ratio_gmc_il10,  id=.$block,  print=F)))

d23_ratio_gmc_il10_N_tr

d23_ratio_il12_il4_N_tr<-lab %>%
  subset(d23_ratio_il12_il4!="NA") %>%
  group_by (tr) %>%
  do(as.data.frame(washb_mean(Y=.$d23_ratio_il12_il4,  id=.$block,  print=F)))

d23_ratio_il12_il4_N_tr

d23_ratio_ifn_il4_N_tr<-lab %>%
  subset(d23_ratio_ifn_il4!="NA") %>%
  group_by (tr) %>%
  do(as.data.frame(washb_mean(Y=.$d23_ratio_ifn_il4,  id=.$block,  print=F)))

d23_ratio_ifn_il4_N_tr

d23_ratio_il12_il5_N_tr<-lab %>%
  subset(d23_ratio_il12_il5!="NA") %>%
  group_by (tr) %>%
  do(as.data.frame(washb_mean(Y=.$d23_ratio_il12_il5, id=.$block,  print=F)))

d23_ratio_il12_il5_N_tr

d23_ratio_ifn_il5_N_tr<-lab %>%
  subset(d23_ratio_ifn_il5!="NA") %>%
  group_by (tr) %>%
  do(as.data.frame(washb_mean(Y=.$d23_ratio_ifn_il5,  id=.$block,  print=F)))

d23_ratio_ifn_il5_N_tr

d23_ratio_il12_il13_N_tr<-lab %>%
  subset(d23_ratio_il12_il13!="NA") %>%
  group_by (tr) %>%
  do(as.data.frame(washb_mean(Y=.$d23_ratio_il12_il13,  id=.$block,  print=F)))

d23_ratio_il12_il13_N_tr

d23_ratio_ifn_il13_N_tr<-lab %>%
  subset(d23_ratio_ifn_il13!="NA") %>%
  group_by (tr) %>%
  do(as.data.frame(washb_mean(Y=.$d23_ratio_ifn_il13, id=.$block,  print=F)))

d23_ratio_ifn_il13_N_tr

d23_ratio_il12_il17_N_tr<-lab %>%
  subset(d23_ratio_il12_il17!="NA") %>%
  group_by (tr) %>%
  do(as.data.frame(washb_mean(Y=.$d23_ratio_il12_il17, id=.$block,  print=F)))

d23_ratio_il12_il17_N_tr

d23_ratio_ifn_il17_N_tr<-lab %>%
  subset(d23_ratio_ifn_il17!="NA") %>%
  group_by (tr) %>%
  do(as.data.frame(washb_mean(Y=.$d23_ratio_ifn_il17, id=.$block,  print=F)))

d23_ratio_ifn_il17_N_tr

d23_ratio_il12_il21_N_tr<-lab %>%
  subset(d23_ratio_il12_il21!="NA") %>%
  group_by (tr) %>%
  do(as.data.frame(washb_mean(Y=.$d23_ratio_il12_il21, id=.$block,  print=F)))

d23_ratio_il12_il21_N_tr

d23_ratio_ifn_il21_N_tr<-lab %>%
  subset(d23_ratio_ifn_il21!="NA") %>%
  group_by (tr) %>%
  do(as.data.frame(washb_mean(Y=.$d23_ratio_ifn_il21,  id=.$block,  print=F)))

d23_ratio_ifn_il21_N_tr

d23_ratio_pro_il10_N_tr<-lab %>%
  subset(d23_ratio_pro_il10!="NA") %>%
  group_by (tr) %>%
  do(as.data.frame(washb_mean(Y=.$d23_ratio_pro_il10, id=.$block,  print=F)))

d23_ratio_pro_il10_N_tr

d23_ratio_th1_il10_N_tr<-lab %>%
  subset(d23_ratio_th1_il10!="NA") %>%
  group_by (tr) %>%
  do(as.data.frame(washb_mean(Y=.$d23_ratio_th1_il10,  id=.$block,  print=F)))

d23_ratio_th1_il10_N_tr

d23_ratio_th2_il10_N_tr<-lab %>%
  subset(d23_ratio_th2_il10!="NA") %>%
  group_by (tr) %>%
  do(as.data.frame(washb_mean(Y=.$d23_ratio_th2_il10, id=.$block,  print=F)))

d23_ratio_th2_il10_N_tr

d23_ratio_th17_il10_N_tr<-lab %>%
  subset(d23_ratio_th17_il10!="NA") %>%
  group_by (tr) %>%
  do(as.data.frame(washb_mean(Y=.$d23_ratio_th17_il10,  id=.$block,  print=F)))

d23_ratio_th17_il10_N_tr

d23_ratio_th1_th2_N_tr<-lab %>%
  subset(d23_ratio_th1_th2!="NA") %>%
  group_by (tr) %>%
  do(as.data.frame(washb_mean(Y=.$d23_ratio_th1_th2,  id=.$block,  print=F)))

d23_ratio_th1_th2_N_tr

d23_ratio_th1_th17_N_tr<-lab %>%
  subset(d23_ratio_th1_th17!="NA") %>%
  group_by (tr) %>%
  do(as.data.frame(washb_mean(Y=.$d23_ratio_th1_th17,  id=.$block,  print=F)))

d23_ratio_th1_th17_N_tr


save(list=as.vector(ls(pattern="N_tr")), file=here('results/immune_tr_means.RData'))

