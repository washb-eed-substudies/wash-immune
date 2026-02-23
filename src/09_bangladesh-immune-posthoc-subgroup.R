#---------------------------------------
# bangladesh-immune-subgroup.R
#
# andrew mertens (amertens@berkeley.edu)
#
#  Firdaus suggested that we conduct an EMM analysis by pathogen status, diarrheal infection, 
#and fever to strengthen the paper and address some of the reviewers' points. 
#---------------------------------------




# d <- d %>% mutate(
#   group=case_when(
#     outcome %in% c("ratio_th1_th2", "ratio_il12_il4", "ratio_ifn_il4", "ratio_il12_il5", "ratio_ifn_il5",  "ratio_il12_il13", "ratio_ifn_il13") ~"one", 
#     outcome %in% c("ratio_pro_il10", "ratio_il1_il10","ratio_il6_il10", "ratio_tnf_il10", 
#                    "ratio_il2_il10",  "ratio_th1_il10",  "ratio_th2_il10", "ratio_il12_il10", "ratio_ifn_il10", "ratio_il4_il10", 
#                    "ratio_il5_il10", "ratio_il13_il10","ratio_th17_il10", "ratio_il17_il10", "ratio_il21_il10",
#                    "ratio_gmc_il10") ~ "two",
#     outcome %in% c("ratio_th1_th17", "ratio_il12_il17", "ratio_ifn_il17", "ratio_il12_il21", "ratio_ifn_il21") ~"three"
#   ),
#   group=factor(group, level=c("one","two","three")),
#   outcome = factor(outcome, levels =c(
#     "ratio_th1_th2", "ratio_il12_il4", "ratio_ifn_il4", "ratio_il12_il5", "ratio_ifn_il5",  "ratio_il12_il13", "ratio_ifn_il13",
#     "ratio_pro_il10", "ratio_il1_il10","ratio_il6_il10", "ratio_tnf_il10", 
#     "ratio_il2_il10",  "ratio_th1_il10",  "ratio_th2_il10", "ratio_il12_il10", "ratio_ifn_il10", "ratio_il4_il10", 
#     "ratio_il5_il10", "ratio_il13_il10","ratio_th17_il10", "ratio_il17_il10", "ratio_il21_il10",
#     "ratio_gmc_il10",
#     "ratio_th1_th17", "ratio_il12_il17", "ratio_ifn_il17", "ratio_il12_il21", "ratio_ifn_il21"))
# )



#---------------------------------------
# preamble
#---------------------------------------

rm(list=ls())
source(here::here("0-config.R"))


setwd(paste0(dropboxDir,"Data/Cleaned/Audrie/")) #Set working directory


#dfull <- readRDS("C:/Users/andre/OneDrive/Documents/washb_substudies/eed-substudy-data/bangladesh-cleaned-master-data.RDS")
dfull <- read.csv("C:/Users/andre/Downloads/bangladesh-merged-eed-mn-data.csv")

Wvars<-c("monsoon_bt2","ageday_bt2", "sex","birthord", "momage","momheight","momedu", "hfiacat", "Nlt18","Ncomp", "watmin", "walls", "floor", "elec", "asset_wardrobe", "asset_table", "asset_chair", "asset_clock","asset_khat", "asset_chouki", "asset_radio", "asset_tv", "asset_refrig", "asset_bike", "asset_moto", "asset_sewmach", "asset_mobile", "n_cattle", "n_goat", "n_chicken")


d_illness <- dfull %>% select(childid,diar7d_t2, diar7d_t3, fever7d_t2, fever7d_t3) %>% mutate(illness_data=1)



#load clean dataset from the diarhhea-pathogen analysis
wbb_sth <- readRDS(paste0("C:/Users/andre/OneDrive/Documents/washb_substudies/eed-substudy-data/sth data/clean_wbb_sth_pathogens.rds"))
# and the public conversion IDs
IDs <- read.csv(paste0("C:/Users/andre/OneDrive/Documents/washb_substudies/eed-substudy-data/public-ids.csv"))
wbb_sth <- wbb_sth %>% rename(block_r=block,	clusterid_r=clusterid,	dataid_r=dataid)
wbb_sth <- left_join(wbb_sth, IDs, by = c("block_r",	"clusterid_r",	"dataid_r"))
head(wbb_sth)

colnames(wbb_sth)

table(wbb_sth$hhid)

table(wbb_sth$tr)
table(wbb_sth$tr, wbb_sth$qpcr.positive.Ac)


ch <- wbb_sth %>% subset(., select = c(
  dataid,hhid,childid, svyweek, svyyear,           agedays,         
  clusterid,       logalepg,        loghwepg,        
  logttepg,       posgi,           poseh,           poscr,           
  posprot,         posmult,         ctgi,            cteh,            
  ctcr,            qpcr.positive.Ac,qpcr.positive.Ad,qpcr.positive.Al,
  qpcr.positive.IAC,qpcr.positive.Na,qpcr.positive.Ss,qpcr.positive.Tt,
  qpcr.positive.Hw,qpcr.positive.Sth,qpcr.CTmean.Ac,  qpcr.CTmean.Ad,  
  qpcr.CTmean.Al,  qpcr.CTmean.IAC, qpcr.CTmean.Na,  qpcr.CTmean.Ss,  
  qpcr.CTmean.Tt)) %>% 
  filter(!is.na(logalepg) | !is.na(loghwepg) | !is.na(logttepg) | 
           !is.na(posgi) | !is.na(poseh) | !is.na(poscr) | !is.na(posprot) | 
           !is.na(posmult) | !is.na(qpcr.positive.Ac) | !is.na(qpcr.positive.Ad) |
           !is.na(qpcr.positive.Al) | !is.na(qpcr.positive.IAC) | !is.na(qpcr.positive.Na) |   
           !is.na(qpcr.positive.Ss) | !is.na(qpcr.positive.Tt) | !is.na(qpcr.positive.Sth))  %>%
  rename(#ch_ascaris=ascaris_yn, 
    # ch_trichuris=trichuris_yn, 
    # ch_hook=hook_yn, 
    # ch_sth=sth_yn, 
    # ch_giardia=posgi, 
    # ch_sth_giar_coinf=posmult,
    aged_pathogen = agedays) %>%
  mutate(month=ceiling(svyweek/52*12),
         pathogen_date = dmy(paste0("15-",month, "-", svyyear))) %>%
  subset(., select = -c(month))  %>% 
  subset(., select = -c(posmult,
                        loghwepg,
                        qpcr.CTmean.Na, #don't have hookworm in soil
                        qpcr.CTmean.Ad,
                        qpcr.positive.Ad,
                        qpcr.positive.Na,
                        qpcr.positive.Hw,
                        qpcr.positive.Sth,
                        qpcr.CTmean.Ss, #no threadworm
                        qpcr.positive.Ss,
                        qpcr.positive.Ac, #Figure out what Ac and IAC are  
                        qpcr.positive.IAC,  
                        qpcr.CTmean.Ac,      
                        qpcr.CTmean.IAC)) %>%
  rename(
    ch_abund_ascaris=logalepg,
    ch_abund_trichuris=logttepg,
    ch_pos_giardia=posgi,
    ch_pos_entamoeba=poseh,
    ch_pos_crypto=poscr,
    ch_abund_giardia=ctgi,
    ch_abund_entamoeba=cteh,
    ch_abund_crypto=ctcr,
    ch_qpcr_pos_trichuris=qpcr.positive.Tt,
    ch_qpcr_abund_trichuris=qpcr.CTmean.Tt,
    ch_qpcr_pos_ascaris=qpcr.positive.Al,
    ch_qpcr_abund_ascaris=qpcr.CTmean.Al) %>%
  mutate(ch_pos_ascaris = 1*(ch_abund_ascaris>0), 
         ch_pos_trichuris = 1*(ch_abund_trichuris>0))


ch$childid <- as.numeric(gsub("T","",ch$childid))
ch$childid <- ch$dataid*10 + ch$childid

ch <- ch %>% mutate(pathogen_data=1)




#---------------------------------------
# Load the analysis dataset,
# the baseline covariate dataset
#---------------------------------------

EE_df<-readRDS(paste0(dropboxDir,"Data/Cleaned/Audrie/bangladesh-immune-analysis-dataset.rds"))

#merge pathogen data
dim(EE_df)
dim(ch)
ch$flag<-1
d_path <- left_join(EE_df, ch, by=c("childid"))
table(is.na(d_path$flag))
dim(d_path)
d_path <- d_path %>% filter(!is.na(pathogen_data)) %>% droplevels()
dim(d_path)

#merge fever and diarrhea data
dim(d_path)
dim(d_illness)
d_illness <- left_join(d_illness,d_path, by=c("childid"))
dim(d_illness)

table(is.na(d_illness$pathogen_data))
table(is.na(d_illness$illness_data))
table(is.na(d_illness$sumscore_t2_Z))
table(is.na(d_illness$pathogen_data),is.na(d_illness$sumscore_t2_Z))


d_illness <- d_illness %>% filter(!is.na(illness_data)) %>% droplevels()
dim(d_illness)


Vvars <- c("diar7d_t2", "diar7d_t3", "fever7d_t2", "fever7d_t3",
           "ch_pos_giardia",
           "ch_pos_entamoeba",
           "ch_pos_crypto",
           "ch_qpcr_pos_trichuris",
           "ch_qpcr_pos_ascaris")




d <- d_illness
dim(d)

table(d$diar7d_t2)
d$diar7d_t2 <- factor(d$diar7d_t2, levels=c("0","1","Missing"))
d$diar7d_t3 <- factor(d$diar7d_t3, levels=c("0","1","Missing"))
table(d$fever7d_t2)

d$fever7d_t2 <- as.character(d$fever7d_t2)
d$fever7d_t2[is.na(d$fever7d_t2)] <- "Missing"
d$fever7d_t2 <- factor(d$fever7d_t2, levels=c("0","1","Missing"))
table(d$fever7d_t2)

d$fever7d_t3 <- as.character(d$fever7d_t3)
d$fever7d_t3[is.na(d$fever7d_t3)] <- "Missing"
d$fever7d_t3 <- factor(d$fever7d_t3, levels=c("0","1","Missing"))
table(d$fever7d_t3)

d$ch_pos_giardia <- as.character(d$ch_pos_giardia)
d$ch_pos_giardia[is.na(d$ch_pos_giardia)] <- "Missing"
d$ch_pos_giardia <- factor(d$ch_pos_giardia, levels=c("0","1","Missing"))
table(d$ch_pos_giardia)

d$ch_pos_entamoeba <- as.character(d$ch_pos_entamoeba)
d$ch_pos_entamoeba[is.na(d$ch_pos_entamoeba)] <- "Missing"
d$ch_pos_entamoeba <- factor(d$ch_pos_entamoeba, levels=c("0","1","Missing"))
table(d$ch_pos_entamoeba)

d$ch_pos_crypto <- as.character(d$ch_pos_crypto)
d$ch_pos_crypto[is.na(d$ch_pos_crypto)] <- "Missing"
d$ch_pos_crypto <- factor(d$ch_pos_crypto, levels=c("0","1","Missing"))
table(d$ch_pos_crypto)

d$ch_qpcr_pos_trichuris <- as.character(d$ch_qpcr_pos_trichuris)
d$ch_qpcr_pos_trichuris[is.na(d$ch_qpcr_pos_trichuris)] <- "Missing"
d$ch_qpcr_pos_trichuris <- factor(d$ch_qpcr_pos_trichuris, levels=c("0","1","Missing"))
table(d$ch_qpcr_pos_trichuris)

d$ch_qpcr_pos_ascaris <- as.character(d$ch_qpcr_pos_ascaris)
d$ch_qpcr_pos_ascaris[is.na(d$ch_qpcr_pos_ascaris)] <- "Missing"
d$ch_qpcr_pos_ascaris <- factor(d$ch_qpcr_pos_ascaris, levels=c("0","1","Missing"))
table(d$ch_qpcr_pos_ascaris)




#outcomes at time 3- 


d$monsoon_bt3<-as.factor(d$monsoon_bt3)
d$monsoon_bt3<-addNA(d$monsoon_bt3)
levels(d$monsoon_bt3)[length(levels(d$monsoon_bt3))]<-"Missing"



colnames(d)
W<- subset(d, select=c(Wvars,"monsoon_bt3","ageday_bt3", "fever7d_t2"))

temp <- washb_glm(Y=d$t3_ln_igf, tr=d$tr, pair=NULL, W=W, V="fever7d_t2", id=d$block, contrast = c("Control","Nutrition + WSH"), family="gaussian", verbose=FALSE)
temp_metric<-as.data.frame(temp$lincom[1:2,])
colnames(temp_metric) <-c("subgroup", "RD","SE","ci.lb","ci.ub","z","P-value")
temp_metric$int_Pval <- c( temp$fit[4,6],NA)

washb_function <- function(df,Y, v) {
  
  W<- subset(df, select=c(Wvars, v))
  
  mean_res=df %>% filter(!is.na(!!sym(Y)),  !!sym(v)!="Missing") %>% group_by(tr, !!sym(v)) %>% summarise(n=n(), mean=mean(!!sym(Y)), sd=sd(!!sym(Y)))
  colnames(mean_res)[2] = "subgroup"
  
  mean_res <- mean_res %>%
    pivot_wider(
      id_cols = subgroup,
      names_from = tr,
      values_from = c(n, mean, sd),
      names_glue = "{.value}_{tr}"
    )
  
  colnames(mean_res) <- gsub("Nutrition \\+ WSH","NWSH", colnames(mean_res))
  
  temp <- washb_glm(Y=d[,Y], tr=d$tr, pair=NULL, W=W, V=v, id=d$block, contrast = c("Control","Nutrition + WSH"), family="gaussian", print=F, verbose=FALSE)
  temp_metric<-as.data.frame(temp$lincom[1:2,])
  colnames(temp_metric) <-c("subgroup", "RD","SE","ci.lb","ci.ub","z","P-value")
  temp_metric$int_Pval <- c( temp$fit[4,6],NA)
  temp_metric$Y <- Y
  temp_metric$V <- v
  temp_metric= left_join(mean_res, temp_metric, by="subgroup")
  temp_metric
  return(temp_metric)
}

t3_Y = names(d)[grep('t3_', names(d))]
t3_Y=paste0("t3_", c(
  "ratio_th1_th2", "ratio_il12_il4", "ratio_ifn_il4", "ratio_il12_il5", "ratio_ifn_il5",  "ratio_il12_il13", "ratio_ifn_il13",
  "ratio_pro_il10", "ratio_il1_il10","ratio_il6_il10", "ratio_tnf_il10", 
  "ratio_il2_il10",  "ratio_th1_il10",  "ratio_th2_il10", "ratio_il12_il10", "ratio_ifn_il10", "ratio_il4_il10", 
  "ratio_il5_il10", "ratio_il13_il10","ratio_th17_il10", "ratio_il17_il10", "ratio_il21_il10",
  "ratio_gmc_il10",
  "ratio_th1_th17", "ratio_il12_il17", "ratio_ifn_il17", "ratio_il12_il21", "ratio_ifn_il21"))

Y=i=t3_Y[1]
v=h="ch_qpcr_pos_trichuris"
df=d

full_res=NULL
for(i in t3_Y){
  for(j in Vvars){
    cat(i,"\n")
    cat(j,"\n")
    temp_res=NULL
    temp_res=washb_function(d, i, j)
    full_res=bind_rows(full_res, temp_res)
  }
}

full_res
dim(full_res)

write.csv(full_res, here("results/bangladesh-immune-posthoc-subgroup-results-t3.csv"), row.names=F)

