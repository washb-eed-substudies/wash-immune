
rm(list=ls())
source(here::here("0-config.R"))
library(boxr)
library(washbgam)
library(tableone)
library(table1)
library(kableExtra)
library(knitr)
library(labelled)
library(lubridate)
library(mediation)

library("survival")
library("multcomp")
library("sandwich")
library("MASS")

box_auth()

dfull <- box_read("871638120165") 

EE_df <- dfull %>% subset(., select=c(tr, childid:d23_ratio_th1_th17))
colnames(EE_df)

#To do:

#Figure out which pathogens to use as the outcome
#Exactly which immune markers as mediatiors
#Check merge and time-ordering
   # I assume we should be using 


#Year 2 pathogens, STH, giardia

#if we can show that the interventions induced a TH2 response/immunoregulatory response at
#Y1 which then reduced pathogens at Y2, that would strengthen the paper
#Also IL-10


#Use significant year 1 from here:
#https://www.medrxiv.org/content/10.1101/2021.11.10.21266206v1.full

#Use the one-year significant in the outcome:
#After one year, TNF-??/IL-10, IL-12/IL-10, and IL-17A/IL-10 ratios were lower in the intervention group compared to the control group 
#(mean difference: -0.12 to -0.19, p<0.05), indicating the intervention promoted IL-10 driven immunoregulation. 
#Similar reductions in ratios of pro-inflammatory cytokines to IL-10 were sustained in the intervention group after two years. 
#After one year, IL-12/IL-4, IL-12/IL-5, IFN-??/IL-5, and IL-12/IL-13 ratios were lower in the intervention group 
#(???0.18 to -0.27, p<0.05), suggesting a shift towards a Th2 cytokine response.


####
# Add in:
# th1/th2
# Pro/IL-10
# Th1/IL-10
# Th2/IL-10



# #--------------------------------------------------------------
# # Clean TAC data to pathogens in the environment
# #--------------------------------------------------------------
# EE_TAC <- readRDS("C:/Users/andre/Dropbox/IPD WASH/Data/WBB/data.UVA.relativeCt_spikeinAdjusted_Caitlin.RDS")
# 
# ch <- EE_TAC %>% 
#   mutate(ch_pos_path_ecoli=1*(ETEC_LT+EAEC+ETEC_ST+ETEC.any +tEPEC+ aEPEC+ EPEC.any+STEC > 0),
#          ch_pos_adenovirus=1*(`Adenovirus 40/41` +  `Adenovirus pan`> 0)) %>%
#   rename(ch_abund_giardia_EE=Giardia,
#          ch_abund_norovirus=Norovirus.any,
#          ch_abund_crypto_EE =Cryptosporidium,
#          ch_abund_salmonella=Salmonella,
#          ch_abund_ascaris_EE=Ascaris,
#          ch_abund_trichuris_EE=Trichuris,
#          ch_abund_entamoeba_EE=E.histolytica,
#          ch_abund_cholera=V.cholerae,
#          ch_abund_cdiff=C.difficile,
#          ch_abund_shigella=Shig_EIEC,
#          ch_abund_rotavirus=Rotavirus,
#          ch_abund_campylobacter=Campy.pan) %>%
#   mutate(
#     ch_pos_giardia_EE=1*(ch_abund_giardia_EE> 0),
#     ch_pos_norovirus=1*(ch_abund_norovirus> 0),
#     ch_pos_crypto_EE =1*(ch_abund_crypto_EE> 0),
#     ch_pos_salmonella=1*(ch_abund_salmonella> 0),
#     ch_pos_ascaris_EE=1*(ch_abund_ascaris_EE> 0),
#     ch_pos_trichuris_EE=1*(ch_abund_trichuris_EE> 0),
#     ch_pos_entamoeba_EE=1*(ch_abund_entamoeba_EE> 0),
#     ch_pos_cholera=1*(ch_abund_cholera> 0),
#     ch_pos_cdiff=1*(ch_abund_cdiff> 0),
#     ch_pos_shigella=1*(ch_abund_shigella> 0),
#     ch_pos_rotavirus= 1*(ch_abund_rotavirus> 0),
#     ch_pos_campylobacter=1*(ch_abund_campylobacter> 0)
#   ) %>%
#   subset(., select = -c(ETEC_LT,EAEC,ETEC_ST,ETEC.any,tEPEC, aEPEC, EPEC.any,STEC,Ancyclostoma,Necator,E.bieneusi,
#                         Blastocystis,Plesiomonas,Aeromonas,PhHV,MS2,Schistosoma,H.nana,`Adenovirus 40/41`, `Adenovirus pan`,
#                         Cyclospora, Isospora, H.pylori, Sapovirus, `Cryptosporidium hominis`, `Cryptosporidium parvum`,
#                         E.intestinalis,`pan Entamoeba`, Astrovirus, Strongyloides, M.tuberculosis, `Norovirus GI`, `Norovirus GII`,
#                         B.fragilis, `Campylobacter jejuni/coli`))
# colnames(ch)



#load clean dataset from the diarhhea-pathogen analysis
wbb_sth <- readRDS(paste0("C:/Users/andre/Dropbox/IPD WASH/Data/WBB/clean/clean_wbb_sth_pathogens.rds"))
# and the public conversion IDs
IDs <- read.csv(paste0("C:/Users/andre/Dropbox/IPD WASH/Data/WBB/public-ids.csv"))
wbb_sth <- wbb_sth %>% rename(block_r=block,	clusterid_r=clusterid,	dataid_r=dataid)
wbb_sth <- left_join(wbb_sth, IDs, by = c("block_r",	"clusterid_r",	"dataid_r"))
head(wbb_sth)

colnames(wbb_sth)

table(wbb_sth$hhid)

table(wbb_sth$tr)
table(wbb_sth$tr, wbb_sth$qpcr.positive.Ac)


ch <- wbb_sth %>% subset(., select = c(
  dataid,hhid,childid, svyweek, svyyear, sex,             agedays,         
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


# #get midline stool dates from ee dataset and merge in
# EEdates <- read.csv("C:/Users/andre/Dropbox/WASHB-EE-analysis/WBB-EE-analysis/Data/Cleaned/Andrew/BD-EE-stool.csv")
# EEdates <- EEdates %>% subset(., select = c(childid, date2)) %>% 
#   rename(child_date_pathogen=date2) %>% mutate(child_date_pathogen=dmy(child_date_pathogen)) %>%
#   filter(!is.na(child_date_pathogen))
# 
# #merge:
# dim(EE_TAC)
# dim(EEdates)
# EE_TAC<-left_join(EE_TAC, EEdates, by = c("childid")) %>% filter(!is.na(child_date_pathogen))
# dim(EE_TAC)
# head(EE_TAC)


colnames(ch)


summary(EE_df$childid)
table(ch$childid)
summary(ch$dataid)

table(as.numeric())

ch$childid <- as.numeric(gsub("T","",ch$childid))
ch$childid <- ch$dataid*10 + ch$childid

ch <- ch %>% mutate(pathogen_data=1)

dim(EE_df)
dim(ch)
d <- left_join(EE_df, ch, by=c("childid"))
dim(d)
d <- d %>% filter(!is.na(pathogen_data)) %>% droplevels()
dim(d)

df <- anti_join(EE_df, ch, by=c("childid"))
dim(df)
df <- df %>% arrange(childid)

df2 <- anti_join( ch, EE_df, by=c("childid"))

unique(df$childid)
unique(df2$childid)


#TEMP
d <- d %>% filter(!is.na(t2_ratio_tnf_il10))



table(is.na(d$ch_pos_ascaris ))

table((d$ch_pos_ascaris ))
table((d$ch_qpcr_pos_ascaris ))
table(d$tr)

table(d$ch_qpcr_pos_ascaris, d$tr)
table(d$ch_pos_ascaris, d$tr)



#Yes - if we can show that the interventions induced a TH2 response/immunoregulatory response at Y1 which then reduced pathogens at Y2, that would strengthen the paper.

#Another reviewer point: "No attempt is made to dissect which part of the intervention might have contributed." If we only have nutritional markers (e.g., vit A, B, iron, anemia) at 
#Y2, does it make sense to look at the association between nutritional markers at Y2 and immune response at Y2 to try to get at this answer? We don't have the temporal ordering so not sure if it's worth it. 
#Thoughts? I suppose we could also go the route of verbally explaining that we didn't have the funding to look at the nutrition only and WSH only arms but that we do have those samples collected




######### Mediate package
# https://towardsdatascience.com/doing-and-reporting-your-first-mediation-analysis-in-r-2fe423b92171


dfull = d
mediator = 't2_ratio_tnf_il10'
outcome = "ch_pos_giardia"

wash_mediate <- function(dfull, mediator, outcome){

  # df <- dfull %>% subset(., select=c(!!(outcome),"tr",!!(ch_pos_giardia)))
  
  Y <- subset(dfull, select = outcome)
  colnames(Y) <- "Y"
  tr <- subset(dfull, select = "tr")
  Z <- subset(d, select = mediator)
  colnames(Z) <- "Z"
  
  df <- bind_cols(Y,tr,Z)
  
  
  df <- df[complete.cases(df),]
  df$tr <- ifelse(df$tr=="Control",0,1)
  
  fit.totaleffect=glm(Y~tr, family="binomial", data=df)
  summary(fit.totaleffect)
  fit.mediator=glm(Z~tr, family="gaussian", data=df)
  summary(fit.mediator)
  fit.dv=glm(Y~tr+Z, family="binomial", data=df)
  summary(fit.dv)
  
  results = mediate(fit.mediator, fit.dv, treat='tr', mediator="Z", boot=T)
  res <- summary(results)
  summary(res)
  
  
  res.tab <- bind_rows(
    data.frame(type="TE", est=res$tau.p, ci.lb=res$tau.ci[1],  ci.ub=res$tau.ci[2], pval=res$tau.p),
    data.frame(type="ACME", est=res$d1,  ci.lb=res$d1.ci[1],  ci.ub=res$d1.ci[2],  pval=res$d1.p),
    data.frame(type="ADE", est=res$z0,  ci.lb=res$z0.ci[1],  ci.ub=res$z0.ci[2],  pval=res$z0.p),
    data.frame(type="Prop", est=res$n0,  ci.lb=res$n0.ci[1],  ci.ub=res$n0.ci[2],  pval=res$n0.p)
  )
  rownames(res.tab) <- NULL
  res.tab$Y <- outcome
  res.tab$Z <- mediator
  
  return(res.tab)
}





# After one year, TNF/IL-10, IL-12/IL-10, and IL-17A/IL-10 ratios were lower in the intervention group compared to the control group (mean difference: -0.12 to -0.19, p<0.05), 
# indicating the intervention promoted IL-10 driven immunoregulation. Similar reductions in ratios of pro-inflammatory cytokines to IL-10 were sustained in the intervention group after two years. 
# After one year, IL-12/IL-4, IL-12/IL-5, IFN/IL-5, and IL-12/IL-13 ratios were lower in the intervention group (???0.18 to -0.27, p<0.05), suggesting a shift towards a Th2 cytokine response.
immune_markers <- c("t2_ratio_tnf_il10", "t2_ratio_il17_il10", "t2_ratio_il12_il10", "t2_ratio_il12_il4", "t2_ratio_il12_il5", "t2_ratio_ifn_il5","t2_ratio_il12_il13")

#any pathogen
d$any_pathogen <- ifelse(d$ch_pos_ascaris==1| d$ch_pos_trichuris==1|d$ch_pos_giardia==1, 1, 0)
table(d$any_pathogen)

res_full <- NULL
for(i in immune_markers){
  for(j in c("any_pathogen","ch_pos_ascaris","ch_pos_trichuris","ch_pos_giardia")){
    res_temp <- NULL
    res_temp <- wash_mediate(dfull = d, mediator = i, outcome = j)
    res_full <- bind_rows(res_full, res_temp)
  }
}

res_full

res_full %>% filter(type=="ACME", Y=="any_pathogen")
res_full %>% filter(type=="ADE", Y=="any_pathogen")
res_full %>% filter(type=="TE")

res_full %>% filter(pval  < 0.05)


saveRDS(res_full, file=here::here("results/mediation.RDS"))
write.csv(res_full, file=here::here("results/mediation_results.csv"))



# d1 point estimate for average causal mediation effects.
# d1.ci confidence intervals for average causal mediation effect. The confidence level is
# set at the value specified in 'conf.level'.
# z0 point estimates for average direct effect.
# z0.ci confidence intervals for average direct effect.
# z0.p two-sided p-values for average causal direct effect.
# n0 the "proportions mediated", or the size of the average causal mediation effect
# relative to the total effect.
# n0.ci confidence intervals for the proportion mediated.
# n0.p two-sided p-values for proportion mediated.
# tau.coef point estimate for total effect.
# tau.ci confidence interval for total effect.
# tau.p two-sided p-values for total effect.
# boot logical, the boot argument used.
# treat a character string indicating the name of the 'treat' variable used.
# mediator a character string indicating the name of the 'mediator' variable used.
# INT a logical value indicating whether the model specification allows the effects to
# differ between the treatment and control conditions.
# conf.level the confidence level used.
# model.y the outcome model used.
# model.m the mediator model used.
# nobs number of observations in the model frame for 'model.m' and 'model.y'. May
# differ from the numbers in the original models input to 'mediate' if 'dropobs'
# was 'TRUE'.
# cluster the clusters used.
