
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


#TEMP
d <- d %>% filter(!is.na(t2_ratio_th2_il10))



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

d2 <- d %>% subset(., select=c(ch_pos_ascaris,tr,t2_ratio_th2_il10))
d2 <- d2[complete.cases(d2),]
d2$tr <- ifelse(d2$tr=="Control",0,1)

library(mediation)
fit.totaleffect=glm(ch_pos_ascaris~tr, family="binomial", data=d2)
summary(fit.totaleffect)
fit.mediator=glm(t2_ratio_th2_il10~tr, family="gaussian", data=d2)
summary(fit.mediator)
fit.dv=glm(ch_pos_ascaris~tr+t2_ratio_th2_il10, family="binomial", data=d2)
summary(fit.dv)

results = mediate(fit.mediator, fit.dv, treat='tr', mediator='t2_ratio_th2_il10', boot=T)
res <- summary(results)
summary(res)

##### TRADITIONAL ANALYSIS FOR TOTAL EFFECT OF TREATMENT
TE.asc <- glm(ch_pos_ascaris ~ tr , family="binomial", data=d)
summary(TE.asc)


##### BARON AND KENNEY METHOD
BK.asc <- glm(ch_pos_ascaris ~ tr + t2_ratio_th2_il10, family="binomial", data=d)
summary(BK.asc)


#t2_ratio_th2_il10e=mediator
d$tr <- ifelse(d$tr=="Control",0,1)

d$Z <- ifelse(d$t2_ratio_th2_il10> mean(d$t2_ratio_th2_il10),1,0)


############# DIRECT/INDIRECT EFFECTS BY METHOD OF LANGE ET AL. 2012 AJE
##### MSM for NDE, NIE
## Step 1: Fit model for exposure conditional on  covariates
fit.A.D <- glm(tr ~ aged_pathogen + sex , data=d, family=binomial)
## Marginal model for stabilized weights
fit.A.N <- glm(tr ~ 1, data=d, family=binomial)
# Calculate exposure weights
w.A.N <- predict(fit.A.N, type="response") # for numerator
w.A.N[d$tr==0] <- 1-w.A.N[d$tr==0]

w.A.D <- predict(fit.A.D, type="response") # for denominator
w.A.D[d$tr==0] <- 1-w.A.D[d$tr==0]

d$w.A <- w.A.N/w.A.D

## Fit model for mediator, conditional on exposure and covariates
d$tr.tmp <- d$tr
fit.M <- glm(Z~tr.tmp + aged_pathogen + sex , data=d, family=binomial )
#fit.M <- glm(Z~tr.tmp + aged_pathogen + sex , data=d, family="gaussian" )

## Step 2: Create new dataset with a*, with all possible values of A
d.1 <- d.2 <- d
d.1$tr.star <- d$tr
d.2$tr.star <- 1-d$tr
d.all <- rbind(d.1, d.2)

d.all <- d.all[order(d.all$childid),]

## Step 3: Predicted probabilities for numerator (using tr.star)
##         and denominator (using tr)
d.all$tr.tmp <- d.all$tr # for denominator use tr for predicted values
w.M.D <- predict(fit.M, type="response", newdata=d.all)
w.M.D[d.all$Z==0] <- 1-w.M.D[d.all$Z==0] 

d.all$tr.tmp <- d.all$tr.star # for numerator use tr.star for predicted values
w.M.N <- predict(fit.M, type="response", newdata=d.all)
w.M.N[d.all$Z==0] <- 1-w.M.N[d.all$Z==0] 

d.all$w.M <- w.M.N/w.M.D

d.all$w <- d.all$w.A*d.all$w.M

summary(d.all[,c("w","w.A","w.M")])

## Step 4: MSM for natural effects
# NE.pathogen <-  coxph(Surv(timestrk,pathogen) ~ tr*tr.star + 
#                       cluster(childid), 
#                     weights=d.all$w, ties="efron", data=d.all)
# summary(NE.pathogen)

NE.asc <- glm(ch_pos_ascaris ~tr*tr.star + cluster(childid), family="binomial", weights=d.all$w, data=d.all)
summary(NE.asc)

table(d.all$tr, d.all$tr.star)

# Total direct/indirect effects
TDE.path <- matrix(c(0, 1, 0, 0, 1), 1)   # TDE of TREATMENT
round(exp(confint(glht(NE.asc, linfct=TDE.path))$confint)[,1:3],2)

TIE.path <- matrix(c(0,0, 1, 0, 1), 1)   # TIE of TREATMENT
round(exp(confint(glht(NE.asc, linfct=TIE.path))$confint)[,1:3],2)

# Pure direct/indirect effects
PDE.path <- matrix(c(0,1, 0, 0, 0), 1)   # PDE of TREATMENT
round(exp(confint(glht(NE.asc, linfct=PDE.path))$confint)[,1:3],2)

PIE.path <- matrix(c(0, 0, 1, 0, 0), 1)   # PIE of TREATMENT
round(exp(confint(glht(NE.asc, linfct=PIE.path))$confint)[,1:3],2)








