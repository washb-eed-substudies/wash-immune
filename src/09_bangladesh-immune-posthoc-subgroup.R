#---------------------------------------
# bangladesh-immune-subgroup.R
#
# andrew mertens (amertens@berkeley.edu)
#
#  Firdaus suggested that we conduct an EMM analysis by pathogen status, diarrheal infection, and fever to strengthen the paper and address some of the reviewers' points. 
#---------------------------------------

#---------------------------------------
# input files:
#	bangladesh-dm-ee-anthro-diar-ee-med-plasma-blind-tr-enrol-covariates-lab.csv
#
# output files:
#	bangladesh-immune-subgroup.Rdata
# 
# 
#---------------------------------------

#---------------------------------------
# preamble
#---------------------------------------

rm(list=ls())
source(here::here("0-config.R"))


setwd(paste0(dropboxDir,"Data/Cleaned/Audrie/")) #Set working directory


dfull <- readRDS("C:/Users/andre/OneDrive/Documents/washb_substudies/eed-substudy-data/bangladesh-cleaned-master-data.RDS")

d_illness <- dfull %>% select(childid, diar7d_t2, diar7d_t3, fever7d_t2, fever7d_t3) %>% mutate(illness_data=1)



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
d_path <- left_join(EE_df, ch, by=c("childid"))
dim(d_path)
d_path <- d_path %>% filter(!is.na(pathogen_data)) %>% droplevels()
dim(d_path)

#merge fever and diarrhea data
dim(d_illness)
d_illness <- left_join(EE_df, d_illness, by=c("childid"))
dim(d_illness)
d_illness <- d_illness %>% filter(!is.na(illness_data)) %>% droplevels()
dim(d_illness)




XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
Keep working from here
XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX








#---------------------------------------
# subset to the relevant measurement
# T2 
#---------------------------------------

dim(d)



# re-order the tr factor for convenience
d$tr <- factor(d$tr,levels=c("Control","Nutrition + WSH"))

# sort the data for perfect replication with andrew on the V-fold cross-validation
d <- d[order(d$block,d$clusterid,d$dataid,d$childid),]


###################
#Subgroup Analyses @T2
###################


#Select adjustment covariates 
Wvars<-c("sex")

#subset the main dataframe and create a new W dataframe
W<- subset(d, select=Wvars)

#check that all the factor variables are set
for(i in 1:ncol(W)){
  cat(colnames(W)[i],"  ",class(W[,i]),"\n")
}

#Looks good. Use this if any need to be changes:


W$sex<-as.factor(W$sex)

# Set up the WASHB function
# df=data frame

# stratified by "sex"

washb_function <- function(df,x) {
  
  temp <- washb_glm(Y=d[,x], tr=d$tr, pair=NULL, W=d["sex"], V="sex", id=d$block, contrast = c("Control","Nutrition + WSH"), family="gaussian", verbose=FALSE)
  
  temp_metric<-as.data.frame(temp$lincom)
  
  colnames(temp_metric) <-c("subgroup", "RD","SE","ci.lb","ci.ub","z","P-value")
  temp_metric$int_Pval <- c( temp$fit[4,6],NA)
  return(temp_metric)
}

#grab the variables with prefix 't2_' from the data frame and then apply the washb_function
list_immune <- lapply(names(d)[grep('t2_', names(d))],  function(x) washb_function(d,x))

list_immune

#put names of each of the variables into the matrix
names(list_immune) <- names(d)[grep('t2_', names(d))]

#resulting matrix
list_immune
