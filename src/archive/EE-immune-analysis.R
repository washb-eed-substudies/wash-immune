

rm(list=ls())
source(here::here("0-config.R"))



#load immune outcomes
imm<-read_dta(paste0(dropboxDir,"Data/Untouched/immune/washb-bangladesh-immun-lab-t2-t3.dta"))
imm <- as.data.frame(imm)
imm$childid <- as.numeric(imm$childid)
head(imm)

#log transform outcomes
table(is.na(imm[,-1]))
imm[,-1] <- log(imm[,-1])
table(is.na(imm[,-1]))

fulld <- read.csv(paste0(dropboxDir,"Data/Cleaned/Andrew/EE-BD_fulldata.csv"))
colnames(fulld)



#Merge in immune outcomes
dim(imm)
dim(fulld)
d <- left_join(imm, fulld, by="childid")
dim(d)

#Drop unneeded variables
d <- subset(d, select = -c(ur_agem1, ur_agem2, ur_agem3, st_agem1, st_agem2, st_agem3))

#Drop real treatment arm
# d <- subset(d, select = -c(tr))
# 
# #Merge in blinded treatment
# blind_tr <- read.csv(paste0(dropboxDir,"Data/Untouched/washb-BD-telo-blind-tr.csv"))
# blind_tr$clusterid <- as.numeric(blind_tr$clusterid)
# d <- left_join(d, blind_tr, by=c("block", "clusterid"))
# dim(d)
table(d$tr)

#Drop observation from one child erroneously sampled from Nutrition arm
d <- d %>% filter(tr != "Nutrition")
d <- droplevels(d)


#Load and merge in ages and dates at blood collection from 
ages <- read.csv(paste0(dropboxDir,"Data/Cleaned/Andrew/BD-EE-telo.csv"))
ages <- ages %>% subset(., select = c(dataid, childNo, aged2, agem2, aged3, agem3, month2, month3))
dim(d)
d <- left_join(d, ages, by=c("dataid", "childNo"))
dim(d)
table(is.na(d$agem2))
table(is.na(d$agem3))

#NOTE: TEMP
#Need to get ages for childid 30061 and 433011 from Audrie
da <- read.csv(paste0(dropboxDir,"Data/Cleaned/Audrie/bangladesh-dm-ee-anthro-diar-ee-med-plasma-blind-tr-enrol-covariates-lab.csv"))
da <- da %>% subset(., select = c(childid,agemth_bt2, agemth_bt3, month_t2, month_t3, ageday_bt2, ageday_bt3)) 

dim(d)
dim(da)
d <- full_join(d, da, by="childid")
dim(d)

#TEMP - need to make my child ages match Audrie's for one missing child
summary(d$agem2 - d$agemth_bt2)
summary(d$agem3 - d$agemth_bt3)
table(is.na(d$agem2))
table(is.na(d$agemth_bt2))
d$agem2 <- d$agemth_bt2
d$agem3 <- d$agemth_bt3

summary(d$aged2 - d$ageday_bt2)
summary(d$aged3 - d$ageday_bt3)
d$aged2 <- d$ageday_bt2
d$aged3 <- d$ageday_bt3

d$month2 <- d$month_t2
d$month3 <- d$month_t3
summary(d$month2 - d$month_t2)
summary(d$month3 - d$month_t3)

#Save data.frame
save(d, file = c(paste0(dropboxDir,"Data/Cleaned/Andrew/BD-EE-immune.Rdata")))

#load in names of Audrie's objects
nm <- list.files(path=here("Results/Audrie/"))
nm


#------------------------------------------------------------------------------------------------
# Age stats
#------------------------------------------------------------------------------------------------


age2_overall <- d %>% filter(!is.na(agem2)) %>%
  summarize(tr="Overall",N=n(), mean=mean(agem2), median=median(agem2), sd=sd(agem2), female=sum(sex==0), male=sum(sex))
age2_tr <- d %>% group_by(tr) %>% filter(tr %in% c("Control","Nutrition + WSH")) %>% filter(!is.na(agem2)) %>%
  summarize(N=n(), mean=mean(agem2), median=median(agem2), sd=sd(agem2), female=sum(sex==0), male=sum(sex))
age_t2_blood_M <- rbind(age2_overall, age2_tr)

age3_overall <- d %>% filter(!is.na(agem3)) %>%
  summarize(tr="Overall",N=n(), mean=mean(agem3), median=median(agem3), sd=sd(agem3), female=sum(sex==0), male=sum(sex))
age3_tr <- d %>% group_by(tr) %>% filter(tr %in% c("Control","Nutrition + WSH")) %>% filter(!is.na(agem3)) %>%
  summarize(N=n(), mean=mean(agem3), median=median(agem3), sd=sd(agem3), female=sum(sex==0), male=sum(sex))
age_t3_blood_M <- rbind(age3_overall, age3_tr)


#Compare to Audrie's
load(here("audrie results/immune-age-stats.RData"))
age_t2_blood_L[,-1] - age_t2_blood_M[,-1]
age_t3_blood_L[,-1] - age_t3_blood_M[,-1]

#Note: Audrie is not tabulating  N's and gender counts for only kids with ages/obs, but totals


#------------------------------------------------------------------------------------------------
# N's and means
#------------------------------------------------------------------------------------------------

#Set outcomes
outcomes <- c("igf_t2",          "igf_t3",          "crp_t2",         
  "agp_t2",          "gmcsf_t2",        "ifng_t2",         "il10_t2",         "il12_t2",        
  "il13_t2",         "il17_t2",         "il1_t2",          "il2_t2",          "il21_t2",        
  "il4_t2",          "il5_t2",          "il6_t2",          "tnfa_t2",         "gmcsf_t3",       
  "ifng_t3",         "il10_t3",         "il12_t3",         "il13_t3",         "il17_t3",        
  "il1_t3",          "il2_t3",          "il21_t3",         "il4_t3",          "il5_t3",         
  "il6_t3",          "tnfa_t3")


mean_sd <- d %>% select(outcomes) %>% summarise_all(funs(mean, sd), na.rm=T) %>% gather()
n <-nrow(mean_sd)/2
#split mean and SD into different columns
mean_sd <- data.frame(Y=gsub("_mean","",mean_sd[1:n,1]), mean=mean_sd[1:n,2], sd=mean_sd[(n+1):(2*n),2]) 

#Compare to Audrie's

load(here("audrie results/immune_N_means.RData"))
ls()

aud_N <- as.data.frame(rbindlist(lapply(ls(pattern="_N_L"), get)))
aud_N$Y = gsub("_N_L","",ls(pattern="_N_L"))

#merge and compare
N_comp <- merge(aud_N, mean_sd, by="Y")
dim(N_comp)
N_comp$mean.diff <- N_comp$mean.x - N_comp$mean.y
N_comp$sd.diff <- N_comp$sd.x - N_comp$sd.y
max(N_comp$mean.diff)
max(N_comp$sd.diff)


#------------------------------------------------------------------------------------------------
# Unadjusted GLMs
#------------------------------------------------------------------------------------------------

#Update names to match Audrie
colnames(d) <- gsub("gmcsf","gmc",colnames(d))
colnames(d) <- gsub("tnfa","tnf",colnames(d))
colnames(d) <- gsub("ifng","ifn",colnames(d))


#Generate ratio outcomes
# [1] "d23_ln_gmc_unadj_L"          "d23_ln_ifn_unadj_L"          "d23_ln_igf_unadj_L"         
# [4] "d23_ln_il1_unadj_L"          "d23_ln_il10_unadj_L"         "d23_ln_il12_unadj_L"        
# [7] "d23_ln_il13_unadj_L"         "d23_ln_il17_unadj_L"         "d23_ln_il2_unadj_L"         
# [10] "d23_ln_il21_unadj_L"         "d23_ln_il4_unadj_L"          "d23_ln_il5_unadj_L"         
# [13] "d23_ln_il6_unadj_L"          "d23_ln_tnf_unadj_L"          "d23_ratio_gmc_il10_unadj_L" 
# [16] "d23_ratio_ifn_il10_unadj_L"  "d23_ratio_ifn_il13_unadj_L"  "d23_ratio_ifn_il17_unadj_L" 
# [19] "d23_ratio_ifn_il21_unadj_L"  "d23_ratio_ifn_il4_unadj_L"   "d23_ratio_ifn_il5_unadj_L"  
# [22] "d23_ratio_il1_il10_unadj_L"  "d23_ratio_il12_il10_unadj_L" "d23_ratio_il12_il13_unadj_L"
# [25] "d23_ratio_il12_il17_unadj_L" "d23_ratio_il12_il21_unadj_L" "d23_ratio_il12_il4_unadj_L" 
# [28] "d23_ratio_il12_il5_unadj_L"  "d23_ratio_il13_il10_unadj_L" "d23_ratio_il17_il10_unadj_L"
# [31] "d23_ratio_il2_il10_unadj_L"  "d23_ratio_il21_il10_unadj_L" "d23_ratio_il4_il10_unadj_L" 
# [34] "d23_ratio_il5_il10_unadj_L"  "d23_ratio_il6_il10_unadj_L"  "d23_ratio_pro_il10_unadj_L" 
# [37] "d23_ratio_th1_il10_unadj_L"  "d23_ratio_th1_th17_unadj_L"  "d23_ratio_th1_th2_unadj_L"  
# [40] "d23_ratio_th17_il10_unadj_L" "d23_ratio_th2_il10_unadj_L"  "d23_ratio_tnf_il10_unadj_L" 
# 
#             
# [58] "t2_ratio_gmc_il10_unadj_L"   "t2_ratio_ifn_il10_unadj_L"   "t2_ratio_ifn_il13_unadj_L"  
# [61] "t2_ratio_ifn_il17_unadj_L"   "t2_ratio_ifn_il21_unadj_L"   "t2_ratio_ifn_il4_unadj_L"   
# [64] "t2_ratio_ifn_il5_unadj_L"    "t2_ratio_il1_il10_unadj_L"   "t2_ratio_il12_il10_unadj_L" 
# [67] "t2_ratio_il12_il13_unadj_L"  "t2_ratio_il12_il17_unadj_L"  "t2_ratio_il12_il21_unadj_L" 
# [70] "t2_ratio_il12_il4_unadj_L"   "t2_ratio_il12_il5_unadj_L"   "t2_ratio_il13_il10_unadj_L" 
# [73] "t2_ratio_il17_il10_unadj_L"  "t2_ratio_il2_il10_unadj_L"   "t2_ratio_il21_il10_unadj_L" 
# [76] "t2_ratio_il4_il10_unadj_L"   "t2_ratio_il5_il10_unadj_L"   "t2_ratio_il6_il10_unadj_L"  
# [79] "t2_ratio_pro_il10_unadj_L"   "t2_ratio_th1_il10_unadj_L"   "t2_ratio_th1_th17_unadj_L"  
# [82] "t2_ratio_th1_th2_unadj_L"    "t2_ratio_th17_il10_unadj_L"  "t2_ratio_th2_il10_unadj_L"  
# 
#            
# [100] "t3_ratio_gmc_il10_unadj_L"   "t3_ratio_ifn_il10_unadj_L"   "t3_ratio_ifn_il13_unadj_L"  
# [103] "t3_ratio_ifn_il17_unadj_L"   "t3_ratio_ifn_il21_unadj_L"   "t3_ratio_ifn_il4_unadj_L"   
# [106] "t3_ratio_ifn_il5_unadj_L"    "t3_ratio_il1_il10_unadj_L"   "t3_ratio_il12_il10_unadj_L" 
# [109] "t3_ratio_il12_il13_unadj_L"  "t3_ratio_il12_il17_unadj_L"  "t3_ratio_il12_il21_unadj_L" 
# [112] "t3_ratio_il12_il4_unadj_L"   "t3_ratio_il12_il5_unadj_L"   "t3_ratio_il13_il10_unadj_L" 
# [115] "t3_ratio_il17_il10_unadj_L"  "t3_ratio_il2_il10_unadj_L"   "t3_ratio_il21_il10_unadj_L" 
# [118] "t3_ratio_il4_il10_unadj_L"   "t3_ratio_il5_il10_unadj_L"   "t3_ratio_il6_il10_unadj_L"  
# [121] "t3_ratio_pro_il10_unadj_L"   "t3_ratio_th1_il10_unadj_L"   "t3_ratio_th1_th17_unadj_L"  
# [124] "t3_ratio_th1_th2_unadj_L"    "t3_ratio_th17_il10_unadj_L"  "t3_ratio_th2_il10_unadj_L"  
# [127] "t3_ratio_tnf_il10_unadj_L"   "t3_tnf_unadj_L"  

d <- d %>% rename(t2_igf=igf_t2,          t3_igf=igf_t3,          t2_crp=crp_t2,          t2_agp=agp_t2,          t2_gmc=gmc_t2,       
                t2_ifn=ifn_t2,         t2_il10=il10_t2,         t2_il12=il12_t2,         t2_il13=il13_t2,         t2_il17=il17_t2,        
                t2_il1=il1_t2,          t2_il2=il2_t2,          t2_il21=il21_t2,         t2_il4=il4_t2,          t2_il5=il5_t2,         
                t2_il6=il6_t2,          t2_tnf=tnf_t2,         t3_gmc=gmc_t3,        t3_ifn=ifn_t3,         t3_il10=il10_t3,        
                t3_il12=il12_t3,         t3_il13=il13_t3,         t3_il17=il17_t3,         t3_il1=il1_t3,          t3_il2=il2_t3,         
                t3_il21=il21_t3,         t3_il4=il4_t3,          t3_il5=il5_t3,          t3_il6=il6_t3,          t3_tnf=tnf_t3) %>%
           mutate(
                t2_ratio_gmc_il10   =   t2_gmc/t2_il10,   
                t2_ratio_ifn_il10   =   t2_ifn/t2_il10,   
                t2_ratio_ifn_il13  =   t2_ifn/t2_il13,  
                t2_ratio_ifn_il17   =   t2_ifn/t2_il17,   
                t2_ratio_ifn_il21   =   t2_ifn/t2_il21,   
                t2_ratio_ifn_il4  =   t2_ifn/t2_il4,  
                t2_ratio_ifn_il5   =   t2_ifn/t2_il5,   
                t2_ratio_il1_il10   =   t2_il1/t2_il10,   
                t2_ratio_il12_il10 =   t2_il12/t2_il10, 
                t2_ratio_il12_il13  =   t2_il12/t2_il13,  
                t2_ratio_il12_il17  =   t2_il12/t2_il17,  
                t2_ratio_il12_il21 =   t2_il12/t2_il21, 
                t2_ratio_il12_il4   =   t2_il12/t2_il4,   
                t2_ratio_il12_il5   =   t2_il12/t2_il5,   
                t2_ratio_il13_il10 =   t2_il13/t2_il10, 
                t2_ratio_il17_il10  =   t2_il17/t2_il10,  
                t2_ratio_il2_il10   =   t2_il2/t2_il10,   
                t2_ratio_il21_il10 =   t2_il21/t2_il10, 
                t2_ratio_il4_il10   =   t2_il4/t2_il10,   
                t2_ratio_il5_il10   =   t2_il5/t2_il10,   
                t2_ratio_il6_il10  =   t2_il6/t2_il10,  
                #t2_ratio_pro_il10   =   t2_pro/t2_il10,   
                # t2_ratio_th1_il10   =   t2_th1/t2_il10,   
                # t2_ratio_th1_th17  =   t2_th1/t2_th17,  
                # t2_ratio_th1_th2   =   t2_th1/t2_th2,   
                # t2_ratio_th17_il10  =   t2_th17/t2_il10,  
                # t2_ratio_th2_il10  =   t2_th2/t2_il10,  
                # t2_ratio_tnf_il10 =   t2_tnf/t2_il10,
                
                t3_ratio_gmc_il10   =   t3_gmc/t3_il10,   
                t3_ratio_ifn_il10   =   t3_ifn/t3_il10,   
                t3_ratio_ifn_il13  =   t3_ifn/t3_il13,  
                t3_ratio_ifn_il17   =   t3_ifn/t3_il17,   
                t3_ratio_ifn_il21   =   t3_ifn/t3_il21,   
                t3_ratio_ifn_il4  =   t3_ifn/t3_il4,  
                t3_ratio_ifn_il5   =   t3_ifn/t3_il5,   
                t3_ratio_il1_il10   =   t3_il1/t3_il10,   
                t3_ratio_il12_il10 =   t3_il12/t3_il10, 
                t3_ratio_il12_il13  =   t3_il12/t3_il13,  
                t3_ratio_il12_il17  =   t3_il12/t3_il17,  
                t3_ratio_il12_il21 =   t3_il12/t3_il21, 
                t3_ratio_il12_il4   =   t3_il12/t3_il4,   
                t3_ratio_il12_il5   =   t3_il12/t3_il5,   
                t3_ratio_il13_il10 =   t3_il13/t3_il10, 
                t3_ratio_il17_il10  =   t3_il17/t3_il10,  
                t3_ratio_il2_il10   =   t3_il2/t3_il10,   
                t3_ratio_il21_il10 =   t3_il21/t3_il10, 
                t3_ratio_il4_il10   =   t3_il4/t3_il10,   
                t3_ratio_il5_il10   =   t3_il5/t3_il10,   
                t3_ratio_il6_il10  =   t3_il6/t3_il10#,  
                #t3_ratio_pro_il10   =   t3_pro/t3_il10,   
                # t3_ratio_th1_il10   =   t3_th1/t3_il10,   
                # t3_ratio_th1_th17  =   t3_th1/t3_th17,  
                # t3_ratio_th1_th2   =   t3_th1/t3_th2,   
                # t3_ratio_th17_il10  =   t3_th17/t3_il10,  
                # t3_ratio_th2_il10  =   t3_th2/t3_il10,  
                # t3_ratio_tnf_il10 =   t3_tnf/t3_il10
                )



#dataframe of urine biomarkers:
colnames(d)
Y<-d %>% select(t2_igf,          t3_igf,          t2_crp,          t2_agp,          t2_gmc,       
         t2_ifn,         t2_il10,         t2_il12,         t2_il13,         t2_il17,        
         t2_il1,          t2_il2,          t2_il21,         t2_il4,          t2_il5,         
         t2_il6,          t2_tnf,         t3_gmc,        t3_ifn,         t3_il10,        
         t3_il12,         t3_il13,         t3_il17,         t3_il1,          t3_il2,         
         t3_il21,         t3_il4,          t3_il5,          t3_il6,          t3_tnf ,
         t2_ratio_gmc_il10,
         t2_ratio_ifn_il10,
         t2_ratio_ifn_il13,
         t2_ratio_ifn_il17,
         t2_ratio_ifn_il21,
         t2_ratio_ifn_il4,
         t2_ratio_ifn_il5,
         t2_ratio_il1_il10,
         t2_ratio_il12_il10,
         t2_ratio_il12_il13,
         t2_ratio_il12_il17,
         t2_ratio_il12_il21,
         t2_ratio_il12_il4,
         t2_ratio_il12_il5,
         t2_ratio_il13_il10,
         t2_ratio_il17_il10,
         t2_ratio_il2_il10,
         t2_ratio_il21_il10,
         t2_ratio_il4_il10,
         t2_ratio_il5_il10,
         t2_ratio_il6_il10,
         # t2_ratio_pro_il10,
         # t2_ratio_th1_il10,
         # t2_ratio_th1_th17,
         # t2_ratio_th1_th2,
         # t2_ratio_th17_il10,
         # t2_ratio_th2_il10,
         # t2_ratio_tnf_il10,
         
         t3_ratio_gmc_il10,
         t3_ratio_ifn_il10,
         t3_ratio_ifn_il13,
         t3_ratio_ifn_il17,
         t3_ratio_ifn_il21,
         t3_ratio_ifn_il4,
         t3_ratio_ifn_il5,
         t3_ratio_il1_il10,
         t3_ratio_il12_il10,
         t3_ratio_il12_il13,
         t3_ratio_il12_il17,
         t3_ratio_il12_il21,
         t3_ratio_il12_il4,
         t3_ratio_il12_il5,
         t3_ratio_il13_il10,
         t3_ratio_il17_il10,
         t3_ratio_il2_il10,
         t3_ratio_il21_il10,
         t3_ratio_il4_il10,
         t3_ratio_il5_il10,
         t3_ratio_il6_il10#,
         # t3_ratio_pro_il10,
         # t3_ratio_th1_il10,
         # t3_ratio_th1_th17,
         # t3_ratio_th1_th2,
         # t3_ratio_th17_il10,
         # t3_ratio_th2_il10,
         # t3_ratio_tnf_il10
         )

#Replace inf with NA
Y <- do.call(data.frame,lapply(Y, function(x) replace(x, is.infinite(x),NA)))


#Unadjusted glm models
res_unadj <- NULL
for(i in 1:ncol(Y)){
    temp<-washb_glm(Y=(Y[,i]), tr=d$tr, W=NULL, id=d$block, pair=NULL, family="gaussian", contrast= c("Control","Nutrition + WSH"), print=F)
    res_unadj<-rbind(res_unadj, as.numeric(temp$TR))
}
res_unadj <- as.data.frame(res_unadj)

colnames(res_unadj)<-c("RD","ci.l","ci.u", "Std. Error", "z value", "Pval")
res_unadj$Y <-colnames(Y)


#Compare to Audrie's objects
load(here("audrie results/immune_unadj_glm.RData"))


#Function to load and compile Audrie's objects
load_aud <- function(name.pattern, object_list, subgroup=F){
  print(object_list)
  df <- lapply(object_list, get)
  for(i in which(sapply(df, is.null))){
    print(object_list[i])
    if(subgroup){
      df[[i]] <- data.frame(subgroup=rep(NA,2), RD=rep(NA,2), ci.lb=rep(NA,2), ci.ub=rep(NA,2), SE=rep(NA,2), z=rep(NA,2), `P-value`=rep(NA,2))
    }else{
      df[[i]] <- data.frame(RD=NA, ci.lb=NA, ci.ub=NA, SE=NA, z=NA, `P-value`=NA)
    }
  }
  if(subgroup){
    df <- data.frame(rbindlist(lapply(df, as.data.frame), fill=TRUE),Y = rep(gsub(name.pattern,"",object_list), each=2))
  }else{
    df <- data.frame(rbindlist(lapply(df, as.data.frame), fill=TRUE),Y = gsub(name.pattern,"",object_list))
  }
  return(df)
}

name.pattern="_unadj_L"
object_list=ls(pattern=name.pattern)
aud_unadj <- load_aud(name.pattern, object_list)

dim(res_unadj)
dim(aud_unadj)
comp_unadj <- full_join(res_unadj, aud_unadj, by="Y")
dim(comp_unadj)

comp_unadj$RD.x - comp_unadj$RD.y


#------------------------------------------------------------------------------------------------
# Age and sex adjusted GLMs
#------------------------------------------------------------------------------------------------

d$sex<-as.factor(d$sex)
d$sex=relevel(d$sex,ref="0")

#Age and sex adjusted glm models
res_sex <- NULL
for(i in 1:ncol(Y)){
  if(grepl("t2_", colnames(Y)[i])){
    temp<-washb_glm(Y=(Y[,i]), tr=d$tr, W=cbind(d$sex, d$agem2), id=d$block, pair=NULL, family="gaussian", contrast= c("Control","Nutrition + WSH"), print=F)
  }else{
    temp<-washb_glm(Y=(Y[,i]), tr=d$tr, W=cbind(d$sex, d$agem3), id=d$block, pair=NULL, family="gaussian", contrast= c("Control","Nutrition + WSH"), print=F)
  }
  res_sex<-rbind(res_sex, as.numeric(temp$TR))
}
res_sex <- as.data.frame(res_sex)

colnames(res_sex)<-c("RD","ci.l","ci.u", "Std. Error", "z value", "Pval")
res_sex$Y <-colnames(Y)

#Compare to Audrie's objects
load(here("audrie results/immune_adj_sex_age_glm.RData"))

aud_sex <- as.data.frame(rbindlist(lapply(lapply(ls(pattern="_adj_sex_age_L"), get), as.data.frame)))
aud_sex$Y = gsub("_adj_sex_age_L","",ls(pattern="_adj_sex_age_L"))

dim(res_sex)
dim(aud_sex)
comp_sex <- full_join(res_sex, aud_sex, by="Y")
dim(comp_sex)

comp_sex$RD.x - comp_sex$RD.y


#------------------
#Adjusted GLM
#------------------

#Set birthorder to 1, >=2, or missing
# class(d$birthord)
# d$birthord[d$birthord>1]<-"2+"
# d$birthord[is.na(d$birthord)]<-"missing"
# d$birthord<-factor(d$birthord)

#Make vectors of adjustment variable names
Wvars<-c('sex', 'birthord',
         'momage', 'momheight','momedu','hfiacat',
         'Nlt18','Ncomp','watmin',
         'walls', 'floor',
         'elec', 'asset_wardrobe', 'asset_table', 'asset_chair', 'asset_clock', 
         'asset_khat', 'asset_chouki', 'asset_radio', 
         'asset_tv', 'asset_refrig', 'asset_bike',
         'asset_moto', 'asset_sewmach', 'asset_mobile',
         'n_cows', 'n_goats', 'n_chickens')



#Add in time varying covariates:
d <- d %>% mutate(monsoon2 = ifelse(month2 > 4 & month2 < 11, "1", "0"),
                  monsoon3 = ifelse(month3 > 4 & month3 < 11, "1", "0"),
                  monsoon2 = ifelse(is.na(month2),"missing", monsoon2),
                  monsoon3 = ifelse(is.na(month3),"missing", monsoon3),
                  monsoon2 = factor(monsoon2),
                  monsoon3 = factor(monsoon3))

Wvars2<-c("aged2", "monsoon2") 
Wvars3<-c("aged3", "monsoon3") 


#subset time-constant W adjustment set
W<- subset(d, select=Wvars)

#Clean adjustment variables 
#Check missingness
for(i in 1:ncol(W)){
  print(colnames(W)[i])
  print(table(is.na(W[,i])))
}

#Replace missingness for factors with new level
#in main dataset 


d$asset_clock[is.na(d$asset_clock)]<-99
d$asset_clock<-factor(d$asset_clock)

#Order data to replicate SL
d <- d[order(d$dataid,d$childNo, d$svy),]



#Re-subset W so new missing categories are included
W<- subset(d, select=Wvars)

#check that all the factor variables are set
for(i in 1:ncol(W)){
  print(colnames(W)[i])
  print(class(W[,i])  )
}

#Truncate unrealistic levels of n_chickens to 60
table(d$n_chickens)
d$n_chickens[d$n_chickens>60]<-60
table(d$n_chickens)

#Relevel all factors
W$sex<-as.factor(W$sex)
d$sex=relevel(d$sex,ref="0")
d$momedu=relevel(d$momedu,ref="No education")
d$hfiacat=relevel(d$hfiacat,ref="Food Secure")
d$hfiacat<-addNA(d$hfiacat)
d$wall<-factor(d$wall)
d$wall<-addNA(d$wall)
levels(d$wall)<-c("No improved wall","Improved wall","Missing")
d$wall=relevel(d$wall,ref="No improved wall")
d$floor<-factor(d$floor)
d$floor<-addNA(d$floor)
levels(d$floor)<-c("No improved floor","Improved floor","Missing")
d$floor=relevel(d$floor,ref="No improved floor")
d$elec<-factor(d$elec)
d$elec<-addNA(d$elec)
levels(d$elec)<-c("No electricity","Electricity","Missing")
d$elec=relevel(d$elec,ref="No electricity")
d$asset_wardrobe<-factor(d$asset_wardrobe)
d$asset_wardrobe<-addNA(d$asset_wardrobe)
levels(d$asset_wardrobe)<-c("No wardrobe","Wardrobe","Missing")
d$asset_wardrobe=relevel(d$asset_wardrobe,ref="No wardrobe")
d$asset_table<-factor(d$asset_table)
d$asset_table<-addNA(d$asset_table)
levels(d$asset_table)<-c("No table","Improved table","Missing")
d$asset_table=relevel(d$asset_table,ref="No table")
d$asset_chair<-factor(d$asset_chair)
d$asset_chair<-addNA(d$asset_chair)
levels(d$asset_chair)<-c("No chair","Chair","Missing")
d$asset_chair=relevel(d$asset_chair,ref="No chair")
d$asset_clock<-factor(d$asset_clock)
d$asset_clock=relevel(d$asset_clock,ref="No clock")
levels(d$asset_clock)<-c("No clock","Clock","Missing", "Missing")
d$asset_clock=relevel(d$asset_clock,ref="No clock")
d$asset_khat<-factor(d$asset_khat)
d$asset_khat<-addNA(d$asset_khat)
levels(d$asset_khat)<-c("No khat","Khat","Missing")
d$asset_khat=relevel(d$asset_khat,ref="No khat")
d$asset_chouki<-factor(d$asset_chouki)
d$asset_chouki<-addNA(d$asset_chouki)
levels(d$asset_chouki)<-c("No chouki","Chouki","Missing")
d$asset_chouki=relevel(d$asset_chouki,ref="No chouki")
d$asset_tv<-factor(d$asset_tv)
d$asset_tv<-addNA(d$asset_tv)
levels(d$asset_tv)<-c("No TV","Improved TV","Missing")
d$asset_tv=relevel(d$asset_tv,ref="No TV")
d$asset_refrig<-factor(d$asset_refrig)
d$asset_refrig<-addNA(d$asset_refrig)
levels(d$asset_refrig)<-c("No refrigerator","Refrigerator","Missing")
d$asset_refrig=relevel(d$asset_refrig,ref="No refrigerator")
d$asset_bike<-factor(d$asset_bike)
d$asset_bike<-addNA(d$asset_bike)
levels(d$asset_bike)<-c("No bicycle","Bicycle","Missing")
d$asset_bike=relevel(d$asset_bike,ref="No bicycle")
d$asset_moto<-factor(d$asset_moto)
d$asset_moto<-addNA(d$asset_moto)
levels(d$asset_moto)<-c("No motorcycle","Motorcycle","Missing")
d$asset_moto=relevel(d$asset_moto,ref="No motorcycle")
d$asset_sewmach<-factor(d$asset_sewmach)
d$asset_sewmach<-addNA(d$asset_sewmach)
levels(d$asset_sewmach)<-c("No sewing machine","Sewing machine","Missing")
d$asset_sewmach=relevel(d$asset_sewmach,ref="No sewing machine")
d$asset_mobile<-factor(d$asset_mobile)
d$asset_mobile<-addNA(d$asset_mobile)
levels(d$asset_mobile)<-c("No mobile phone","Mobile phone","Missing")
d$asset_mobile=relevel(d$asset_mobile,ref="No mobile phone")    

#Re-subset W so new re-leveled factors are included
W<- subset(d, select=Wvars)




#Check that prevalence >5% for all binary variables
for(i in 1:ncol(W)){
  if(class(W[,i])=="factor"){
    for(j in 1:dim(table(W[,i]))){
      flag<-0
      if(sum(W[,i]==levels(W[,i])[j], na.rm=T)/nrow(W)*100<5){
        perc<-sum(W[,i]==levels(W[,i])[j], na.rm=T)/nrow(W)*100
        cat("\n>95% missing: ",colnames(W)[i]," level:",levels(W[,i])[j],"perc:",perc,"\n")
        flag<-1
      }
    }
    if(flag==1){
      print(table(W[,i]))
    }
  }else{
    if(sum(is.na(W[,i]))/nrow(W)*100>95){
      cat("\n>95% missing: ",colnames(W)[i],"\n")
    }
  }
}



#Add in time-varying covariates
W2<- cbind(W, subset(d, select=Wvars2))
W3<- cbind(W, subset(d, select=Wvars3))




#Tabulate missingness
for(i in 1:ncol(W)){
  print(colnames(W)[i])
  print(table(is.na(W[,i])))
}


#Print means for continious, Ns for factors
for(i in 1:ncol(W)){
  print(colnames(W)[i])
  if(class(W[,i])=="factor"){
    print(table(W[,i]))
  }else{print(mean(W[,i], na.rm=T))}
}



for(i in 1:ncol(W3)){
  print(colnames(W3)[i])
  if(class(W3[,i])=="factor"){
    print(table(W3[,i]))
  }else{print(mean(W3[,i], na.rm=T))}
}




##############################################
#Run GLMs for the adjusted parameter estimates
##############################################


#dataframe of urine biomarkers:
colnames(d)
Y<-d %>% select(t2_igf,          t3_igf,          t2_crp,          t2_agp,          t2_gmc,       
                t2_ifn,         t2_il10,         t2_il12,         t2_il13,         t2_il17,        
                t2_il1,          t2_il2,          t2_il21,         t2_il4,          t2_il5,         
                t2_il6,          t2_tnf,         t3_gmc,        t3_ifn,         t3_il10,        
                t3_il12,         t3_il13,         t3_il17,         t3_il1,          t3_il2,         
                t3_il21,         t3_il4,          t3_il5,          t3_il6,          t3_tnf ,
                t2_ratio_gmc_il10,
                t2_ratio_ifn_il10,
                t2_ratio_ifn_il13,
                t2_ratio_ifn_il17,
                t2_ratio_ifn_il21,
                t2_ratio_ifn_il4,
                t2_ratio_ifn_il5,
                t2_ratio_il1_il10,
                t2_ratio_il12_il10,
                t2_ratio_il12_il13,
                t2_ratio_il12_il17,
                t2_ratio_il12_il21,
                t2_ratio_il12_il4,
                t2_ratio_il12_il5,
                t2_ratio_il13_il10,
                t2_ratio_il17_il10,
                t2_ratio_il2_il10,
                t2_ratio_il21_il10,
                t2_ratio_il4_il10,
                t2_ratio_il5_il10,
                t2_ratio_il6_il10,
                # t2_ratio_pro_il10,
                # t2_ratio_th1_il10,
                # t2_ratio_th1_th17,
                # t2_ratio_th1_th2,
                # t2_ratio_th17_il10,
                # t2_ratio_th2_il10,
                # t2_ratio_tnf_il10,
                
                t3_ratio_gmc_il10,
                t3_ratio_ifn_il10,
                t3_ratio_ifn_il13,
                t3_ratio_ifn_il17,
                t3_ratio_ifn_il21,
                t3_ratio_ifn_il4,
                t3_ratio_ifn_il5,
                t3_ratio_il1_il10,
                t3_ratio_il12_il10,
                t3_ratio_il12_il13,
                t3_ratio_il12_il17,
                t3_ratio_il12_il21,
                t3_ratio_il12_il4,
                t3_ratio_il12_il5,
                t3_ratio_il13_il10,
                t3_ratio_il17_il10,
                t3_ratio_il2_il10,
                t3_ratio_il21_il10,
                t3_ratio_il4_il10,
                t3_ratio_il5_il10,
                t3_ratio_il6_il10#,
                # t3_ratio_pro_il10,
                # t3_ratio_th1_il10,
                # t3_ratio_th1_th17,
                # t3_ratio_th1_th2,
                # t3_ratio_th17_il10,
                # t3_ratio_th2_il10,
                # t3_ratio_tnf_il10
)

#Replace inf with NA
Y <- do.call(data.frame,lapply(Y, function(x) replace(x, is.infinite(x),NA)))

#Fully adjusted glm models
res_adj <- NULL
for(i in 1:ncol(Y)){
  if(grepl("t2_", colnames(Y)[i])){
    temp<-washb_glm(Y=(Y[,i]), tr=d$tr, W=W2, id=d$block, pair=NULL, family="gaussian", contrast= c("Control","Nutrition + WSH"), print=F)
  }else{
    temp<-washb_glm(Y=(Y[,i]), tr=d$tr, W=W3, id=d$block, pair=NULL, family="gaussian", contrast= c("Control","Nutrition + WSH"), print=F)
  }
  res_adj<-rbind(res_adj, as.numeric(temp$TR))
}
res_adj <- as.data.frame(res_adj)

colnames(res_adj)<-c("RD","ci.l","ci.u", "Std. Error", "z value", "Pval")
res_adj$Y <-colnames(Y)

#Compare to Audrie's objects
load(here("audrie results/immune_adj_glm.RData"))
name.pattern="_adj_L"
object_list=ls(pattern=name.pattern)
aud_adj <- load_aud(name.pattern, object_list)


dim(res_adj)
dim(aud_adj)
comp_adj <- full_join(res_adj, aud_adj, by="Y")
dim(comp_adj)
comp_adj$RD.x - comp_adj$RD.y
comp_adj$Pval - comp_adj$P.value

#Save intermediate R objects for replication comparison
dm <- d
save(res_adj, W, W2, W3, dm, comp_adj, Y, file = here("replication objects/andrew_immune_W.rdata"))


 
##############################################
#Run GLMs for the sex-stratified subgroup analysis
##############################################

#sex stratified glm models
res_sub <- NULL
for(i in 1:ncol(Y)){
  temp<-washb_glm(Y=(Y[,i]), tr=d$tr, W=data.frame(sex=d$sex), V="sex", id=d$block, pair=NULL, family="gaussian", contrast= c("Control","Nutrition + WSH"), print=F)
  res_sub<-rbind(res_sub, temp$lincom)
}
res_sub <- as.data.frame(res_sub)

colnames(res_sub)<-c("sex","RD","ci.l","ci.u", "Std. Error", "z value", "Pval")
res_sub$Y <-rep(colnames(Y), each=2)
res_sub <- res_sub %>% mutate(subgroup = case_when(sex==1 ~ "male", sex==0 ~ "female", TRUE~""), subgroup=factor(subgroup))

#Compare to Audrie's objects
load(here("audrie results/immune_subgroup.RData"))

name.pattern="_subgroup_L"
object_list=ls(pattern=name.pattern)
aud_sub <- load_aud(name.pattern, object_list, subgroup = T)

dim(res_sub)
dim(aud_sub)
comp_sub <- full_join(res_sub, aud_sub, by=c("Y","subgroup"))
dim(comp_sub)

comp_sub <- filter(comp_sub, !is.na(RD.x))

comp_sub$RD.x - comp_sub$RD.y

il6_t3_subgroup_L

# ##############################################
# #Plot results
# ##############################################
# 
# res_adj$marker <- str_split(res_adj$Y, "_t", simplify=T)[,1]
# res_adj$round <- as.numeric(str_split(res_adj$Y, "_t", simplify=T)[,2])
# res_adj <- res_adj %>% arrange(round, RD) %>%
#   mutate(Y=factor(Y, levels=unique(Y)))
# 
# ggplot(res_adj, aes(x=Y, y=RD)) + geom_point() +
#   geom_pointrange(aes(ymin=ci.l, ymax=ci.u)) + geom_hline(yintercept = 0) + 
#   facet_wrap(~round, scales="free_x") + theme_bw()
# 
# 
# ##############################################
# #Examine baseline data
# ##############################################
# 
# load(paste0(dropboxDir,"Results/Audrie/immune_enrol_baseline_char.RData"))
# load(paste0(dropboxDir,"Results/Audrie/immune_enrol_supp_baseline_char_t2.RData"))
# load(paste0(dropboxDir,"Results/Audrie/immune_enrol_supp_baseline_char_lost_t3.RData"))
# 
