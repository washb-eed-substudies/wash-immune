




rm(list=ls())
source(here::here("0-config.R"))
library(ggrepel)

######################
###Load in data
######################

#Load in enrollment data,blinded tr data, stool data for adjusted analysis. Use read.dta() to read the .dta files, or read.csv() to 
#read .csv files. Use stringAsFactors=TRUE so that any character-based variable will be read in as a factor.
d<-readRDS(paste0(dropboxDir,"Data/Cleaned/Audrie/bangladesh-immune-analysis-dataset.rds"))


#function to create composite score
create_numscore <- function(d, numerator_vars=c("il1_t2", "il6_t2", "tnfa_t2"),  varname="t2_ratio_pro_il10"){
  for(i in numerator_vars){
    if(i==numerator_vars[1]){
      x = as.vector(scale(d[,i], center = FALSE, scale = apply(as.matrix(d[,i]), 2, sd, na.rm = TRUE)))
    }else{
      x = x + as.vector(scale(d[,i], center = FALSE, scale = apply(as.matrix(d[,i]), 2, sd, na.rm = TRUE)))
    }
  }
  
  d$score <- x
  colnames(d)[ncol(d)] <- varname
  return(d)
}

# d$x1=as.vector(scale(d$il12_t2, center = FALSE, scale = apply(as.matrix(d$il12_t2), 2, sd, na.rm = TRUE))) 
# d$x2=as.vector(scale(d$ifng_t2, center = FALSE, scale = apply(as.matrix(d$ifng_t2), 2, sd, na.rm = TRUE)))
# 
# d %>% group_by(tr) %>% summarize(mean(x1, na.rm=T), mean(x2, na.rm=T), mean(t2_th1, na.rm=T), mean(il12_t2, na.rm=T), mean(ifng_t2, na.rm=T))



d <- create_numscore(d, numerator_vars=c("il12_t2", "ifng_t2"), varname="t2_th1")
d <- create_numscore(d, numerator_vars=c("il4_t2", "il5_t2", "il13_t2"), varname="t2_th2")
d <- create_numscore(d, numerator_vars=c("il12_t3", "ifng_t3"), varname="t3_th1")
d <- create_numscore(d, numerator_vars=c("il4_t3", "il5_t3", "il13_t3"), varname="t3_th2")


df1 <- d %>% subset(., select = c(childid, tr, t2_th1, t2_th2,  t2_ratio_th1_th2, il12_t2, ifng_t2, il4_t2, il5_t2, il13_t2)) %>%
  mutate(time="t2") %>%
  rename(th1=t2_th1, 
         th2=t2_th2, 
         il12=il12_t2, 
         ifng=ifng_t2, 
         il4=il4_t2, 
         il5=il5_t2, 
         il13=il13_t2,
         ratio_th1_th2=t2_ratio_th1_th2)
df2 <- d %>% subset(., select = c(childid, tr, t3_th1, t3_th2, t3_ratio_th1_th2, il12_t3, ifng_t3, il4_t3, il5_t3, il13_t3)) %>%
  mutate(time="t3") %>%
  rename(th1=t3_th1, 
         th2=t3_th2, 
         il12=il12_t3, 
         ifng=ifng_t3, 
         il4=il4_t3, 
         il5=il5_t3, 
         il13=il13_t3,
         ratio_th1_th2=t3_ratio_th1_th2)

dfull<- bind_rows(df1, df2)
head(dfull)

d <- dfull %>% 
  pivot_longer(!c(childid, tr,time), names_to = "biomarker", values_to = "Y") %>%
  filter(!is.na(Y))
head(d)

d$combined <- ifelse(d$biomarker %in% c("th1","th2","ratio_th1_th2"),1,0)
d$group <- ifelse(d$biomarker %in% c("th1","il12", "ifng"),"th1","th2")
d$group[d$biomarker=="ratio_th1_th2"] <- "ratio"

df <- d %>% filter(biomarker=="il12" & time=="t2")


x = as.vector(scale(df$Y, center = FALSE, scale = apply(as.matrix(df$Y), 2, sd, na.rm = TRUE)))
mean(x, na.rm=T)


d <- d %>% group_by(time, biomarker) %>% 
  mutate(Z=ifelse(combined==1,Y,scale(Y)),
         scale=ifelse(combined==1,Y,scale(Y, center=FALSE, scale = apply(as.matrix(Y), 2, sd, na.rm = TRUE))))
mean(d$scale[d$biomarker=="il12" & d$time=="t2"])

d$scale[d$biomarker=="ratio_th1_th2"] <- exp(d$scale[d$biomarker=="ratio_th1_th2"])

d$label <- d$biomarker
d$label[d$time=="t3"] <- NA

d <- d %>% arrange(childid, time, group)
head(d, 10)

# ggplot(d, aes(x=time, y=scale, color=biomarker, group=biomarker)) +
#   facet_wrap(~group) + geom_point(stat='summary', fun=mean) +
#   stat_summary(fun=mean, geom="line") +
#   geom_text(x="t3", aes( y=scale, label=label))

df <- d %>% ungroup() %>% group_by(time, biomarker, group) %>% 
  summarise(scale=mean(scale), Y=mean(Y), Z=mean(Z))

df$label <- df$biomarker
df$label[df$time=="t3"] <- NA

df <- df %>% arrange(time, group, biomarker)

ggplot(df, aes(x=time, y=scale, color=biomarker, group=biomarker, label=label)) +
  facet_wrap(~group) + geom_point(stat='summary', fun=mean) +
  stat_summary(fun=mean, geom="line") +
  geom_text_repel()



#Make density plots with mean/median as vertical line
#t2 and t3 as groups
#th1 and th2 as left and right columns
#ratio at the bottom




df <- d %>% filter(biomarker!="ratio_th1_th2")
p <- ggplot(data=df,aes(x=scale,group=time,color=time,fill=time)) +
  facet_wrap(biomarker~group,ncol=2, scales ="free") +
  geom_density(aes(y=..density..),color=NA,alpha=0.7)



df$time <- factor(df$time, levels = c("t2","t3"))
df1 <- df %>% filter(group=="th1") %>% group_by(time,biomarker) %>% mutate(meanY=mean(scale,na.rm=T))

p <- ggplot(data=df1,aes(x=scale,group=time,color=time,fill=time)) +
  facet_wrap(~biomarker, ncol=1, scales="fixed") +
  geom_density(aes(y=..density..),color=NA,alpha=0.7) +
  geom_vline(aes(xintercept= meanY)) +
  coord_cartesian(xlim=c(0,8)) +
  theme(legend.position = "right")
p





