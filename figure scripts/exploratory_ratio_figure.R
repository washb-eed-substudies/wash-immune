




rm(list=ls())
source(here::here("0-config.R"))

######################
###Load in data
######################

#Load in enrollment data,blinded tr data, stool data for adjusted analysis. Use read.dta() to read the .dta files, or read.csv() to 
#read .csv files. Use stringAsFactors=TRUE so that any character-based variable will be read in as a factor.
d<-readRDS(paste0(dropboxDir,"Data/Cleaned/Audrie/bangladesh-immune-analysis-dataset.rds"))


#function to create composite score
create_numscore <- function(d, numerator_vars=c("il1_t2", "il6_t2", "tnfa_t2"),  varname="t2_ratio_pro_il10"){
  for(i in numerator_vars){
    if(i==i[1]){
      x = as.vector(scale(d[,i], center = FALSE, scale = apply(as.matrix(d[,i]), 2, sd, na.rm = TRUE)))
    }else{
      x = x + as.vector(scale(d[,i], center = FALSE, scale = apply(as.matrix(d[,i]), 2, sd, na.rm = TRUE)))
    }
  }

  d$score <- x
  colnames(d)[ncol(d)] <- varname
  return(d)
}



# *Th1 / IL-10
# *(IL-12 + IFN) / IL-10
# gen t2_ratio_th1_il10 = (il12_t2 + ifng_t2) / il10_t2
d <- create_numscore(d, numerator_vars=c("il12_t2", "ifng_t2"), varname="t2_th1")
summary(d$t2_th1)
ggplot(d, aes(x=t2_th1)) + geom_density()

# *Th2 / IL-10 
# *(IL-4 + IL-5 + IL-13) / IL-10
# gen t2_ratio_th2_il10 = (il4_t2 + il5_t2 + il13_t2) / il10_t2
d <- create_numscore(d, numerator_vars=c("il4_t2", "il5_t2", "il13_t2"), varname="t2_th2")
summary(d$t2_th1)
ggplot(d, aes(x=t2_th1)) + geom_density()


d <- create_numscore(d, numerator_vars=c("il12_t3", "ifng_t3"), varname="t3_th1")
d <- create_numscore(d, numerator_vars=c("il4_t3", "il5_t3", "il13_t3"), varname="t3_th2")

df1 <- d %>% subset(., select = c(tr, t2_th1, t2_th2, il10_t2, t2_ratio_th1_il10, t2_ratio_th2_il10, t2_ratio_th1_th2)) %>%
              mutate(time="t2") %>%
              rename(th1=t2_th1, 
                      th2=t2_th2, 
                      il10=il10_t2, 
                      ratio_th1_il10=t2_ratio_th1_il10, 
                      ratio_th2_il10=t2_ratio_th2_il10,
                      ratio_th1_th2=t2_ratio_th1_th2)
df2 <- d %>% subset(., select = c(tr, t3_th1, t3_th2, il10_t3, t3_ratio_th1_il10, t3_ratio_th2_il10, t3_ratio_th1_th2)) %>%
                mutate(time="t3") %>%
                rename(th1=t3_th1, 
                       th2=t3_th2, 
                       il10=il10_t3, 
                       ratio_th1_il10=t3_ratio_th1_il10, 
                       ratio_th2_il10=t3_ratio_th2_il10,
                       ratio_th1_th2=t3_ratio_th1_th2)


d <- bind_rows(df1, df2) %>% 
  pivot_longer(!c(tr,time), names_to = "biomarker", values_to = "Y") %>%
  filter(!is.na(Y))
head(d)

d <- d %>% group_by(biomarker) %>% mutate(Z=scale(Y))



ggplot(d, aes(x=time, y=Y, group=biomarker, color=biomarker)) + geom_point(alpha=0.1) +
  facet_wrap(~biomarker)

ggplot(d, aes(x=time, y=Y, color=biomarker, group=1)) + geom_point(alpha=0.1) +
  facet_wrap(~biomarker) + geom_point(stat='summary', fun=mean) +
  stat_summary(fun=mean, geom="line")


ggplot(d, aes(x=time, y=Y, color=biomarker, group=biomarker)) + 
  geom_point(stat='summary', fun=mean) +
  stat_summary(fun=mean, geom="line") + 
  theme(legend.position = "right")


ggplot(d, aes(x=time, y=Z, color=biomarker, group=biomarker)) + 
  geom_point(stat='summary', fun=mean) +
  stat_summary(fun=mean, geom="line") + 
  theme(legend.position = "right")







