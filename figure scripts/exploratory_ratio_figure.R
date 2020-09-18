




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

d$x1=as.vector(scale(d$il12_t2, center = FALSE, scale = apply(as.matrix(d$il12_t2), 2, sd, na.rm = TRUE))) 
d$x2=as.vector(scale(d$ifng_t2, center = FALSE, scale = apply(as.matrix(d$ifng_t2), 2, sd, na.rm = TRUE)))

d %>% group_by(tr) %>% summarize(mean(x1, na.rm=T), mean(x2, na.rm=T), mean(t2_th1, na.rm=T), mean(il12_t2, na.rm=T), mean(ifng_t2, na.rm=T))


table(x1+x2==d$t2_th1)

x1[1]
x2[1]
d$t2_th1[1]

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

# df1 <- d %>% subset(., select = c(tr, t2_th1, t2_th2, il10_t2, t2_ratio_th1_il10, t2_ratio_th2_il10, t2_ratio_th1_th2)) %>%
#               mutate(time="t2") %>%
#               rename(th1=t2_th1, 
#                       th2=t2_th2, 
#                       il10=il10_t2, 
#                       ratio_th1_il10=t2_ratio_th1_il10, 
#                       ratio_th2_il10=t2_ratio_th2_il10,
#                       ratio_th1_th2=t2_ratio_th1_th2)
# df2 <- d %>% subset(., select = c(tr, t3_th1, t3_th2, il10_t3, t3_ratio_th1_il10, t3_ratio_th2_il10, t3_ratio_th1_th2)) %>%
#                 mutate(time="t3") %>%
#                 rename(th1=t3_th1, 
#                        th2=t3_th2, 
#                        il10=il10_t3, 
#                        ratio_th1_il10=t3_ratio_th1_il10, 
#                        ratio_th2_il10=t3_ratio_th2_il10,
#                        ratio_th1_th2=t3_ratio_th1_th2)
# 
# 
# d <- bind_rows(df1, df2) %>% 
#   pivot_longer(!c(tr,time), names_to = "biomarker", values_to = "Y") %>%
#   filter(!is.na(Y))
# head(d)
# 
# d <- d %>% group_by(biomarker) %>% mutate(Z=scale(Y))



# ggplot(d, aes(x=time, y=Y, group=biomarker, color=biomarker)) + geom_point(alpha=0.1) +
#   facet_wrap(~biomarker)
# 
# ggplot(d, aes(x=time, y=Y, color=biomarker, group=1)) + geom_point(alpha=0.1) +
#   facet_wrap(~biomarker) + geom_point(stat='summary', fun=mean) +
#   stat_summary(fun=mean, geom="line")
# 
# 
# ggplot(d, aes(x=time, y=Y, color=biomarker, group=biomarker)) + 
#   geom_point(stat='summary', fun=mean) +
#   stat_summary(fun=mean, geom="line") + 
#   theme(legend.position = "right")
# 
# ggplot(d, aes(x=time, y=Y, color=biomarker, group=biomarker)) + 
#   geom_point(alpha=0.05) +
#   facet_wrap(~biomarker, scales = "free") +
#   geom_point(stat='summary', fun=mean) +
#   stat_summary(fun=mean, geom="line") + 
#   theme(legend.position = "right")
# 
# 
# ggplot(d, aes(x=time, y=Z, color=biomarker, group=biomarker)) + 
#   geom_point(stat='summary', fun=mean) +
#   stat_summary(fun=mean, geom="line") + 
#   theme(legend.position = "right")
# 
# 
# ggplot(d, aes(x=time, y=Y, color=tr, group=tr)) + 
#   #geom_point(alpha=0.05) +
#   facet_wrap(~biomarker, scales = "free") +
#   geom_point(stat='summary', fun=mean) +
#   stat_summary(fun=mean, geom="line") + 
#   theme(legend.position = "right")
# 
# 
# ggplot(d, aes(x=time, y=Z, color=biomarker, group=biomarker)) + 
#   geom_point(stat='summary', fun=mean) +
#   stat_summary(fun=mean, geom="line") + 
#   theme(legend.position = "right")
# 
# 
# 
# ggplot(d, aes(x=time, y=Z, color=biomarker, group=biomarker)) + 
#   geom_point(stat='summary', fun=mean) +
#   stat_summary(fun=mean, geom="line") + 
#   theme(legend.position = "right") +
#   facet_wrap(~tr)




#Make a th1/th2 focused plot with individual trajectories




df1 <- d %>% subset(., select = c(tr, t2_th1, t2_th2,  t2_ratio_th1_th2, il12_t2, ifng_t2, il4_t2, il5_t2, il13_t2)) %>%
  mutate(time="t2") %>%
  rename(th1=t2_th1, 
         th2=t2_th2, 
         il12=il12_t2, 
         ifng=ifng_t2, 
         il4=il4_t2, 
         il5=il5_t2, 
         il13=il13_t2,
         ratio_th1_th2=t2_ratio_th1_th2)
df2 <- d %>% subset(., select = c(tr, t3_th1, t3_th2, t3_ratio_th1_th2, il12_t3, ifng_t3, il4_t3, il5_t3, il13_t3)) %>%
  mutate(time="t3") %>%
  rename(th1=t3_th1, 
         th2=t3_th2, 
         il12=il12_t3, 
         ifng=ifng_t3, 
         il4=il4_t3, 
         il5=il5_t3, 
         il13=il13_t3,
         ratio_th1_th2=t3_ratio_th1_th2)


d <- bind_rows(df1, df2) %>% 
  pivot_longer(!c(tr,time), names_to = "biomarker", values_to = "Y") %>%
  filter(!is.na(Y))
head(d)

d <- d %>% group_by(biomarker) %>% mutate(Z=scale(Y))


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

ggplot(d, aes(x=time, y=Y, color=biomarker, group=biomarker)) + 
  geom_point(alpha=0.05) +
  facet_wrap(~biomarker, scales = "free") +
  geom_point(stat='summary', fun=mean) +
  stat_summary(fun=mean, geom="line") + 
  theme(legend.position = "right")


ggplot(d, aes(x=time, y=Z, color=biomarker, group=biomarker)) + 
  geom_point(stat='summary', fun=mean) +
  stat_summary(fun=mean, geom="line") + 
  theme(legend.position = "right")


ggplot(d, aes(x=time, y=Y, color=tr, group=tr)) + 
  #geom_point(alpha=0.05) +
  facet_wrap(~biomarker, scales = "free") +
  geom_point(stat='summary', fun=mean) +
  stat_summary(fun=mean, geom="line") + 
  theme(legend.position = "right")


ggplot(d, aes(x=time, y=Z, color=biomarker, group=biomarker)) + 
  geom_point(stat='summary', fun=mean) +
  stat_summary(fun=mean, geom="line") + 
  theme(legend.position = "right")



ggplot(d, aes(x=time, y=Z, color=biomarker, group=biomarker)) + 
  geom_point(stat='summary', fun=mean) +
  stat_summary(fun=mean, geom="line") + 
  theme(legend.position = "right") +
  facet_wrap(~tr)


ggplot(d, aes(x=time, y=Z, color=tr, group=tr)) + 
  geom_point(stat='summary', fun=mean) +
  stat_summary(fun=mean, geom="line") + 
  theme(legend.position = "right") +
  facet_wrap(~biomarker)




ggplot(d[d$biomarker %in% c("il12", "ifng", "th1"),], aes(x=time, y=Y, color=tr, group=tr)) + 
  geom_point(stat='summary', fun=mean) +
  stat_summary(fun=mean, geom="line") + 
  theme(legend.position = "right") +
  facet_wrap(~biomarker)

ggplot(d[d$biomarker %in% c("il4", "il5", "il13", "th1"),], aes(x=time, y=Z, color=tr, group=tr)) + 
  geom_point(stat='summary', fun=mean) +
  stat_summary(fun=mean, geom="line") + 
  theme(legend.position = "right") +
  facet_wrap(~biomarker)





df <- d 
df$Z[df$biomarker %in% c("th1","th2")] <- df$Y[df$biomarker %in% c("th1","th2")]


ggplot(df[df$biomarker %in% c("il4", "il5", "il13", "th1"),], aes(x=time, y=Z, color=tr, group=tr)) + 
  geom_point(stat='summary', fun=mean) +
  stat_summary(fun=mean, geom="line") + 
  theme(legend.position = "right") +
  facet_wrap(~biomarker)


ggplot(df[df$biomarker %in% c("il4", "il5", "il13", "th1"),], aes(x=time, y=Y, color=tr, group=tr)) + 
  geom_point(stat='summary', fun=mean) +
  stat_summary(fun=mean, geom="line") + 
  theme(legend.position = "right") +
  facet_wrap(~biomarker, scales = "free")


df2<- d %>% group_by(biomarker) %>% do(as.data.frame(scale(.$Y, center = FALSE, scale = apply(as.matrix(.$Y), 2, sd, na.rm = TRUE))))
df$scale <- df2$V1
df$scale[df$biomarker %in% c("th1","th2")] <- df$Y[df$biomarker %in% c("th1","th2")]





ggplot(df[df$biomarker %in% c("il4", "il5", "il13", "th2"),], aes(x=time, y=scale, color=tr, group=tr)) + 
  geom_point(stat='summary', fun=mean) +
  stat_summary(fun=mean, geom="line") + 
  theme(legend.position = "right") +
  facet_wrap(~biomarker, nrow=1)


table(df$scale[df$biomarker=="il4"])




ggplot(df[df$biomarker %in% c("il12", "ifng", "th1"),], aes(x=time, y=scale, color=tr, group=tr)) + 
  geom_point(stat='summary', fun=mean) +
  stat_summary(fun=mean, geom="line") + 
  theme(legend.position = "right") +
  facet_wrap(~biomarker, nrow=1)
