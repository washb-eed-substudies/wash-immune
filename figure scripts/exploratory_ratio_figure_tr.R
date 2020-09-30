




rm(list=ls())
source(here::here("0-config.R"))
library(ggrepel)
library(directlabels)


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

d <- create_numscore(d, numerator_vars=c("il1_t2", "il6_t2", "tnfa_t2"), varname="t2_pro")
d <- create_numscore(d, numerator_vars=c("il1_t2", "il6_t2", "tnfa_t2"), varname="t3_pro")


il10_t3 <- scale(d$il10_t3, center=FALSE, scale = apply(as.matrix(d$il10_t3), 2, sd, na.rm = TRUE))
summary(il10_t3)

df1 <- d %>% subset(., select = c(childid, clusterid, tr, t2_ratio_pro_il10, t2_ratio_th1_il10, t2_ratio_th2_il10,t2_pro, il10_t2, il1_t2, il6_t2, tnfa_t2, t2_th1, t2_th2,  t2_ratio_th1_th2, il12_t2, ifng_t2, il4_t2, il5_t2, il13_t2)) %>%
  mutate(time="14 months") %>%
  rename(Th1=t2_th1, 
         Th2=t2_th2, 
         IL12=il12_t2, 
         `IFN-` =ifng_t2, 
         IL4=il4_t2, 
         IL5=il5_t2, 
         IL13=il13_t2,
         IL1=il1_t2, 
         IL6=il6_t2, 
         `TNF-a`=tnfa_t2,
         IL10=il10_t2,
         Pro=t2_pro,
         ratio_pro_il10=t2_ratio_pro_il10,
         ratio_th1_il10=t2_ratio_th1_il10, 
         ratio_th2_il10=t2_ratio_th2_il10,
         ratio_th1_th2=t2_ratio_th1_th2)
df2 <- d %>% subset(., select = c(childid, clusterid, tr, t3_ratio_pro_il10, t3_ratio_th1_il10, t3_ratio_th2_il10, t3_pro, il10_t3, il1_t3, il6_t3, tnfa_t3, t3_th1, t3_th2, t3_ratio_th1_th2, il12_t3, ifng_t3, il4_t3, il5_t3, il13_t3)) %>%
  mutate(time="28 months") %>%
  rename(Th1=t3_th1, 
         Th2=t3_th2, 
         IL12=il12_t3, 
         `IFN-`=ifng_t3, 
         IL4=il4_t3, 
         IL5=il5_t3, 
         IL13=il13_t3,
         IL1=il1_t3, 
         IL6=il6_t3, 
         `TNF-a`=tnfa_t3,
         IL10=il10_t3,
         Pro=t3_pro,
         ratio_pro_il10=t3_ratio_pro_il10,
         ratio_th1_il10=t3_ratio_th1_il10, 
         ratio_th2_il10=t3_ratio_th2_il10,
         ratio_th1_th2=t3_ratio_th1_th2)

dfull<- bind_rows(df1, df2)
head(dfull)


d <- dfull %>% 
  pivot_longer(!c(childid, clusterid, tr,time), names_to = "biomarker", values_to = "Y") %>%
  filter(!is.na(Y))
head(d)

unique(d$biomarker)

d$combined <- ifelse(d$biomarker %in% c("Th1","Th2","Pro","ratio_th1_th2","ratio_th1_il10","ratio_th2_il10","ratio_pro_il10"),1,0)
d$group <- ifelse(d$biomarker %in% c("Th1","IL12", "IFN-"),"Th1","Th2")
d$group[d$biomarker=="ratio_th1_th2"] <- "Th1/Th2 ratio"
d$group[d$biomarker=="ratio_th1_il10"] <- "Th1/IL10 ratio"
d$group[d$biomarker=="ratio_th2_il10"] <- "Th2/IL10 ratio"
d$group[d$biomarker=="ratio_pro_il10"] <- "Pro/IL10 ratio"
d$group[d$biomarker=="IL10"] <- "IL10"
d$group <- ifelse(d$biomarker %in% c("TNF-a","IL1", "IL6", "Pro"),"Pro-inflammatory",d$group)


#Drop kids with just some of the constituent biomarkers
d <- d %>% group_by(childid, time, tr, group) %>% 
  mutate(N=n()) %>% 
  filter((N==1 & (group=="Th1/Th2 ratio"|group=="Th1/IL10 ratio"|group=="Th2/IL10 ratio"|group=="Pro/IL10 ratio"|group=="IL10"))|(N==3 & group=="Th1")|(N==4 & (group=="Th2"|group=="Pro-inflammatory"))) %>%
  ungroup()

table(d$group)

d <- d %>% group_by(time, biomarker) %>% 
  mutate(Z=ifelse(combined==1,Y,scale(Y)),
         scale=ifelse(combined==1,Y,scale(Y, center=FALSE, scale = apply(as.matrix(Y), 2, sd, na.rm = TRUE))))
mean(d$scale[d$biomarker=="IL-12" & d$time=="14 months"])

d$scale[d$biomarker=="ratio_th1_th2"] <- exp(d$scale[d$biomarker=="ratio_th1_th2"])
d$scale[d$biomarker=="ratio_th1_il10"] <- exp(d$scale[d$biomarker=="ratio_th1_il10"])
d$scale[d$biomarker=="ratio_th2_il10"] <- exp(d$scale[d$biomarker=="ratio_th2_il10"])
d$scale[d$biomarker=="ratio_pro_il10"] <- exp(d$scale[d$biomarker=="ratio_pro_il10"])
d$biomarker[d$biomarker=="ratio_th1_th2"] <- "Th1/Th2\nratio"


# df <- d %>% ungroup() %>% 
#   group_by(time, tr, biomarker, group, combined) %>% 
#   summarise(scale=mean(scale), med_scale=median(scale), Y=mean(Y), Z=mean(Z))
# head(df)


df<-d %>% group_by(biomarker, group, combined, time, tr) %>%
  do(as.data.frame(washb_mean(Y=.$scale, id=.$clusterid, print = F))) %>%
  rename(scale=Mean, ci.lb=`Lower 95%CI`, ci.ub=`Upper 95%CI`) %>%
  mutate(biomarker_tr = paste0(biomarker,tr))


df$label2 <- df$label1 <- df$biomarker
df$label2[df$biomarker=="IFN-"] <- df$label1[df$biomarker=="IFN-"]  <- "IFN-g"
#df$label2[df$biomarker=="IFN-"] <- df$label1[df$biomarker=="IFN-"]  <- "IFN-\u03B3"
df$label1[!(df$time=="14 months" & df$tr=="Control")] <- NA
df$label2[!(df$time=="28 months" & df$tr=="Nutrition + WSH")] <- NA

df$tr <- as.character(df$tr)

#duplicate the th1-th2 trajectories for the ratio panel
df2 <- df %>% filter(biomarker=="Th1"|biomarker=="Th2") %>%
  mutate(group="Th1/Th2 ratio", combined=0) 
  
df <- bind_rows(df, df2)
df$combined <- factor(df$combined)
df$group <- factor(df$group, levels=c("Th1","Th2","Th1/Th2 ratio","Pro-inflammatory","IL10", "Th1/IL10 ratio", "Th2/IL10 ratio", "Pro/IL10 ratio"))



#---------------------------------------------
# Th1/Th2 ratio
#---------------------------------------------

plotdf <- df[df$group %in% c("Th1", "Th2", "Th1/Th2 ratio"),] %>% droplevels()

dodge_width=0.4
p <- ggplot(plotdf, aes(x=time, y=scale, color=tr, group=biomarker_tr, shape=combined, type=combined)) +
  facet_wrap(~group) + 
  geom_point(position=position_dodge(dodge_width)) +
  geom_line(aes(linetype = combined), show.legend=FALSE,  position=position_dodge(dodge_width)) +
  geom_linerange(aes(ymin=ci.lb, ymax=ci.ub), position=position_dodge(dodge_width)) +
  #stat_summary(fun=mean, geom="line", aes(linetype = combined), show.legend=FALSE) + 
  scale_color_manual(name = "Treatment\nArm", values = tableau10) +
  scale_shape_manual(values = c(1, 19), guide=F) +
  scale_linetype_manual(values = c("dashed", "solid")) +
  theme(legend.position = "bottom") +
  #geom_label_repel(nudge_x = 1,
    #nudge_x=.1, nudge_y=-0.1, 
    #show.legend = FALSE) +
  #geom_text_repel(aes(label = label1), position=position_dodge(dodge_width), point.padding=0.2, cex = 2.5, segment.color = NA) +
  geom_dl(aes(label = label1),  position=position_dodge(dodge_width), method = list(dl.trans(x = x - 0.3), "first.points", cex = 0.6)) +
  geom_dl(aes(label = label2),  position=position_dodge(dodge_width), method = list(dl.trans(x = x + 0.2), "last.points", cex = 0.6)) +
  xlab("Measurement time") + ylab("Scaled individual and\ncombined biomarker values")

# geom_dl(aes(label = State), method = list(dl.trans(x = x + 0.2), "last.points", cex = 0.8)) +
#   geom_dl(aes(label = State), method = list(dl.trans(x = x - 0.2), "first.points", cex = 0.8)) 

ggsave(p, file = here("figures/Th_component_traj.tiff"), height = 6, width = 7, dpi=300)

#confirm that the scaling is x/sd(x) and send language for the figure blurb


#---------------------------------------------
# ratio_pro_il10 
#---------------------------------------------
plotdf <- df[df$group %in% c("Pro-inflammatory","Pro/IL10 ratio"),] %>% droplevels()
plotdf2 <- df %>% filter(biomarker=="Pro"|biomarker=="IL10") %>%
  mutate(group="Pro/IL10 ratio", combined="0") %>% droplevels()

plotdf <- bind_rows(plotdf, plotdf2)
plotdf$combined <- factor(plotdf$combined)

dodge_width=0.4
p <- ggplot(plotdf, aes(x=time, y=scale, color=tr, group=biomarker_tr, shape=combined, type=combined)) +
  facet_wrap(~group) + 
  geom_point(position=position_dodge(dodge_width)) +
  geom_line(aes(linetype = combined), show.legend=FALSE,  position=position_dodge(dodge_width)) +
  geom_linerange(aes(ymin=ci.lb, ymax=ci.ub), position=position_dodge(dodge_width)) +
  scale_color_manual(name = "Treatment\nArm", values = tableau10) +
  scale_shape_manual(values = c(1, 19), guide=F) +
  scale_linetype_manual(values = c("dashed", "solid")) +
  theme(legend.position = "bottom") +
  geom_dl(aes(label = label1),  position=position_dodge(dodge_width), method = list(dl.trans(x = x - 0.3), "first.points", cex = 0.6)) +
  geom_dl(aes(label = label2),  position=position_dodge(dodge_width), method = list(dl.trans(x = x + 0.2), "last.points", cex = 0.6)) +
  xlab("Measurement time") + ylab("Scaled individual and\ncombined biomarker values")
p

ggsave(p, file = here("figures/pro_il10_component_traj.tiff"), height = 6, width = 7, dpi=300)


#---------------------------------------------
# ratio_th1_il10  
#---------------------------------------------
plotdf <- df[df$group %in% c("Th1","Th1/IL10 ratio"),] %>% droplevels()
plotdf2 <- df %>% filter(biomarker=="Th1"|biomarker=="IL10") %>%
  mutate(group="Th1/IL10 ratio", combined="0") %>% droplevels()

plotdf <- bind_rows(plotdf, plotdf2)
plotdf$combined <- factor(plotdf$combined)

dodge_width=0.4
p <- ggplot(plotdf, aes(x=time, y=scale, color=tr, group=biomarker_tr, shape=combined, type=combined)) +
  facet_wrap(~group) + 
  geom_point(position=position_dodge(dodge_width)) +
  geom_line(aes(linetype = combined), show.legend=FALSE,  position=position_dodge(dodge_width)) +
  geom_linerange(aes(ymin=ci.lb, ymax=ci.ub), position=position_dodge(dodge_width)) +
  scale_color_manual(name = "Treatment\nArm", values = tableau10) +
  scale_shape_manual(values = c(1, 19), guide=F) +
  scale_linetype_manual(values = c("dashed", "solid")) +
  theme(legend.position = "bottom") +
  geom_dl(aes(label = label1),  position=position_dodge(dodge_width), method = list(dl.trans(x = x - 0.3), "first.points", cex = 0.6)) +
  geom_dl(aes(label = label2),  position=position_dodge(dodge_width), method = list(dl.trans(x = x + 0.2), "last.points", cex = 0.6)) +
  xlab("Measurement time") + ylab("Scaled individual and\ncombined biomarker values")
p

ggsave(p, file = here("figures/Th1_il10_component_traj.tiff"), height = 6, width = 7, dpi=300)



#---------------------------------------------
# ratio_th2_il10
#---------------------------------------------

plotdf <- df[df$group %in% c("Th2","Th2/IL10 ratio"),] %>% droplevels()
plotdf2 <- df %>% filter(biomarker=="Th2"|biomarker=="IL10") %>%
  mutate(group="Th2/IL10 ratio", combined="0") %>% droplevels()

plotdf <- bind_rows(plotdf, plotdf2)
plotdf$combined <- factor(plotdf$combined)

dodge_width=0.4
p <- ggplot(plotdf, aes(x=time, y=scale, color=tr, group=biomarker_tr, shape=combined, type=combined)) +
  facet_wrap(~group) + 
  geom_point(position=position_dodge(dodge_width)) +
  geom_line(aes(linetype = combined), show.legend=FALSE,  position=position_dodge(dodge_width)) +
  geom_linerange(aes(ymin=ci.lb, ymax=ci.ub), position=position_dodge(dodge_width)) +
  scale_color_manual(name = "Treatment\nArm", values = tableau10) +
  scale_shape_manual(values = c(1, 19), guide=F) +
  scale_linetype_manual(values = c("dashed", "solid")) +
  theme(legend.position = "bottom") +
  geom_dl(aes(label = label1),  position=position_dodge(dodge_width), method = list(dl.trans(x = x - 0.3), "first.points", cex = 0.6)) +
  geom_dl(aes(label = label2),  position=position_dodge(dodge_width), method = list(dl.trans(x = x + 0.2), "last.points", cex = 0.6)) +
  xlab("Measurement time") + ylab("Scaled individual and\ncombined biomarker values")
p

ggsave(p, file = here("figures/Th2_il10_component_traj.tiff"), height = 6, width = 7, dpi=300)


