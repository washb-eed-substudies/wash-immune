rm(list=ls())
source(here::here("0-config.R"))
library(cowplot)
library(ggpubr)


#Load tmle results
load(here::here("results/immune_adj_glm.RData"))

readjustfunc <- function(data){
  df <- data.frame(outcome=deparse(substitute(data)), select(as.data.frame(data), RD,ci.lb,ci.ub,`P-value`))
  df$outcome <- gsub("_adj_L","", df$outcome)
  df$outcome <- gsub("t2_","", df$outcome)
  df$outcome <- gsub("t3_","", df$outcome)
  return(df)
}

d <- rbind(
  data.frame(readjustfunc(t2_ratio_il1_il10_adj_L), name="Interleukin-1β/Interleukin-10", age=14),
  data.frame(readjustfunc(t3_ratio_il1_il10_adj_L), name="Interleukin-1β/Interleukin-10", age=28),
  data.frame(readjustfunc(t2_ratio_il6_il10_adj_L), name="Interleukin-6/Interleukin-10", age=14),
  data.frame(readjustfunc(t3_ratio_il6_il10_adj_L), name="Interleukin-6/Interleukin-10", age=28),
  data.frame(readjustfunc(t2_ratio_tnf_il10_adj_L), name="Tumor necrosis factor-α/Interleukin-10", age=14),
  data.frame(readjustfunc(t3_ratio_tnf_il10_adj_L), name="Tumor necrosis factor-α/Interleukin-10", age=28),
  data.frame(readjustfunc(t2_ratio_il12_il10_adj_L), name="Interleukin-12/Interleukin-10", age=14),
  data.frame(readjustfunc(t3_ratio_il12_il10_adj_L), name="Interleukin-12/Interleukin-10", age=28),
  data.frame(readjustfunc(t2_ratio_ifn_il10_adj_L), name="Interferon-γ/Interleukin-10", age=14),
  data.frame(readjustfunc(t3_ratio_ifn_il10_adj_L), name="Interferon-γ/Interleukin-10", age=28),
  data.frame(readjustfunc(t2_ratio_il4_il10_adj_L), name="Interleukin-4/Interleukin-10", age=14),
  data.frame(readjustfunc(t3_ratio_il4_il10_adj_L), name="Interleukin-4/Interleukin-10", age=28),
  data.frame(readjustfunc(t2_ratio_il5_il10_adj_L), name="Interleukin-5/Interleukin-10", age=14),
  data.frame(readjustfunc(t3_ratio_il5_il10_adj_L), name="Interleukin-5/Interleukin-10", age=28),
  data.frame(readjustfunc(t2_ratio_il13_il10_adj_L), name="Interleukin-13/Interleukin-10", age=14),
  data.frame(readjustfunc(t3_ratio_il13_il10_adj_L), name="Interleukin-13/Interleukin-10", age=28),
  data.frame(readjustfunc(t2_ratio_il17_il10_adj_L), name="Interleukin-17/Interleukin-10", age=14),
  data.frame(readjustfunc(t3_ratio_il17_il10_adj_L), name="Interleukin-17/Interleukin-10", age=28),
  data.frame(readjustfunc(t2_ratio_il21_il10_adj_L), name="Interleukin-21/Interleukin-10", age=14),
  data.frame(readjustfunc(t3_ratio_il21_il10_adj_L), name="Interleukin-21/Interleukin-10", age=28),
  data.frame(readjustfunc(t2_ratio_il2_il10_adj_L), name="Interleukin-2/Interleukin-10", age=14),
  data.frame(readjustfunc(t3_ratio_il2_il10_adj_L), name="Interleukin-2/Interleukin-10", age=28),
  data.frame(readjustfunc(t2_ratio_gmc_il10_adj_L), name="Granulocyte-macrophage \n colony-stimulating factor/Interleukin-10", age=14),
  data.frame(readjustfunc(t3_ratio_gmc_il10_adj_L), name="Granulocyte-macrophage \n colony-stimulating factor/Interleukin-10", age=28),
  data.frame(readjustfunc(t2_ratio_il12_il4_adj_L), name="Interleukin-12/Interleukin-4", age=14),
  data.frame(readjustfunc(t3_ratio_il12_il4_adj_L), name="Interleukin-12/Interleukin-4", age=28),
  data.frame(readjustfunc(t2_ratio_ifn_il4_adj_L), name="Interferon-γ/Interleukin-4", age=14),
  data.frame(readjustfunc(t3_ratio_ifn_il4_adj_L), name="Interferon-γ/Interleukin-4", age=28),
  data.frame(readjustfunc(t2_ratio_il12_il5_adj_L), name="Interleukin-12/Interleukin-5", age=14),
  data.frame(readjustfunc(t3_ratio_il12_il5_adj_L), name="Interleukin-12/Interleukin-5", age=28),
  data.frame(readjustfunc(t2_ratio_ifn_il5_adj_L), name="Interferon-γ/Interleukin-5", age=14),
  data.frame(readjustfunc(t3_ratio_ifn_il5_adj_L), name="Interferon-γ/Interleukin-5", age=28),
  data.frame(readjustfunc(t2_ratio_il12_il13_adj_L), name="Interleukin-12/Interleukin-13", age=14),
  data.frame(readjustfunc(t3_ratio_il12_il13_adj_L), name="Interleukin-12/Interleukin-13", age=28),
  data.frame(readjustfunc(t2_ratio_ifn_il13_adj_L), name="Interferon-γ/Interleukin-13", age=14),
  data.frame(readjustfunc(t3_ratio_ifn_il13_adj_L), name="Interferon-γ/Interleukin-13", age=28),
  data.frame(readjustfunc(t2_ratio_il12_il17_adj_L), name="Interleukin-12/Interleukin-17", age=14),
  data.frame(readjustfunc(t3_ratio_il12_il17_adj_L), name="Interleukin-12/Interleukin-17", age=28),
  data.frame(readjustfunc(t2_ratio_ifn_il17_adj_L), name="Interferon-γ/Interleukin-17", age=14),
  data.frame(readjustfunc(t3_ratio_ifn_il17_adj_L), name="Interferon-γ/Interleukin-17", age=28),
  data.frame(readjustfunc(t2_ratio_il12_il21_adj_L), name="Interleukin-12/Interleukin-21", age=14),
  data.frame(readjustfunc(t3_ratio_il12_il21_adj_L), name="Interleukin-12/Interleukin-21", age=28),
  data.frame(readjustfunc(t2_ratio_ifn_il21_adj_L), name="Interferon-γ/Interleukin-21", age=14),
  data.frame(readjustfunc(t3_ratio_ifn_il21_adj_L), name="Interferon-γ/Interleukin-21", age=28),
  data.frame(readjustfunc(t2_ratio_pro_il10_adj_L), name="Pro-inflammatory cytokines/Interleukin-10", age=14),
  data.frame(readjustfunc(t3_ratio_pro_il10_adj_L), name="Pro-inflammatory cytokines/Interleukin-10", age=28),
  data.frame(readjustfunc(t2_ratio_th1_il10_adj_L), name="Th1/Interleukin-10", age=14),
  data.frame(readjustfunc(t3_ratio_th1_il10_adj_L), name="Th1/Interleukin-10", age=28),
  data.frame(readjustfunc(t2_ratio_th2_il10_adj_L), name="Th2/Interleukin-10", age=14),
  data.frame(readjustfunc(t3_ratio_th2_il10_adj_L), name="Th2/Interleukin-10", age=28),
  data.frame(readjustfunc(t2_ratio_th17_il10_adj_L), name="Th17/Interleukin-10", age=14),
  data.frame(readjustfunc(t3_ratio_th17_il10_adj_L), name="Th17/Interleukin-10", age=28),
  data.frame(readjustfunc(t2_ratio_th1_th2_adj_L), name="Th1/Th2", age=14),
  data.frame(readjustfunc(t3_ratio_th1_th2_adj_L), name="Th1/Th2", age=28),
  data.frame(readjustfunc(t2_ratio_th1_th17_adj_L), name="Th1/Th17", age=14),
  data.frame(readjustfunc(t3_ratio_th1_th17_adj_L), name="Th1/Th17", age=28)
)


head(d)
d$Age <- factor(paste0(d$age, " months"))

dodge=0.6
ggplot(d, aes(x=(name))) + 
  geom_point(aes(shape=Age, y=RD, fill=Age, color=Age, group=Age), size = 4, position= position_dodge(width=dodge)) +
  geom_errorbar(aes(ymin=ci.lb, ymax=ci.ub, color=Age, group=Age), position= position_dodge(width=dodge)) +
  coord_flip(ylim=range(-0.5,0.5)) +
  #labs(x = "Study-specific results stratified by risk factor level\nwith reference category N's and cases printed", y = Ylab) +
  ylab("Adjusted mean difference (reference: control arm)") +
  xlab("Biomarker") +
  geom_hline(yintercept = 0) +
  #scale_x_discrete(labels= df$studyid2) +
  scale_shape_manual(values=c(21, 23)) +
  scale_colour_manual(values=tableau10[c(1:5)]) +
  scale_fill_manual(values=tableau10[c(1:5)]) +
  theme(strip.background = element_blank(),
        legend.position="right",
        strip.text.x = element_text(size=12),
        axis.text.x = element_text(size=12),
        axis.title.x = element_text(size=12)) 



#Facet by age and color/arrange by group:
head(d)
unique(d$outcome)
d$Category<- "Other"
d$Category[grepl("ifn",d$outcome)] <-"Interferon"
d$Category[grepl("_il13",d$outcome)] <-"IL 13"
d$Category[grepl("_il10",d$outcome)] <-"IL 10"
d$Category[grepl("th1_",d$outcome)] <- "TH1"
d$Category <- factor(d$Category, levels=rev(c("TH1", "IL 10", "IL 13", "Interferon", "Other")))
d <- d %>% arrange(Category, age, RD) %>% mutate(name=factor(name, levels=unique(name)))


ggplot(d, aes(x=(name))) + 
  geom_point(aes(y=RD, fill=Category, color=Category, group=Category), size = 4) +
  geom_errorbar(aes(ymin=ci.lb, ymax=ci.ub, color=Category, group=Category)) +
  coord_flip(ylim=range(-0.5,0.5)) +
  facet_wrap(~Age) +
  #labs(x = "Study-specific results stratified by risk factor level\nwith reference category N's and cases printed", y = Ylab) +
  ylab("Adjusted mean difference (reference: control arm)") +
  xlab("Biomarker") +
  geom_hline(yintercept = 0) +
  #scale_x_discrete(labels= df$studyid2) +
  #scale_shape_manual(values=c(21, 23)) +
  scale_colour_manual(values=tableau10[c(5:1)]) +
  scale_fill_manual(values=tableau10[c(5:1)]) +
  theme(strip.background = element_blank(),
        legend.position="right",
        strip.text.x = element_text(size=12),
        axis.text.x = element_text(size=12),
        axis.title.x = element_text(size=12)) 

