
rm(list=ls())
source(here::here("0-config.R"))
library(cowplot)
library(ggpubr)


#Load tmle results
load(here::here("results/immune_adj_glm.RData"))
ls()

readjustfunc <- function(data){
  df <- data.frame(outcome=deparse(substitute(data)), select(as.data.frame(data), RD,ci.lb,ci.ub,`P-value`))
  df$outcome <- gsub("_adj_L","", df$outcome)
  df$outcome <- gsub("t2_","", df$outcome)
  df$outcome <- gsub("t3_","", df$outcome)
  return(df)
}



d <- bind_rows(
  data.frame(readjustfunc(t2_il1_adj_L), name="IL-1b", age="14"),
  data.frame(readjustfunc(t2_il6_adj_L), name="IL-6", age="14"),
  data.frame(readjustfunc(t2_tnf_adj_L), name="TNF-a", age="14"),
  data.frame(readjustfunc(t2_il2_adj_L), name="IL-2", age="14"),
  data.frame(readjustfunc(t2_il12_adj_L), name="IL-12", age="14"),
  data.frame(readjustfunc(t2_ifn_adj_L), name="IFN-g", age="14"),
  data.frame(readjustfunc(t2_il4_adj_L), name="IL-4", age="14"),  
  data.frame(readjustfunc(t2_il5_adj_L), name="IL-5", age="14"),
  data.frame(readjustfunc(t2_il13_adj_L), name="IL-13", age="14"),
  data.frame(readjustfunc(t2_il17_adj_L), name="IL-17", age="14"),
  data.frame(readjustfunc(t2_il21_adj_L), name="IL-21", age="14"),
  data.frame(readjustfunc(t2_il10_adj_L), name="IL-10", age="14"),
  data.frame(readjustfunc(t2_gmc_adj_L), name="GMCSF", age="14"),
  data.frame(readjustfunc(t2_crp_adj_L), name="CRP", age="14"),
  data.frame(readjustfunc(t2_agp_adj_L), name="AGP", age="14"),
  data.frame(readjustfunc(t2_igf_adj_L), name="IGF-1", age="14"),
  
  data.frame(readjustfunc(t3_il1_adj_L), name="IL-1b", age="28"),
  data.frame(readjustfunc(t3_il6_adj_L), name="IL-6", age="28"),
  data.frame(readjustfunc(t3_tnf_adj_L), name="TNF-a", age="28"),
  data.frame(readjustfunc(t3_il2_adj_L), name="IL-2", age="28"),
  data.frame(readjustfunc(t3_il12_adj_L), name="IL-12", age="28"),
  data.frame(readjustfunc(t3_ifn_adj_L), name="IFN-g", age="28"),
  data.frame(readjustfunc(t3_il4_adj_L), name="IL-4", age="28"),  
  data.frame(readjustfunc(t3_il5_adj_L), name="IL-5", age="28"),
  data.frame(readjustfunc(t3_il13_adj_L), name="IL-13", age="28"),
  data.frame(readjustfunc(t3_il17_adj_L), name="IL-17", age="28"),
  data.frame(readjustfunc(t3_il21_adj_L), name="IL-21", age="28"),
  data.frame(readjustfunc(t3_il10_adj_L), name="IL-10", age="28"),
  data.frame(readjustfunc(t3_gmc_adj_L), name="GMCSF", age="28"),
  #data.frame(readjustfunc(t3_crp_adj_L), name="CRP", age="28"),
  #data.frame(readjustfunc(t3_agp_adj_L), name="AGP", age="28"),
  data.frame(readjustfunc(t3_igf_adj_L), name="IGF-1", age="28"),
  
  data.frame(readjustfunc(d23_ln_il1_adj_L), name="IL-1b", age="14-28"),
  data.frame(readjustfunc(d23_ln_il6_adj_L), name="IL-6", age="14-28"),
  data.frame(readjustfunc(d23_ln_tnf_adj_L), name="TNF-a", age="14-28"),
  data.frame(readjustfunc(d23_ln_il2_adj_L), name="IL-2", age="14-28"),
  data.frame(readjustfunc(d23_ln_il12_adj_L), name="IL-12", age="14-28"),
  data.frame(readjustfunc(d23_ln_ifn_adj_L), name="IFN-g", age="14-28"),
  data.frame(readjustfunc(d23_ln_il4_adj_L), name="IL-4", age="14-28"),  
  data.frame(readjustfunc(d23_ln_il5_adj_L), name="IL-5", age="14-28"),
  data.frame(readjustfunc(d23_ln_il13_adj_L), name="IL-13", age="14-28"),
  data.frame(readjustfunc(d23_ln_il17_adj_L), name="IL-17", age="14-28"),
  data.frame(readjustfunc(d23_ln_il21_adj_L), name="IL-21", age="14-28"),
  data.frame(readjustfunc(d23_ln_il10_adj_L), name="IL-10", age="14-28"),
  data.frame(readjustfunc(d23_ln_gmc_adj_L), name="GMCSF", age="14-28"),
  #data.frame(readjustfunc(d23_ln_crp_adj_L), name="CRP", age="14-28"),
  #data.frame(readjustfunc(d23_ln_agp_adj_L), name="AGP", age="14-28"),
  data.frame(readjustfunc(d23_ln_igf_adj_L), name="IGF-1", age="14-28")
)



head(d)
d$Age <- factor(paste0(d$age, " months"), levels = c("14 months", "28 months", "14-28 months"))
d$outcome <- gsub("d23_ln_","",d$outcome)
d$outcome <- factor(d$outcome, levels=unique(d$outcome))
d <- d %>% arrange(Age, outcome)
d$name <- factor(d$name, levels=rev(unique(d$name)))


#' Clustering suggested by Firdaus: 
#' IL-1b, IL-6, TNF-a, IL-2, IL-12, IFN-g, IL-4, 
#' IL-5, IL-13, IL-17, IL-21,IL-10, GMCSF, CRP, AGP, IGF-1.
#' 
#' Here's what I'm interpreting the clustering to be based on Firdaus's suggestions:
#' First cluster, 1 color and listed in this order:
#'  Th1/Th2, IL-12/IL-4, IFN-g/IL-4, IL-12/IL-5, IFN-g/IL-5, IL-12/IL-13, IFN-g/IL-13
#' 
#' I think all the "X cytokines / IL-10" should be 1 color and listed in the 
#' following order:   Pro/IL-10, IL-1b/IL-10, IL-6/IL-10, TNF-a/IL-10, 
#' IL-2/IL-10, Th1/IL-10, Th2/IL-10, IL-12/IL-10, IFN-g/IL-10, IL-4/IL-10,
#'  IL-5/IL-10, IL-13/IL-10, Th17/IL-10, IL-17/IL-10, IL-21/IL-10,
#'   GMCSF/IL-10.
#' 
#' These should be another cluster and color: Th1/IL-17, IL-12/IL-17,
#'  IFN-g/IL-17, I hope that covers all the cytokine ratios, 
#'  let me know if there are any leftover that need a cluster!



d$group<-"one"
dodge=0.6

p <- ggplot(d, aes(x=(name))) + 
  geom_point(aes(shape=group, y=RD, fill=group, color=group, group=group), size = 4, position= position_dodge(width=dodge)) +
  geom_errorbar(aes(ymin=ci.lb, ymax=ci.ub, color=group, group=group), position= position_dodge(width=dodge)) +
  facet_grid(~Age, scales="free_y", labeller = labeller(group = groups), switch = "y") +
  coord_flip(ylim=range(-0.6,0.8)) +
  ylab("Adjusted mean difference (reference: control arm)") +
  xlab("Biomarker") +
  geom_hline(yintercept = 0) +
  scale_shape_manual(values=c(21, 23, 25)) +
  scale_colour_manual(values=tableau11[1]) +
  scale_fill_manual(values=tableau11[1]) +
  ggtitle("Adjusted difference in cytokines\nbetween the WASH+N and control arms") +
  theme(strip.background = element_blank(),
        legend.position="none",
        plot.title = element_text(size = 20, face = "bold"),
        axis.text.y = element_text(size=8, hjust = 1),
        strip.text.x = element_text(size=8, face = "bold"),
        strip.text.y = element_text(size=8, angle = 180, face = "bold"),
        strip.placement = "outside",
        axis.text.x = element_text(size=10, vjust = 0.5),
        panel.spacing.y = unit(0, "lines"),
        legend.box.background = element_rect(colour = "black"), 
        title = element_text(margin=margin(0,0,0,0))) 





ggsave(p, file=here("figures/immune_forest_adj.png"), width = 8, height = 10)

