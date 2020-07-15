
rm(list=ls())
source(here::here("0-config.R"))
library(cowplot)
library(ggpubr)


#Load tmle results
load(here::here("results/immune_tr_means.RData"))
ls()

readjustfunc <- function(data){
  df <- data.frame(outcome=deparse(substitute(data)), as.data.frame(data))
  colnames(df)[3] <-"N"
  df$outcome <- gsub("_N_tr","", df$outcome)
  df$outcome <- gsub("","", df$outcome)
  df$outcome <- gsub("t3_","", df$outcome)
  return(df)
}



d <- bind_rows(
  data.frame(readjustfunc(il1_t2_N_tr), name="IL-1b", age="14"),
  data.frame(readjustfunc(il6_t2_N_tr), name="IL-6", age="14"),
  data.frame(readjustfunc(tnf_t2_N_tr), name="TNF-a", age="14"),
  data.frame(readjustfunc(il2_t2_N_tr), name="IL-2", age="14"),
  data.frame(readjustfunc(il12_t2_N_tr), name="IL-12", age="14"),
  data.frame(readjustfunc(ifn_t2_N_tr), name="IFN-g", age="14"),
  data.frame(readjustfunc(il4_t2_N_tr), name="IL-4", age="14"),  
  data.frame(readjustfunc(il5_t2_N_tr), name="IL-5", age="14"),
  data.frame(readjustfunc(il13_t2_N_tr), name="IL-13", age="14"),
  data.frame(readjustfunc(il17_t2_N_tr), name="IL-17", age="14"),
  data.frame(readjustfunc(il21_t2_N_tr), name="IL-21", age="14"),
  data.frame(readjustfunc(il10_t2_N_tr), name="IL-10", age="14"),
  data.frame(readjustfunc(gmc_t2_N_tr), name="GMCSF", age="14"),
  data.frame(readjustfunc(crp_t2_N_tr), name="CRP", age="14"),
  data.frame(readjustfunc(agp_t2_N_tr), name="AGP", age="14"),
  data.frame(readjustfunc(igf_t2_N_tr), name="IGF-1", age="14"),
  
  data.frame(readjustfunc(il1_t3_N_tr), name="IL-1b", age="28"),
  data.frame(readjustfunc(il6_t3_N_tr), name="IL-6", age="28"),
  data.frame(readjustfunc(tnf_t3_N_tr), name="TNF-a", age="28"),
  data.frame(readjustfunc(il2_t3_N_tr), name="IL-2", age="28"),
  data.frame(readjustfunc(il12_t3_N_tr), name="IL-12", age="28"),
  data.frame(readjustfunc(ifn_t3_N_tr), name="IFN-g", age="28"),
  data.frame(readjustfunc(il4_t3_N_tr), name="IL-4", age="28"),  
  data.frame(readjustfunc(il5_t3_N_tr), name="IL-5", age="28"),
  data.frame(readjustfunc(il13_t3_N_tr), name="IL-13", age="28"),
  data.frame(readjustfunc(il17_t3_N_tr), name="IL-17", age="28"),
  data.frame(readjustfunc(il21_t3_N_tr), name="IL-21", age="28"),
  data.frame(readjustfunc(il10_t3_N_tr), name="IL-10", age="28"),
  data.frame(readjustfunc(gmc_t3_N_tr), name="GMCSF", age="28"),
  #data.frame(readjustfunc(crp_t3_N_tr), name="CRP", age="28"),
  #data.frame(readjustfunc(agp_t3_N_tr), name="AGP", age="28"),
  data.frame(readjustfunc(igf_t3_N_tr), name="IGF-1", age="28")
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
d$Arm <- d$tr
p <- ggplot(d, aes(x=(name))) + 
  geom_point(aes(y=Mean, fill=Arm, color=Arm, group=Arm), size = 4, position= position_dodge(width=dodge)) +
  geom_errorbar(aes(ymin=Lower.95.CI, ymax=Upper.95.CI, color=Arm, group=Arm), position= position_dodge(width=dodge)) +
  facet_grid(~Age, scales="free_y", switch = "y") +
  coord_flip(ylim=range(-0.4,7.2)) +
  ylab("Treatment-arm specific log-transformed concentrations") +
  xlab("Biomarker") +
  geom_hline(yintercept = 0) +
  scale_colour_manual(values=tableau10[2:1]) +
  scale_fill_manual(values=tableau10[2:1]) +
  ggtitle("Mean concentrations of cytokines\nin the WASH+N and control arms") +
  theme(strip.background = element_blank(),
        legend.position="bottom",
        plot.title = element_text(size = 20, face = "bold"),
        axis.text.y = element_text(size=8, hjust = 1),
        strip.text.x = element_text(size=8, face = "bold"),
        strip.text.y = element_text(size=8, angle = 180, face = "bold"),
        strip.placement = "outside",
        axis.text.x = element_text(size=10, vjust = 0.5),
        panel.spacing = unit(0, "lines"),
        legend.box.background = element_rect(colour = "black"), 
        title = element_text(margin=margin(0,0,0,0))) 





ggsave(p, file=here("figures/immune_forest_means.png"), width = 8, height = 10)

