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
unique(d$outcome)

d <- d %>% mutate(
  group=case_when(
    outcome %in% c("ratio_th1_th2", "ratio_il12_il4", "ratio_ifn_il4", "ratio_il12_il5", "ratio_ifn_il5",  "ratio_il12_il13", "ratio_ifn_il13") ~"one", 
    outcome %in% c("ratio_pro_il10", "ratio_il1_il10","ratio_il6_il10", "ratio_tnf_il10", 
                   "ratio_il2_il10",  "ratio_th1_il10",  "ratio_th2_il10", "ratio_il12_il10", "ratio_ifn_il10", "ratio_il4_il10", 
                   "ratio_il5_il10", "ratio_il13_il10","ratio_th17_il10", "ratio_il17_il10", "ratio_il21_il10",
                   "ratio_gmc_il10") ~ "two",
    outcome %in% c("ratio_th1_th17", "ratio_il12_il17", "ratio_ifn_il17", "ratio_il12_il21", "ratio_ifn_il21") ~"three"
    ),
  group=factor(group, level=c("one","two","three")),
  outcome = factor(outcome, levels =c(
    "ratio_th1_th2", "ratio_il12_il4", "ratio_ifn_il4", "ratio_il12_il5", "ratio_ifn_il5",  "ratio_il12_il13", "ratio_ifn_il13",
    "ratio_pro_il10", "ratio_il1_il10","ratio_il6_il10", "ratio_tnf_il10", 
    "ratio_il2_il10",  "ratio_th1_il10",  "ratio_th2_il10", "ratio_il12_il10", "ratio_ifn_il10", "ratio_il4_il10", 
    "ratio_il5_il10", "ratio_il13_il10","ratio_th17_il10", "ratio_il17_il10", "ratio_il21_il10",
    "ratio_gmc_il10",
    "ratio_th1_th17", "ratio_il12_il17", "ratio_ifn_il17", "ratio_il12_il21", "ratio_ifn_il21"))
)


groups <- c(
  "one" = "",
  "two" = "",
  "three" = "")

d$Age <- NA
d$Age[d$age==14] <- "14 months"
d$Age[d$age==28] <- "28 months"




dodge=0.6
ggplot(d, aes(x=(name))) + 
  geom_point(aes(shape=group, y=RD, fill=group, color=group, group=group), size = 4, position= position_dodge(width=dodge)) +
  geom_errorbar(aes(ymin=ci.lb, ymax=ci.ub, color=group, group=group), position= position_dodge(width=dodge)) +
  facet_grid(group~Age, scales="free_y", labeller = labeller(group = groups), switch = "y") +
  coord_flip(ylim=range(-0.5,0.5)) +
  #labs(x = "Study-specific results stratified by risk factor level\nwith reference category N's and cases printed", y = Ylab) +
  ylab("Adjusted mean difference (reference: control arm)") +
  xlab("Biomarker") +
  geom_hline(yintercept = 0) +
  #scale_x_discrete(labels= df$studyid2) +
  scale_shape_manual(values=c(21, 23, 25)) +
  scale_colour_manual(values=tableau10[c(1:5)]) +
  scale_fill_manual(values=tableau10[c(1:5)]) +
  ggtitle("Adjusted difference in Cytokine ratios\nbetween the WASH+N and the control arm") +
  theme(strip.background = element_blank(),
        legend.position="right",
        plot.title = element_text(size = 20, face = "bold"),
        # strip.text.x = element_text(size=12),
        # axis.text.x = element_text(size=12),
        # axis.title.x = element_text(size=12)) ,
        axis.text.y = element_text(size=8, hjust = 1),
        strip.text.x = element_text(size=8, face = "bold"),
        strip.text.y = element_text(size=8, angle = 180, face = "bold"),
        strip.placement = "outside",
        axis.text.x = element_text(size=10, vjust = 0.5),
        panel.spacing = unit(0, "lines"),
        legend.box.background = element_rect(colour = "black"), 
        title = element_text(margin=margin(0,0,0,0))) 




