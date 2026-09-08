rm(list=ls())
library("xtable")
source(here::here("0-config.R"))

source(here::here("table scripts/immune-main-tables.R"))
load(here::here("results/immune_ipcw.RData"))

setwd(paste0(dropboxDir,"Data/Cleaned/Audrie/"))
lab<-read.csv("bangladesh-dm-ee-anthro-diar-ee-med-plasma-blind-tr-enrol-covariates-lab.csv", stringsAsFactors = TRUE)

## Tables S3/S5/S7 borrow the unadjusted and fully adjusted p-value
## columns from main Tables 2-4, which show only the corrected value.
## Relabel the copy so the supplement can put a raw column next to it
## without disturbing the main table.
as_bonf_col <- function(v){
  v[1] <- "Bonferroni P-value"
  v
}

## IMPORTANT FUNCTIONS ##
#to be used for formatting ipcw variables for table
maketblvalue<-function(ipcwadjvar){
  rounded<-round(ipcwadjvar, 2)
  paste(rounded[1], " (", rounded[3], ", ", rounded[4], ")", sep="")
}


#### TABLE S1 ####
# filtering for children with no t3 measurements
lost<-y1 %>% filter_at(vars(igf_t3, gmcsf_t3, ifng_t3, il10_t3, il12_t3, il13_t3, il17_t3,
                              il1_t3, il2_t3, il21_t3, il4_t3, il5_t3, il6_t3, tnfa_t3), all_vars(is.na(.)))

#calculating overall N by arm
Nlostctrl<-nrow(lost[lost$tr=="Control",])
Nlostwsh<-nrow(lost[lost$tr=="Nutrition + WSH",])

imomage<-meansdfunc(lost, lost$momage)
imomeduy<-meansdfunc(lost, lost$momeduy)
idadeduy<-meansdfunc(lost, lost$dadeduy)
idadagri<-npercfunc(lost, lost$dadagri)
iNhh<-meansdfunc(lost, lost$Nhh)
ielec<-npercfunc(lost, lost$elec)
icement<-npercfunc(lost, lost$cement)

iacresmctrl<-round(mean(lost$landacre[lost$tr=="Control"], na.rm=TRUE), 2)
iacressdctrl<-round(sd(lost$landacre[lost$tr=="Control"], na.rm=TRUE), 2)
iacresmwsh<-round(mean(lost$landacre[lost$tr=="Nutrition + WSH"], na.rm=TRUE), 2)
iacressdwsh<-round(mean(lost$landacre[lost$tr=="Nutrition + WSH"], na.rm=TRUE), 2)
iacres<-c(iacresmctrl, iacressdctrl, iacresmwsh, iacressdwsh)

itubewell<-npercfunc(lost, lost$tubewell)
istorewater<-npercfunc(lost, lost$storewat)
itreatwater<-npercfunc(lost, lost$treatwat)
iwaterdis<-meansdfunc(lost, lost$watmin)
iodmen<-npercfunc(lost, lost$odmen)
iodwomen<-npercfunc(lost, lost$odwom)
iodchild815<-npercfunc(lost, lost$odch815)
iodchild38<-npercfunc(lost, lost$odch38)
iodchild03<-npercfunc(lost, lost$odchu3)
ilatowned<-npercfunc(lost, lost$latown)
ilatslab<-npercfunc(lost, lost$latslab)
ilatseal<-npercfunc(lost, lost$latseal)

ilatfctrln<-sum(lost$latfeces[lost$tr=="Control"]==0, na.rm=T)
ilatfctrlperc<-round(ilatfctrln/sum(!is.na(lost$latfeces[lost$tr=="Control"]), na.rm=T)*100)
ilatfwshn<-sum(lost$latfeces[lost$tr=="Nutrition + WSH"]==0, na.rm=T)
ilatfwshperc<-round(ilatfwshn/sum(!is.na(lost$latfeces[lost$tr=="Nutrition + WSH"]), na.rm=T)*100)
ilatfeces<-c(ilatfctrln, ilatfctrlperc, ilatfwshn, ilatfwshperc)

ipotty<-npercfunc(lost, lost$potty)
ifeceshouse<-npercfunc(lost, lost$humfeces)
ifeceschildarea<-npercfunc(lost, lost$humfecesch)
ihandlatwater<-npercfunc(lost, lost$hwlatwat)
ihandlatsoap<-npercfunc(lost, lost$hwlatsoap)
ihandkitwater<-npercfunc(lost, lost$hwkitwat)
ihandkitsoap<-npercfunc(lost, lost$hwkitsoap)

ifsnctrl<-length(lost$hfiacat[lost$tr=="Control" & lost$hfiacat=="Food Secure"])
ifspercctrl<-round(ifsnctrl/length(lost$hfiacat[lost$tr=="Control"])*100)
ifsnwsh<-length(lost$hfiacat[lost$tr=="Nutrition + WSH" & lost$hfiacat=="Food Secure"])
ifspercwsh<-round(ifsnwsh/length(lost$hfiacat[lost$tr=="Nutrition + WSH"])*100)
ifoodsecure<-c(ifsnctrl, ifspercctrl, ifsnwsh, ifspercwsh)

#make vectors to put in table
#function combines n and percent or mean and sd for vectors created from npercfunc or meansdfunc
#num is 1 if ctrl group, 3 if wsh
ctrl<-c(paste("Control (N=", Nlostctrl, ")", sep=""), " ", charobject(imomage, 1),charobject(imomeduy, 1), " ", charobject(idadeduy, 1), charobjectperc(idadagri, 1),
               " ", charobject(iNhh, 1), charobjectperc(ielec, 1), charobjectperc(icement, 1), charobject(iacres, 1),
               " ", charobjectperc(itubewell, 1), charobjectperc(istorewater, 1), charobjectperc(itreatwater, 1), charobject(iwaterdis, 1), 
               " ", " ", charobjectperc(iodmen, 1), charobjectperc(iodwomen, 1), charobjectperc(iodchild815, 1), charobjectperc(iodchild38, 1), charobjectperc(iodchild03, 1), 
               " ", charobjectperc(ilatowned, 1), charobjectperc(ilatslab, 1), charobjectperc(ilatseal, 1), charobjectperc(ilatfeces, 1),
               charobjectperc(ipotty, 1), 
               " ", charobjectperc(ifeceshouse, 1), charobjectperc(ifeceschildarea, 1), 
               " ", " ", charobjectperc(ihandlatwater, 1), charobjectperc(ihandlatsoap, 1), 
               " ", charobjectperc(ihandkitwater, 1), charobjectperc(ihandkitsoap, 1), 
               " ", charobjectperc(ifoodsecure, 1))
wsh<-c(paste("N + WSH Intervention (N=", Nlostwsh, ")", sep=""), " ", charobject(imomage, 3),charobject(imomeduy, 3), " ", charobject(idadeduy, 3), charobjectperc(idadagri, 3),
       " ", charobject(iNhh, 3), charobjectperc(ielec, 3), charobjectperc(icement, 3), charobject(iacres, 3),
       " ", charobjectperc(itubewell, 3), charobjectperc(istorewater, 3), charobjectperc(itreatwater, 3), charobject(iwaterdis, 3), 
       " ", " ", charobjectperc(iodmen, 3), charobjectperc(iodwomen, 3), charobjectperc(iodchild815, 3), charobjectperc(iodchild38, 3), charobjectperc(iodchild03, 3), 
       " ", charobjectperc(ilatowned, 3), charobjectperc(ilatslab, 3), charobjectperc(ilatseal, 3), charobjectperc(ilatfeces, 3),
       charobjectperc(ipotty, 3), 
       " ", charobjectperc(ifeceshouse, 3), charobjectperc(ifeceschildarea, 3), 
       " ", " ", charobjectperc(ihandlatwater, 3), charobjectperc(ihandlatsoap, 3), 
       " ", charobjectperc(ihandkitwater, 3), charobjectperc(ihandkitsoap, 3), 
       " ", charobjectperc(ifoodsecure, 3))

# Table S1/S2
tbls1<-data.table(" "=c("No. of compounds:", "Maternal", "Age(years)", "Years of education", 
                          "Paternal", "Years of education", "Works in agriculture", 
                          "Household", "Number of people", "Has electricity", "Has a cement floor", "Acres of agricultural land owned", 
                          "Drinking Water", "Shallow tubewell primary water source", "Stored water observed at home", "Reported treating water yesterday", "Distance (mins) to primary water source",
                          "Sanitation", "Reported daily open defecation", "Adult men", "Adult women", "Children: 8 to <15 years", "Children: 3 to <8 years", "Children: 0 to <3 years", 
                          "Latrine", "Owned", "Concrete Slab", "Functional water seal", "Visible stool on slab or floor",
                          "Owned a child potty",
                          "Human feces observed in the", "House", "Child's play area",
                          "Handwashing location", "Within six steps of latrine", "Has water", "Has soap", "Within six steps of kitchen", "Has water", "Has soap", 
                          "Nutrition", "Household is food secure"), 
                  "WASH Benefits Main Trial"=c("Control (N=1382)"," ", "24 (5)", "6 (3)", " ", "5 (4)", "414 (30%)", 
                                               " ", "5 (2)", "784 (57%)", "145 (10%)", "0.15 (0.21)",
                                               " ", "1038 (75%)", "666 (48%)", "4 (0%)", "1 (3)", 
                                               " ", " ", "97 (7%)", "62 (4%)", "53 (10%)", "267 (38%)", "245 (82%)",
                                               " ", "750 (54%)", "1251 (95%)", "358 (31%)", "625 (48%)", "61 (4%)", 
                                               " ", "114 (8%)", "21 (2%)", 
                                               " ", " ", "178 (14%)", "88 (7%)", " ", "118 (9%)", "33 (3%)", " ", "932 (67%)"),
                  " "=c ("N + WSH Intervention (N=686)", " ", "24 (6)", "6 (3)", " ", "5 (4)", "207 (30%)", 
                         " ", "5 (2)", "412 (60%)", "72 (10%)", "0.14 (0.38)",
                         " ", "504 (73%)", "331 (48%)", "2 (0%)", "1 (2)", 
                         " ", " ", "50 (7%)", "24 (4%)", "28 (10%)", "134 (37%)", "123 (88%)",
                         " ", "367 (53%)", "621 (94%)", "155 (27%)", "298 (46%)", "30 (4%)", 
                         " ", "49 (8%)", "7 (1%)", 
                         " ", " ", "72 (11%)", "36 (6%)", " ", "60 (9%)", "18 (3%)", " ", "485 (71%)"), 
                  "Immune Status Study: Had outcomes at Year 1"=ctrly1,
                  " "=wshy1, 
                  "Immune Status Study: Lost to follow-up at Year 2"=ctrl,
                  " "=wsh
)

write.csv(tbls1, file=here('tables/supplementary/immune_supptable1.csv'))
print(xtable(tbls1), type="html", file=here("tables/supplementary/immune_supptable1.html"))




#### TABLE S2 ####

outcomes2 <- c("Outcome, Arm", paste("Ln IL-1", "\u03b2", " (pg/ml)", sep=""), "Control", "Nutrition + WSH", 
             "Ln IL-6 (pg/ml)", "Control", "Nutrition + WSH", 
             paste("Ln TNF-", "\u03b1", " (pg/ml)", sep=""), "Control", "Nutrition + WSH",
             "Ln CRP (mg/L)", "Control", "Nutrition + WSH", 
             "Ln IL-12 (pg/ml)", "Control", "Nutrition + WSH", 
             paste("Ln IFN-", "\u03b3", " (pg/ml)", sep=""), "Control", "Nutrition + WSH", 
             "Ln IL-4 (pg/ml)", "Control", "Nutrition + WSH", 
             "Ln IL-5 (pg/ml)", "Control", "Nutrition + WSH", 
             "Ln IL-13 (pg/ml)", "Control", "Nutrition + WSH", 
             "Ln IL-17A (pg/ml)", "Control", "Nutrition + WSH", 
             "Ln IL-21 (pg/ml)", "Control", "Nutrition + WSH", 
             "Ln IL-10 (pg/ml)", "Control", "Nutrition + WSH", 
             "Ln IL-2 (pg/ml)", "Control", "Nutrition + WSH", 
             "Ln GM-CSF (pg/ml)", "Control", "Nutrition + WSH", 
             "Ln AGP (g/L)", "Control", "Nutrition + WSH",
             paste("Ln IGF-1 (", "\u03bc", "g/L)", sep=""), "Control", "Nutrition + WSH")

Ns2 <- c("N", " ", as.character(il1_t2_N_tr$t2_ln_il1_N_tr[1]), as.character(il1_t2_N_tr$t2_ln_il1_N_tr[2]),
       " ", as.character(il6_t2_N_tr$t2_ln_il6_N_tr[1]), as.character(il6_t2_N_tr$t2_ln_il6_N_tr[2]),
       " ", as.character(tnf_t2_N_tr$t2_ln_tnf_N_tr[1]), as.character(tnf_t2_N_tr$t2_ln_tnf_N_tr[2]),
       " ", as.character(crp_t2_N_tr$t2_ln_crp_N_tr[1]), as.character(crp_t2_N_tr$t2_ln_crp_N_tr[2]),
       " ", as.character(il12_t2_N_tr$t2_ln_il12_N_tr[1]), as.character(il12_t2_N_tr$t2_ln_il12_N_tr[2]),
       " ", as.character(ifn_t2_N_tr$t2_ln_ifn_N_tr[1]), as.character(ifn_t2_N_tr$t2_ln_ifn_N_tr[2]),
       " ", as.character(il4_t2_N_tr$t2_ln_il4_N_tr[1]), as.character(il4_t2_N_tr$t2_ln_il4_N_tr[2]),
       " ", as.character(il5_t2_N_tr$t2_ln_il5_N_tr[1]), as.character(il5_t2_N_tr$t2_ln_il5_N_tr[2]),
       " ", as.character(il13_t2_N_tr$t2_ln_il13_N_tr[1]), as.character(il13_t2_N_tr$t2_ln_il13_N_tr[2]),
       " ", as.character(il17_t2_N_tr$t2_ln_il17_N_tr[1]), as.character(il17_t2_N_tr$t2_ln_il17_N_tr[2]),
       " ", as.character(il21_t2_N_tr$t2_ln_il21_N_tr[1]), as.character(il21_t2_N_tr$t2_ln_il21_N_tr[2]),
       " ", as.character(il10_t2_N_tr$t2_ln_il10_N_tr[1]), as.character(il10_t2_N_tr$t2_ln_il10_N_tr[2]),
       " ", as.character(il2_t2_N_tr$t2_ln_il2_N_tr[1]), as.character(il2_t2_N_tr$t2_ln_il2_N_tr[2]),
       " ", as.character(gmc_t2_N_tr$t2_ln_gmc_N_tr[1]), as.character(gmc_t2_N_tr$t2_ln_gmc_N_tr[2]),
       " ", as.character(agp_t2_N_tr$t2_ln_agp_N_tr[1]), as.character(agp_t2_N_tr$t2_ln_agp_N_tr[2]),
       " ", as.character(igf_t2_N_tr$t2_ln_igf_N_tr[1]), as.character(igf_t2_N_tr$t2_ln_igf_N_tr[2]))

abssds2 <-c("Absolute SD", " ", as.character(round(abs_il1_t2_N_tr$sd[1], 2)), as.character(round(abs_il1_t2_N_tr$sd[2], 2)),
          " ", as.character(round(abs_il6_t2_N_tr$sd[1], 2)), as.character(round(abs_il6_t2_N_tr$sd[2], 2)),
          " ", as.character(round(abs_tnf_t2_N_tr$sd[1], 2)), as.character(round(abs_tnf_t2_N_tr$sd[2], 2)),
          " ", as.character(round(abs_crp_t2_N_tr$sd[1], 2)), as.character(round(abs_crp_t2_N_tr$sd[2], 2)),
          " ", as.character(round(abs_il12_t2_N_tr$sd[1], 2)), as.character(round(abs_il12_t2_N_tr$sd[2], 2)),
          " ", as.character(round(abs_ifn_t2_N_tr$sd[1], 2)), as.character(round(abs_ifn_t2_N_tr$sd[2], 2)),
          " ", as.character(round(abs_il4_t2_N_tr$sd[1], 2)), as.character(round(abs_il4_t2_N_tr$sd[2], 2)),
          " ", as.character(round(abs_il5_t2_N_tr$sd[1], 2)), as.character(round(abs_il5_t2_N_tr$sd[2], 2)),
          " ", as.character(round(abs_il13_t2_N_tr$sd[1], 2)), as.character(round(abs_il13_t2_N_tr$sd[2], 2)),
          " ", as.character(round(abs_il17_t2_N_tr$sd[1], 2)), as.character(round(abs_il17_t2_N_tr$sd[2], 2)),
          " ", as.character(round(abs_il21_t2_N_tr$sd[1], 2)), as.character(round(abs_il21_t2_N_tr$sd[2], 2)),
          " ", as.character(round(abs_il10_t2_N_tr$sd[1], 2)), as.character(round(abs_il10_t2_N_tr$sd[2], 2)),
          " ", as.character(round(abs_il2_t2_N_tr$sd[1], 2)), as.character(round(abs_il2_t2_N_tr$sd[2], 2)),
          " ", as.character(round(abs_gmc_t2_N_tr$sd[1], 2)), as.character(round(abs_gmc_t2_N_tr$sd[2], 2)),
          " ", as.character(round(abs_agp_t2_N_tr$sd[1], 2)), as.character(round(abs_agp_t2_N_tr$sd[2], 2)),
          " ", as.character(round(abs_igf_t2_N_tr$sd[1], 2)), as.character(round(abs_igf_t2_N_tr$sd[2], 2)))

absmeans2 <- c("Absolute Mean", " ", as.character(round(abs_il1_t2_N_tr$mean[1], 2)), as.character(round(abs_il1_t2_N_tr$mean[2], 2)),
             " ", as.character(round(abs_il6_t2_N_tr$mean[1], 2)), as.character(round(abs_il6_t2_N_tr$mean[2], 2)),
             " ", as.character(round(abs_tnf_t2_N_tr$mean[1], 2)), as.character(round(abs_tnf_t2_N_tr$mean[2], 2)),
             " ", as.character(round(abs_crp_t2_N_tr$mean[1], 2)), as.character(round(abs_crp_t2_N_tr$mean[2], 2)),
             " ", as.character(round(abs_il12_t2_N_tr$mean[1], 2)), as.character(round(abs_il12_t2_N_tr$mean[2], 2)),
             " ", as.character(round(abs_ifn_t2_N_tr$mean[1], 2)), as.character(round(abs_ifn_t2_N_tr$mean[2], 2)),
             " ", as.character(round(abs_il4_t2_N_tr$mean[1], 2)), as.character(round(abs_il4_t2_N_tr$mean[2], 2)),
             " ", as.character(round(abs_il5_t2_N_tr$mean[1], 2)), as.character(round(abs_il5_t2_N_tr$mean[2], 2)),
             " ", as.character(round(abs_il13_t2_N_tr$mean[1], 2)), as.character(round(abs_il13_t2_N_tr$mean[2], 2)),
             " ", as.character(round(abs_il17_t2_N_tr$mean[1], 2)), as.character(round(abs_il17_t2_N_tr$mean[2], 2)),
             " ", as.character(round(abs_il21_t2_N_tr$mean[1], 2)), as.character(round(abs_il21_t2_N_tr$mean[2], 2)),
             " ", as.character(round(abs_il10_t2_N_tr$mean[1], 2)), as.character(round(abs_il10_t2_N_tr$mean[2], 2)),
             " ", as.character(round(abs_il2_t2_N_tr$mean[1], 2)), as.character(round(abs_il2_t2_N_tr$mean[2], 2)),
             " ", as.character(round(abs_gmc_t2_N_tr$mean[1], 2)), as.character(round(abs_gmc_t2_N_tr$mean[2], 2)),
             " ", as.character(round(abs_agp_t2_N_tr$mean[1], 2)), as.character(round(abs_agp_t2_N_tr$mean[2], 2)),
             " ", as.character(round(abs_igf_t2_N_tr$mean[1], 2)), as.character(round(abs_igf_t2_N_tr$mean[2], 2)))

means2 <- c("Mean", " ", as.character(round(il1_t2_N_tr$mean[1], 2)), as.character(round(il1_t2_N_tr$mean[2], 2)),
          " ", as.character(round(il6_t2_N_tr$mean[1], 2)), as.character(round(il6_t2_N_tr$mean[2], 2)),
          " ", as.character(round(tnf_t2_N_tr$mean[1], 2)), as.character(round(tnf_t2_N_tr$mean[2], 2)),
          " ", as.character(round(crp_t2_N_tr$mean[1], 2)), as.character(round(crp_t2_N_tr$mean[2], 2)),
          " ", as.character(round(il12_t2_N_tr$mean[1], 2)), as.character(round(il12_t2_N_tr$mean[2], 2)),
          " ", as.character(round(ifn_t2_N_tr$mean[1], 2)), as.character(round(ifn_t2_N_tr$mean[2], 2)),
          " ", as.character(round(il4_t2_N_tr$mean[1], 2)), as.character(round(il4_t2_N_tr$mean[2], 2)),
          " ", as.character(round(il5_t2_N_tr$mean[1], 2)), as.character(round(il5_t2_N_tr$mean[2], 2)),
          " ", as.character(round(il13_t2_N_tr$mean[1], 2)), as.character(round(il13_t2_N_tr$mean[2], 2)),
          " ", as.character(round(il17_t2_N_tr$mean[1], 2)), as.character(round(il17_t2_N_tr$mean[2], 2)),
          " ", as.character(round(il21_t2_N_tr$mean[1], 2)), as.character(round(il21_t2_N_tr$mean[2], 2)),
          " ", as.character(round(il10_t2_N_tr$mean[1], 2)), as.character(round(il10_t2_N_tr$mean[2], 2)),
          " ", as.character(round(il2_t2_N_tr$mean[1], 2)), as.character(round(il2_t2_N_tr$mean[2], 2)),
          " ", as.character(round(gmc_t2_N_tr$mean[1], 2)), as.character(round(gmc_t2_N_tr$mean[2], 2)),
          " ", as.character(round(agp_t2_N_tr$mean[1], 2)), as.character(round(agp_t2_N_tr$mean[2], 2)),
          " ", as.character(round(igf_t2_N_tr$mean[1], 2)), as.character(round(igf_t2_N_tr$mean[2], 2)))

sds2 <- c("SD", " ", as.character(round(il1_t2_N_tr$sd[1], 2)), as.character(round(il1_t2_N_tr$sd[2], 2)),
        " ", as.character(round(il6_t2_N_tr$sd[1], 2)), as.character(round(il6_t2_N_tr$sd[2], 2)),
        " ", as.character(round(tnf_t2_N_tr$sd[1], 2)), as.character(round(tnf_t2_N_tr$sd[2], 2)),
        " ", as.character(round(crp_t2_N_tr$sd[1], 2)), as.character(round(crp_t2_N_tr$sd[2], 2)),
        " ", as.character(round(il12_t2_N_tr$sd[1], 2)), as.character(round(il12_t2_N_tr$sd[2], 2)),
        " ", as.character(round(ifn_t2_N_tr$sd[1], 2)), as.character(round(ifn_t2_N_tr$sd[2], 2)),
        " ", as.character(round(il4_t2_N_tr$sd[1], 2)), as.character(round(il4_t2_N_tr$sd[2], 2)),
        " ", as.character(round(il5_t2_N_tr$sd[1], 2)), as.character(round(il5_t2_N_tr$sd[2], 2)),
        " ", as.character(round(il13_t2_N_tr$sd[1], 2)), as.character(round(il13_t2_N_tr$sd[2], 2)),
        " ", as.character(round(il17_t2_N_tr$sd[1], 2)), as.character(round(il17_t2_N_tr$sd[2], 2)),
        " ", as.character(round(il21_t2_N_tr$sd[1], 2)), as.character(round(il21_t2_N_tr$sd[2], 2)),
        " ", as.character(round(il10_t2_N_tr$sd[1], 2)), as.character(round(il10_t2_N_tr$sd[2], 2)),
        " ", as.character(round(il2_t2_N_tr$sd[1], 2)), as.character(round(il2_t2_N_tr$sd[2], 2)),
        " ", as.character(round(gmc_t2_N_tr$sd[1], 2)), as.character(round(gmc_t2_N_tr$sd[2], 2)),
        " ", as.character(round(agp_t2_N_tr$sd[1], 2)), as.character(round(agp_t2_N_tr$sd[2], 2)),
        " ", as.character(round(igf_t2_N_tr$sd[1], 2)), as.character(round(igf_t2_N_tr$sd[2], 2)))

rawpvals2 <- c("P-value", " ", " ", rawpval(t2_il1_unadj_L[6]), " ", " ", rawpval(t2_il6_unadj_L[6]), " ", " ", rawpval(t2_tnf_unadj_L[6]),
            " ", " ", rawpval(t2_crp_unadj_L[6]), " ", " ", rawpval(t2_il12_unadj_L[6]), " ", " ", rawpval(t2_ifn_unadj_L[6]),
            " ", " ", rawpval(t2_il4_unadj_L[6]), " ", " ", rawpval(t2_il5_unadj_L[6]), " ", " ", rawpval(t2_il13_unadj_L[6]),
            " ", " ", rawpval(t2_il17_unadj_L[6]), " ", " ", rawpval(t2_il21_unadj_L[6]), " ", " ", rawpval(t2_il10_unadj_L[6]),
            " ", " ", rawpval(t2_il2_unadj_L[6]), " ", " ", rawpval(t2_gmc_unadj_L[6]), " ", " ", rawpval(t2_agp_unadj_L[6]),
            " ", " ", rawpval(t2_igf_unadj_L[6]))

pvals2 <- c("Bonferroni P-value", " ", " ", bonpval(t2_il1_unadj_L[6]), " ", " ", bonpval(t2_il6_unadj_L[6]), " ", " ", bonpval(t2_tnf_unadj_L[6]),
            " ", " ", bonpval(t2_crp_unadj_L[6]), " ", " ", bonpval(t2_il12_unadj_L[6]), " ", " ", bonpval(t2_ifn_unadj_L[6]),
            " ", " ", bonpval(t2_il4_unadj_L[6]), " ", " ", bonpval(t2_il5_unadj_L[6]), " ", " ", bonpval(t2_il13_unadj_L[6]),
            " ", " ", bonpval(t2_il17_unadj_L[6]), " ", " ", bonpval(t2_il21_unadj_L[6]), " ", " ", bonpval(t2_il10_unadj_L[6]),
            " ", " ", bonpval(t2_il2_unadj_L[6]), " ", " ", bonpval(t2_gmc_unadj_L[6]), " ", " ", bonpval(t2_agp_unadj_L[6]),
            " ", " ", bonpval(t2_igf_unadj_L[6]))

t2_agp_unadj_L <- round(t2_agp_unadj_L, 2)
t2_crp_unadj_L <- round(t2_crp_unadj_L, 2)
t2_gmc_unadj_L <- round(t2_gmc_unadj_L, 2)
t2_ifn_unadj_L <- round(t2_ifn_unadj_L, 2)
t2_igf_unadj_L <- round(t2_igf_unadj_L, 2)
t2_il1_unadj_L <- round(t2_il1_unadj_L, 2)
t2_il10_unadj_L <- round(t2_il10_unadj_L, 2)
t2_il12_unadj_L <- round(t2_il12_unadj_L, 2)
t2_il13_unadj_L <- round(t2_il13_unadj_L, 2)
t2_il17_unadj_L <- round(t2_il17_unadj_L, 2)
t2_il2_unadj_L <- round(t2_il2_unadj_L, 2)
t2_il21_unadj_L <- round(t2_il21_unadj_L, 2)
t2_il4_unadj_L <- round(t2_il4_unadj_L, 2)
t2_il5_unadj_L <- round(t2_il5_unadj_L, 2)
t2_il6_unadj_L <- round(t2_il6_unadj_L, 2)
t2_tnf_unadj_L <- round(t2_tnf_unadj_L, 2)

unadjs2 <- c("95% CI", " ", " ", paste(t2_il1_unadj_L[1], " (", t2_il1_unadj_L[2], ", ", t2_il1_unadj_L[3], ")", sep=""),
           " ", " ", paste(t2_il6_unadj_L[1], " (", t2_il6_unadj_L[2], ", ", t2_il6_unadj_L[3], ")", sep=""),
           " ", " ", paste(t2_tnf_unadj_L[1], " (", t2_tnf_unadj_L[2], ", ", t2_tnf_unadj_L[3], ")", sep=""),
           " ", " ", paste(t2_crp_unadj_L[1], " (", t2_crp_unadj_L[2], ", ", t2_crp_unadj_L[3], ")", sep=""),
           " ", " ", paste(t2_il12_unadj_L[1], " (", t2_il12_unadj_L[2], ", ", t2_il12_unadj_L[3], ")", sep=""),
           " ", " ", paste(t2_ifn_unadj_L[1], " (", t2_ifn_unadj_L[2], ", ", t2_ifn_unadj_L[3], ")", sep=""),
           " ", " ", paste(t2_il4_unadj_L[1], " (", t2_il4_unadj_L[2], ", ", t2_il4_unadj_L[3], ")", sep=""),
           " ", " ", paste(t2_il5_unadj_L[1], " (", t2_il5_unadj_L[2], ", ", t2_il5_unadj_L[3], ")", sep=""),
           " ", " ", paste(t2_il13_unadj_L[1], " (", t2_il13_unadj_L[2], ", ", t2_il13_unadj_L[3], ")", sep=""),
           " ", " ", paste(t2_il17_unadj_L[1], " (", t2_il17_unadj_L[2], ", ", t2_il17_unadj_L[3], ")", sep=""),
           " ", " ", paste(t2_il21_unadj_L[1], " (", t2_il21_unadj_L[2], ", ", t2_il21_unadj_L[3], ")", sep=""),
           " ", " ", paste(t2_il10_unadj_L[1], " (", t2_il10_unadj_L[2], ", ", t2_il10_unadj_L[3], ")", sep=""),
           " ", " ", paste(t2_il2_unadj_L[1], " (", t2_il2_unadj_L[2], ", ", t2_il2_unadj_L[3], ")", sep=""),
           " ", " ", paste(t2_gmc_unadj_L[1], " (", t2_gmc_unadj_L[2], ", ", t2_gmc_unadj_L[3], ")", sep=""),
           " ", " ", paste(t2_agp_unadj_L[1], " (", t2_agp_unadj_L[2], ", ", t2_agp_unadj_L[3], ")", sep=""),
           " ", " ", paste(t2_igf_unadj_L[1], " (", t2_igf_unadj_L[2], ", ", t2_igf_unadj_L[3], ")", sep="")) 

rawasapvals2 <- c("P-value", " ", " ", rawpval(t2_il1_adj_sex_age_L[6]), " ", " ", rawpval(t2_il6_adj_sex_age_L[6]), " ", " ", rawpval(t2_tnf_adj_sex_age_L[6]),
               " ", " ", rawpval(t2_crp_adj_sex_age_L[6]), " ", " ", rawpval(t2_il12_adj_sex_age_L[6]), " ", " ", rawpval(t2_ifn_adj_sex_age_L[6]),
               " ", " ", rawpval(t2_il4_adj_sex_age_L[6]), " ", " ", rawpval(t2_il5_adj_sex_age_L[6]), " ", " ", rawpval(t2_il13_adj_sex_age_L[6]),
               " ", " ", rawpval(t2_il17_adj_sex_age_L[6]), " ", " ", rawpval(t2_il21_adj_sex_age_L[6]), " ", " ", rawpval(t2_il10_adj_sex_age_L[6]),
               " ", " ", rawpval(t2_il2_adj_sex_age_L[6]), " ", " ", rawpval(t2_gmc_adj_sex_age_L[6]), " ", " ", rawpval(t2_agp_adj_sex_age_L[6]),
               " ", " ", rawpval(t2_igf_adj_sex_age_L[6]))

asapvals2 <- c("Bonferroni P-value", " ", " ", bonpval(t2_il1_adj_sex_age_L[6]), " ", " ", bonpval(t2_il6_adj_sex_age_L[6]), " ", " ", bonpval(t2_tnf_adj_sex_age_L[6]),
               " ", " ", bonpval(t2_crp_adj_sex_age_L[6]), " ", " ", bonpval(t2_il12_adj_sex_age_L[6]), " ", " ", bonpval(t2_ifn_adj_sex_age_L[6]),
               " ", " ", bonpval(t2_il4_adj_sex_age_L[6]), " ", " ", bonpval(t2_il5_adj_sex_age_L[6]), " ", " ", bonpval(t2_il13_adj_sex_age_L[6]),
               " ", " ", bonpval(t2_il17_adj_sex_age_L[6]), " ", " ", bonpval(t2_il21_adj_sex_age_L[6]), " ", " ", bonpval(t2_il10_adj_sex_age_L[6]),
               " ", " ", bonpval(t2_il2_adj_sex_age_L[6]), " ", " ", bonpval(t2_gmc_adj_sex_age_L[6]), " ", " ", bonpval(t2_agp_adj_sex_age_L[6]),
               " ", " ", bonpval(t2_igf_adj_sex_age_L[6]))

t2_agp_adj_sex_age_L <- round(t2_agp_adj_sex_age_L, 2)
t2_crp_adj_sex_age_L <- round(t2_crp_adj_sex_age_L, 2)
t2_gmc_adj_sex_age_L <- round(t2_gmc_adj_sex_age_L, 2)
t2_ifn_adj_sex_age_L <- round(t2_ifn_adj_sex_age_L, 2)
t2_igf_adj_sex_age_L <- round(t2_igf_adj_sex_age_L, 2)
t2_il1_adj_sex_age_L <- round(t2_il1_adj_sex_age_L, 2)
t2_il10_adj_sex_age_L <- round(t2_il10_adj_sex_age_L, 2)
t2_il12_adj_sex_age_L <- round(t2_il12_adj_sex_age_L, 2)
t2_il13_adj_sex_age_L <- round(t2_il13_adj_sex_age_L, 2)
t2_il17_adj_sex_age_L <- round(t2_il17_adj_sex_age_L, 2)
t2_il2_adj_sex_age_L <- round(t2_il2_adj_sex_age_L, 2)
t2_il21_adj_sex_age_L <- round(t2_il21_adj_sex_age_L, 2)
t2_il4_adj_sex_age_L <- round(t2_il4_adj_sex_age_L, 2)
t2_il5_adj_sex_age_L <- round(t2_il5_adj_sex_age_L, 2)
t2_il6_adj_sex_age_L <- round(t2_il6_adj_sex_age_L, 2)
t2_tnf_adj_sex_age_L <- round(t2_tnf_adj_sex_age_L, 2)

asadjs2 <- c("95% CI", " ", " ", paste(t2_il1_adj_sex_age_L[1], " (", t2_il1_adj_sex_age_L[2], ", ", t2_il1_adj_sex_age_L[3], ")", sep=""),
           " ", " ", paste(t2_il6_adj_sex_age_L[1], " (", t2_il6_adj_sex_age_L[2], ", ", t2_il6_adj_sex_age_L[3], ")", sep=""),
           " ", " ", paste(t2_tnf_adj_sex_age_L[1], " (", t2_tnf_adj_sex_age_L[2], ", ", t2_tnf_adj_sex_age_L[3], ")", sep=""),
           " ", " ", paste(t2_crp_adj_sex_age_L[1], " (", t2_crp_adj_sex_age_L[2], ", ", t2_crp_adj_sex_age_L[3], ")", sep=""),
           " ", " ", paste(t2_il12_adj_sex_age_L[1], " (", t2_il12_adj_sex_age_L[2], ", ", t2_il12_adj_sex_age_L[3], ")", sep=""),
           " ", " ", paste(t2_ifn_adj_sex_age_L[1], " (", t2_ifn_adj_sex_age_L[2], ", ", t2_ifn_adj_sex_age_L[3], ")", sep=""),
           " ", " ", paste(t2_il4_adj_sex_age_L[1], " (", t2_il4_adj_sex_age_L[2], ", ", t2_il4_adj_sex_age_L[3], ")", sep=""),
           " ", " ", paste(t2_il5_adj_sex_age_L[1], " (", t2_il5_adj_sex_age_L[2], ", ", t2_il5_adj_sex_age_L[3], ")", sep=""),
           " ", " ", paste(t2_il13_adj_sex_age_L[1], " (", t2_il13_adj_sex_age_L[2], ", ", t2_il13_adj_sex_age_L[3], ")", sep=""),
           " ", " ", paste(t2_il17_adj_sex_age_L[1], " (", t2_il17_adj_sex_age_L[2], ", ", t2_il17_adj_sex_age_L[3], ")", sep=""),
           " ", " ", paste(t2_il21_adj_sex_age_L[1], " (", t2_il21_adj_sex_age_L[2], ", ", t2_il21_adj_sex_age_L[3], ")", sep=""),
           " ", " ", paste(t2_il10_adj_sex_age_L[1], " (", t2_il10_adj_sex_age_L[2], ", ", t2_il10_adj_sex_age_L[3], ")", sep=""),
           " ", " ", paste(t2_il2_adj_sex_age_L[1], " (", t2_il2_adj_sex_age_L[2], ", ", t2_il2_adj_sex_age_L[3], ")", sep=""),
           " ", " ", paste(t2_gmc_adj_sex_age_L[1], " (", t2_gmc_adj_sex_age_L[2], ", ", t2_gmc_adj_sex_age_L[3], ")", sep=""),
           " ", " ", paste(t2_agp_adj_sex_age_L[1], " (", t2_agp_adj_sex_age_L[2], ", ", t2_agp_adj_sex_age_L[3], ")", sep=""),
           " ", " ", paste(t2_igf_adj_sex_age_L[1], " (", t2_igf_adj_sex_age_L[2], ", ", t2_igf_adj_sex_age_L[3], ")", sep=""))  

rawadjpvals2 <- c("P-value", " ", " ", rawpval(t2_il1_adj_L[6]), " ", " ", rawpval(t2_il6_adj_L[6]), " ", " ", rawpval(t2_tnf_adj_L[6]),
               " ", " ", rawpval(t2_crp_adj_L[6]), " ", " ", rawpval(t2_il12_adj_L[6]), " ", " ", rawpval(t2_ifn_adj_L[6]),
               " ", " ", rawpval(t2_il4_adj_L[6]), " ", " ", rawpval(t2_il5_adj_L[6]), " ", " ", rawpval(t2_il13_adj_L[6]),
               " ", " ", rawpval(t2_il17_adj_L[6]), " ", " ", rawpval(t2_il21_adj_L[6]), " ", " ", rawpval(t2_il10_adj_L[6]),
               " ", " ", rawpval(t2_il2_adj_L[6]), " ", " ", rawpval(t2_gmc_adj_L[6]), " ", " ", rawpval(t2_agp_adj_L[6]),
               " ", " ", rawpval(t2_igf_adj_L[6]))

adjpvals2 <- c("Bonferroni P-value", " ", " ", bonpval(t2_il1_adj_L[6]), " ", " ", bonpval(t2_il6_adj_L[6]), " ", " ", bonpval(t2_tnf_adj_L[6]),
               " ", " ", bonpval(t2_crp_adj_L[6]), " ", " ", bonpval(t2_il12_adj_L[6]), " ", " ", bonpval(t2_ifn_adj_L[6]),
               " ", " ", bonpval(t2_il4_adj_L[6]), " ", " ", bonpval(t2_il5_adj_L[6]), " ", " ", bonpval(t2_il13_adj_L[6]),
               " ", " ", bonpval(t2_il17_adj_L[6]), " ", " ", bonpval(t2_il21_adj_L[6]), " ", " ", bonpval(t2_il10_adj_L[6]),
               " ", " ", bonpval(t2_il2_adj_L[6]), " ", " ", bonpval(t2_gmc_adj_L[6]), " ", " ", bonpval(t2_agp_adj_L[6]),
               " ", " ", bonpval(t2_igf_adj_L[6]))

t2_agp_adj_L <- round(t2_agp_adj_L, 2)
t2_crp_adj_L <- round(t2_crp_adj_L, 2)
t2_gmc_adj_L <- round(t2_gmc_adj_L, 2)
t2_ifn_adj_L <- round(t2_ifn_adj_L, 2)
t2_igf_adj_L <- round(t2_igf_adj_L, 2)
t2_il1_adj_L <- round(t2_il1_adj_L, 2)
t2_il10_adj_L <- round(t2_il10_adj_L, 2)
t2_il12_adj_L <- round(t2_il12_adj_L, 2)
t2_il13_adj_L <- round(t2_il13_adj_L, 2)
t2_il17_adj_L <- round(t2_il17_adj_L, 2)
t2_il2_adj_L <- round(t2_il2_adj_L, 2)
t2_il21_adj_L <- round(t2_il21_adj_L, 2)
t2_il4_adj_L <- round(t2_il4_adj_L, 2)
t2_il5_adj_L <- round(t2_il5_adj_L, 2)
t2_il6_adj_L <- round(t2_il6_adj_L, 2)
t2_tnf_adj_L <- round(t2_tnf_adj_L, 2)

adjs2 <- c("95% CI", " ", " ", paste(t2_il1_adj_L[1], " (", t2_il1_adj_L[2], ", ", t2_il1_adj_L[3], ")", sep=""),
         " ", " ", paste(t2_il6_adj_L[1], " (", t2_il6_adj_L[2], ", ", t2_il6_adj_L[3], ")", sep=""),
         " ", " ", paste(t2_tnf_adj_L[1], " (", t2_tnf_adj_L[2], ", ", t2_tnf_adj_L[3], ")", sep=""),
         " ", " ", paste(t2_crp_adj_L[1], " (", t2_crp_adj_L[2], ", ", t2_crp_adj_L[3], ")", sep=""),
         " ", " ", paste(t2_il12_adj_L[1], " (", t2_il12_adj_L[2], ", ", t2_il12_adj_L[3], ")", sep=""),
         " ", " ", paste(t2_ifn_adj_L[1], " (", t2_ifn_adj_L[2], ", ", t2_ifn_adj_L[3], ")", sep=""),
         " ", " ", paste(t2_il4_adj_L[1], " (", t2_il4_adj_L[2], ", ", t2_il4_adj_L[3], ")", sep=""),
         " ", " ", paste(t2_il5_adj_L[1], " (", t2_il5_adj_L[2], ", ", t2_il5_adj_L[3], ")", sep=""),
         " ", " ", paste(t2_il13_adj_L[1], " (", t2_il13_adj_L[2], ", ", t2_il13_adj_L[3], ")", sep=""),
         " ", " ", paste(t2_il17_adj_L[1], " (", t2_il17_adj_L[2], ", ", t2_il17_adj_L[3], ")", sep=""),
         " ", " ", paste(t2_il21_adj_L[1], " (", t2_il21_adj_L[2], ", ", t2_il21_adj_L[3], ")", sep=""),
         " ", " ", paste(t2_il10_adj_L[1], " (", t2_il10_adj_L[2], ", ", t2_il10_adj_L[3], ")", sep=""),
         " ", " ", paste(t2_il2_adj_L[1], " (", t2_il2_adj_L[2], ", ", t2_il2_adj_L[3], ")", sep=""),
         " ", " ", paste(t2_gmc_adj_L[1], " (", t2_gmc_adj_L[2], ", ", t2_gmc_adj_L[3], ")", sep=""),
         " ", " ", paste(t2_agp_adj_L[1], " (", t2_agp_adj_L[2], ", ", t2_agp_adj_L[3], ")", sep=""),
         " ", " ", paste(t2_igf_adj_L[1], " (", t2_igf_adj_L[2], ", ", t2_igf_adj_L[3], ")", sep="")) 

ipcws2 <-c("95% CI", " ", " ", maketblvalue(il1_t2_adj_ipcw_L$`unlist(il1_t2_adj_ipcw$estimates$ATE)`),
          " ", " ", maketblvalue(il6_t2_adj_ipcw_L$`unlist(il6_t2_adj_ipcw$estimates$ATE)`),
          " ", " ", maketblvalue(tnf_t2_adj_ipcw_L$`unlist(tnf_t2_adj_ipcw$estimates$ATE)`), 
          " ", " ", maketblvalue(crp_t2_adj_ipcw_L$`unlist(crp_t2_adj_ipcw$estimates$ATE)`),
          " ", " ", maketblvalue(il12_t2_adj_ipcw_L$`unlist(il12_t2_adj_ipcw$estimates$ATE)`),
          " ", " ", maketblvalue(ifn_t2_adj_ipcw_L$`unlist(ifn_t2_adj_ipcw$estimates$ATE)`),
          " ", " ", maketblvalue(il4_t2_adj_ipcw_L$`unlist(il4_t2_adj_ipcw$estimates$ATE)`),
          " ", " ", maketblvalue(il5_t2_adj_ipcw_L$`unlist(il5_t2_adj_ipcw$estimates$ATE)`),
          " ", " ", maketblvalue(il13_t2_adj_ipcw_L$`unlist(il13_t2_adj_ipcw$estimates$ATE)`),
          " ", " ", maketblvalue(il17_t2_adj_ipcw_L$`unlist(il17_t2_adj_ipcw$estimates$ATE)`),
          " ", " ", maketblvalue(il21_t2_adj_ipcw_L$`unlist(il21_t2_adj_ipcw$estimates$ATE)`), 
          " ", " ", maketblvalue(il10_t2_adj_ipcw_L$`unlist(il10_t2_adj_ipcw$estimates$ATE)`),
          " ", " ", maketblvalue(il2_t2_adj_ipcw_L$`unlist(il2_t2_adj_ipcw$estimates$ATE)`),
          " ", " ", maketblvalue(gmc_t2_adj_ipcw_L$`unlist(gmc_t2_adj_ipcw$estimates$ATE)`),
          " ", " ", maketblvalue(agp_t2_adj_ipcw_L$`unlist(agp_t2_adj_ipcw$estimates$ATE)`),
          " ", " ", maketblvalue(igf_t2_adj_ipcw_L$`unlist(igf_t2_adj_ipcw$estimates$ATE)`))

rawipcwpvals2 <- c("P-value", " ", " ", rawpval(il1_t2_adj_ipcw_L[5,1]),
                " ", " ", rawpval(il6_t2_adj_ipcw_L[5,1]),
                " ", " ", rawpval(tnf_t2_adj_ipcw_L[5,1]), 
                " ", " ", rawpval(crp_t2_adj_ipcw_L[5,1]),
                " ", " ", rawpval(il12_t2_adj_ipcw_L[5,1]),
                " ", " ", rawpval(ifn_t2_adj_ipcw_L[5,1]),
                " ", " ", rawpval(il4_t2_adj_ipcw_L[5,1]),
                " ", " ", rawpval(il5_t2_adj_ipcw_L[5,1]),
                " ", " ", rawpval(il13_t2_adj_ipcw_L[5,1]),
                " ", " ", rawpval(il17_t2_adj_ipcw_L[5,1]),
                " ", " ", rawpval(il21_t2_adj_ipcw_L[5,1]), 
                " ", " ", rawpval(il10_t2_adj_ipcw_L[5,1]),
                " ", " ", rawpval(il2_t2_adj_ipcw_L[5,1]),
                " ", " ", rawpval(gmc_t2_adj_ipcw_L[5,1]),
                " ", " ", rawpval(agp_t2_adj_ipcw_L[5,1]),
                " ", " ", rawpval(igf_t2_adj_ipcw_L[5,1]))

ipcwpvals2 <- c("Bonferroni P-value", " ", " ", bonpval(il1_t2_adj_ipcw_L[5,1]),
                " ", " ", bonpval(il6_t2_adj_ipcw_L[5,1]),
                " ", " ", bonpval(tnf_t2_adj_ipcw_L[5,1]), 
                " ", " ", bonpval(crp_t2_adj_ipcw_L[5,1]),
                " ", " ", bonpval(il12_t2_adj_ipcw_L[5,1]),
                " ", " ", bonpval(ifn_t2_adj_ipcw_L[5,1]),
                " ", " ", bonpval(il4_t2_adj_ipcw_L[5,1]),
                " ", " ", bonpval(il5_t2_adj_ipcw_L[5,1]),
                " ", " ", bonpval(il13_t2_adj_ipcw_L[5,1]),
                " ", " ", bonpval(il17_t2_adj_ipcw_L[5,1]),
                " ", " ", bonpval(il21_t2_adj_ipcw_L[5,1]), 
                " ", " ", bonpval(il10_t2_adj_ipcw_L[5,1]),
                " ", " ", bonpval(il2_t2_adj_ipcw_L[5,1]),
                " ", " ", bonpval(gmc_t2_adj_ipcw_L[5,1]),
                " ", " ", bonpval(agp_t2_adj_ipcw_L[5,1]),
                " ", " ", bonpval(igf_t2_adj_ipcw_L[5,1]))

# Table S2: Effect of intervention on individual immune status and growth factor measurements at age 14 months
tbls2 <- data.table(
  " " = outcomes2,
  " " = Ns2, 
  " " = absmeans2,
  " " = abssds2,
  " " = means2, 
  " " = sds2,
  "Unadjusted difference: Intervention vs. Control" = unadjs2,
  " " = rawpvals2,
  " " = pvals2,
  "Age- and sex- adjusted difference: Intervention vs. Control" = asadjs2, 
  " " = rawasapvals2,
  " " = asapvals2,
  "Fully adjusted difference: Intervention vs. Control" = adjs2,
  " " = rawadjpvals2,
  " " = adjpvals2,
  "IPCW adjusted difference: Intervention vs. Control" = ipcws2,
  " " = rawipcwpvals2,
  " " = ipcwpvals2)

write.csv(tbls2, file=here('tables/supplementary/immune_supptable2.csv'))
print(xtable(tbls2), type="html", file=here("tables/supplementary/immune_supptable2.html"))



#### TABLE S3 ####

rawpvals3 <- c("P-value", " ", " ", rawpval(t2_ratio_il1_il10_adj_sex_age_L[6]), " ", " ", rawpval(t2_ratio_il6_il10_adj_sex_age_L[6]), " ", " ", rawpval(t2_ratio_tnf_il10_adj_sex_age_L[6]), 
            " ", " ", rawpval(t2_ratio_il12_il10_adj_sex_age_L[6]), " ", " ", rawpval(t2_ratio_ifn_il10_adj_sex_age_L[6]), " ", " ", rawpval(t2_ratio_il4_il10_adj_sex_age_L[6]), 
            " ", " ", rawpval(t2_ratio_il5_il10_adj_sex_age_L[6]), " ", " ", rawpval(t2_ratio_il13_il10_adj_sex_age_L[6]), " ", " ", rawpval(t2_ratio_il17_il10_adj_sex_age_L[6]), 
            " ", " ", rawpval(t2_ratio_il21_il10_adj_sex_age_L[6]), " ", " ", rawpval(t2_ratio_il2_il10_adj_sex_age_L[6]), " ", " ", rawpval(t2_ratio_gmc_il10_adj_sex_age_L[6]), 
            " ", " ", rawpval(t2_ratio_il12_il4_adj_sex_age_L[6]), " ", " ", rawpval(t2_ratio_ifn_il4_adj_sex_age_L[6]), " ", " ", rawpval(t2_ratio_il12_il5_adj_sex_age_L[6]), 
            " ", " ", rawpval(t2_ratio_ifn_il5_adj_sex_age_L[6]), " ", " ", rawpval(t2_ratio_il12_il13_adj_sex_age_L[6]), " ", " ", rawpval(t2_ratio_ifn_il13_adj_sex_age_L[6]), 
            " ", " ", rawpval(t2_ratio_il12_il17_adj_sex_age_L[6]), " ", " ", rawpval(t2_ratio_ifn_il17_adj_sex_age_L[6]), " ", " ", rawpval(t2_ratio_il12_il21_adj_sex_age_L[6]), 
            " ", " ", rawpval(t2_ratio_ifn_il21_adj_sex_age_L[6]), " ", " ", rawpval(t2_ratio_pro_il10_adj_sex_age_L[6]), " ", " ", rawpval(t2_ratio_th1_il10_adj_sex_age_L[6]), 
            " ", " ", rawpval(t2_ratio_th2_il10_adj_sex_age_L[6]), " ", " ", rawpval(t2_ratio_th17_il10_adj_sex_age_L[6]), " ", " ", rawpval(t2_ratio_th1_th2_adj_sex_age_L[6]), " ", " ", rawpval(t2_ratio_th1_th17_adj_sex_age_L[6])
)

pvals3 <- c("Bonferroni P-value", " ", " ", bonpval(t2_ratio_il1_il10_adj_sex_age_L[6]), " ", " ", bonpval(t2_ratio_il6_il10_adj_sex_age_L[6]), " ", " ", bonpval(t2_ratio_tnf_il10_adj_sex_age_L[6]), 
            " ", " ", bonpval(t2_ratio_il12_il10_adj_sex_age_L[6]), " ", " ", bonpval(t2_ratio_ifn_il10_adj_sex_age_L[6]), " ", " ", bonpval(t2_ratio_il4_il10_adj_sex_age_L[6]), 
            " ", " ", bonpval(t2_ratio_il5_il10_adj_sex_age_L[6]), " ", " ", bonpval(t2_ratio_il13_il10_adj_sex_age_L[6]), " ", " ", bonpval(t2_ratio_il17_il10_adj_sex_age_L[6]), 
            " ", " ", bonpval(t2_ratio_il21_il10_adj_sex_age_L[6]), " ", " ", bonpval(t2_ratio_il2_il10_adj_sex_age_L[6]), " ", " ", bonpval(t2_ratio_gmc_il10_adj_sex_age_L[6]), 
            " ", " ", bonpval(t2_ratio_il12_il4_adj_sex_age_L[6]), " ", " ", bonpval(t2_ratio_ifn_il4_adj_sex_age_L[6]), " ", " ", bonpval(t2_ratio_il12_il5_adj_sex_age_L[6]), 
            " ", " ", bonpval(t2_ratio_ifn_il5_adj_sex_age_L[6]), " ", " ", bonpval(t2_ratio_il12_il13_adj_sex_age_L[6]), " ", " ", bonpval(t2_ratio_ifn_il13_adj_sex_age_L[6]), 
            " ", " ", bonpval(t2_ratio_il12_il17_adj_sex_age_L[6]), " ", " ", bonpval(t2_ratio_ifn_il17_adj_sex_age_L[6]), " ", " ", bonpval(t2_ratio_il12_il21_adj_sex_age_L[6]), 
            " ", " ", bonpval(t2_ratio_ifn_il21_adj_sex_age_L[6]), " ", " ", bonpval(t2_ratio_pro_il10_adj_sex_age_L[6]), " ", " ", bonpval(t2_ratio_th1_il10_adj_sex_age_L[6]), 
            " ", " ", bonpval(t2_ratio_th2_il10_adj_sex_age_L[6]), " ", " ", bonpval(t2_ratio_th17_il10_adj_sex_age_L[6]), " ", " ", bonpval(t2_ratio_th1_th2_adj_sex_age_L[6]), " ", " ", bonpval(t2_ratio_th1_th17_adj_sex_age_L[6])
)

t2_ratio_il1_il10_adj_sex_age_L <- round(t2_ratio_il1_il10_adj_sex_age_L, 2)
t2_ratio_il6_il10_adj_sex_age_L <- round(t2_ratio_il6_il10_adj_sex_age_L, 2)
t2_ratio_tnf_il10_adj_sex_age_L <- round(t2_ratio_tnf_il10_adj_sex_age_L, 2)
t2_ratio_il12_il10_adj_sex_age_L <- round(t2_ratio_il12_il10_adj_sex_age_L, 2)
t2_ratio_ifn_il10_adj_sex_age_L <- round(t2_ratio_ifn_il10_adj_sex_age_L, 2)
t2_ratio_il4_il10_adj_sex_age_L <- round(t2_ratio_il4_il10_adj_sex_age_L, 2)
t2_ratio_il5_il10_adj_sex_age_L <- round(t2_ratio_il5_il10_adj_sex_age_L, 2)
t2_ratio_il13_il10_adj_sex_age_L <- round(t2_ratio_il13_il10_adj_sex_age_L, 2)
t2_ratio_il17_il10_adj_sex_age_L <- round(t2_ratio_il17_il10_adj_sex_age_L, 2)
t2_ratio_il21_il10_adj_sex_age_L <- round(t2_ratio_il21_il10_adj_sex_age_L, 2)
t2_ratio_il2_il10_adj_sex_age_L <- round(t2_ratio_il2_il10_adj_sex_age_L, 2)
t2_ratio_gmc_il10_adj_sex_age_L <- round(t2_ratio_gmc_il10_adj_sex_age_L, 2)
t2_ratio_il12_il4_adj_sex_age_L <- round(t2_ratio_il12_il4_adj_sex_age_L, 2)
t2_ratio_ifn_il4_adj_sex_age_L <- round(t2_ratio_ifn_il4_adj_sex_age_L, 2)
t2_ratio_il12_il5_adj_sex_age_L <- round(t2_ratio_il12_il5_adj_sex_age_L, 2)
t2_ratio_ifn_il5_adj_sex_age_L <- round(t2_ratio_ifn_il5_adj_sex_age_L, 2)
t2_ratio_il12_il13_adj_sex_age_L <- round(t2_ratio_il12_il13_adj_sex_age_L, 2)
t2_ratio_ifn_il13_adj_sex_age_L <- round(t2_ratio_ifn_il13_adj_sex_age_L, 2)
t2_ratio_il12_il17_adj_sex_age_L <- round(t2_ratio_il12_il17_adj_sex_age_L, 2)
t2_ratio_ifn_il17_adj_sex_age_L <- round(t2_ratio_ifn_il17_adj_sex_age_L, 2)
t2_ratio_il12_il21_adj_sex_age_L <- round(t2_ratio_il12_il21_adj_sex_age_L, 2)
t2_ratio_ifn_il21_adj_sex_age_L <- round(t2_ratio_ifn_il21_adj_sex_age_L, 2)
t2_ratio_pro_il10_adj_sex_age_L <- round(t2_ratio_pro_il10_adj_sex_age_L, 2)
t2_ratio_th1_il10_adj_sex_age_L <- round(t2_ratio_th1_il10_adj_sex_age_L, 2)
t2_ratio_th2_il10_adj_sex_age_L <- round(t2_ratio_th2_il10_adj_sex_age_L, 2)
t2_ratio_th17_il10_adj_sex_age_L <- round(t2_ratio_th17_il10_adj_sex_age_L, 2)
t2_ratio_th1_th2_adj_sex_age_L <- round(t2_ratio_th1_th2_adj_sex_age_L, 2)
t2_ratio_th1_th17_adj_sex_age_L <- round(t2_ratio_th1_th17_adj_sex_age_L, 2)

asadjtbls3 <- c("95% CI", " ", " ", paste(t2_ratio_il1_il10_adj_sex_age_L[1], " (", t2_ratio_il1_il10_adj_sex_age_L[2], ", ", t2_ratio_il1_il10_adj_sex_age_L[3], ")", sep=""),
               " ", " ", paste(t2_ratio_il6_il10_adj_sex_age_L[1], " (", t2_ratio_il6_il10_adj_sex_age_L[2], ", ", t2_ratio_il6_il10_adj_sex_age_L[3], ")", sep=""),
               " ", " ", paste(t2_ratio_tnf_il10_adj_sex_age_L[1], " (", t2_ratio_tnf_il10_adj_sex_age_L[2], ", ", t2_ratio_tnf_il10_adj_sex_age_L[3], ")", sep=""),
               " ", " ", paste(t2_ratio_il12_il10_adj_sex_age_L[1], " (", t2_ratio_il12_il10_adj_sex_age_L[2], ", ", t2_ratio_il12_il10_adj_sex_age_L[3], ")", sep=""),
               " ", " ", paste(t2_ratio_ifn_il10_adj_sex_age_L[1], " (", t2_ratio_ifn_il10_adj_sex_age_L[2], ", ", t2_ratio_ifn_il10_adj_sex_age_L[3], ")", sep=""),
               " ", " ", paste(t2_ratio_il4_il10_adj_sex_age_L[1], " (", t2_ratio_il4_il10_adj_sex_age_L[2], ", ", t2_ratio_il4_il10_adj_sex_age_L[3], ")", sep=""),
               " ", " ", paste(t2_ratio_il5_il10_adj_sex_age_L[1], " (", t2_ratio_il5_il10_adj_sex_age_L[2], ", ", t2_ratio_il5_il10_adj_sex_age_L[3], ")", sep=""),
               " ", " ", paste(t2_ratio_il13_il10_adj_sex_age_L[1], " (", t2_ratio_il13_il10_adj_sex_age_L[2], ", ", t2_ratio_il13_il10_adj_sex_age_L[3], ")", sep=""),
               " ", " ", paste(t2_ratio_il17_il10_adj_sex_age_L[1], " (", t2_ratio_il17_il10_adj_sex_age_L[2], ", ", t2_ratio_il17_il10_adj_sex_age_L[3], ")", sep=""),
               " ", " ", paste(t2_ratio_il21_il10_adj_sex_age_L[1], " (", t2_ratio_il21_il10_adj_sex_age_L[2], ", ", t2_ratio_il21_il10_adj_sex_age_L[3], ")", sep=""),
               " ", " ", paste(t2_ratio_il2_il10_adj_sex_age_L[1], " (", t2_ratio_il2_il10_adj_sex_age_L[2], ", ", t2_ratio_il2_il10_adj_sex_age_L[3], ")", sep=""),
               " ", " ", paste(t2_ratio_gmc_il10_adj_sex_age_L[1], " (", t2_ratio_gmc_il10_adj_sex_age_L[2], ", ", t2_ratio_gmc_il10_adj_sex_age_L[3], ")", sep=""),
               " ", " ", paste(t2_ratio_il12_il4_adj_sex_age_L[1], " (", t2_ratio_il12_il4_adj_sex_age_L[2], ", ", t2_ratio_il12_il4_adj_sex_age_L[3], ")", sep=""),
               " ", " ", paste(t2_ratio_ifn_il4_adj_sex_age_L[1], " (", t2_ratio_ifn_il4_adj_sex_age_L[2], ", ", t2_ratio_ifn_il4_adj_sex_age_L[3], ")", sep=""),
               " ", " ", paste(t2_ratio_il12_il5_adj_sex_age_L[1], " (", t2_ratio_il12_il5_adj_sex_age_L[2], ", ", t2_ratio_il12_il5_adj_sex_age_L[3], ")", sep=""),
               " ", " ", paste(t2_ratio_ifn_il5_adj_sex_age_L[1], " (", t2_ratio_ifn_il5_adj_sex_age_L[2], ", ", t2_ratio_ifn_il5_adj_sex_age_L[3], ")", sep=""),
               " ", " ", paste(t2_ratio_il12_il13_adj_sex_age_L[1], " (", t2_ratio_il12_il13_adj_sex_age_L[2], ", ", t2_ratio_il12_il13_adj_sex_age_L[3], ")", sep=""),
               " ", " ", paste(t2_ratio_ifn_il13_adj_sex_age_L[1], " (", t2_ratio_ifn_il13_adj_sex_age_L[2], ", ", t2_ratio_ifn_il13_adj_sex_age_L[3], ")", sep=""),
               " ", " ", paste(t2_ratio_il12_il17_adj_sex_age_L[1], " (", t2_ratio_il12_il17_adj_sex_age_L[2], ", ", t2_ratio_il12_il17_adj_sex_age_L[3], ")", sep=""),
               " ", " ", paste(t2_ratio_ifn_il17_adj_sex_age_L[1], " (", t2_ratio_ifn_il17_adj_sex_age_L[2], ", ", t2_ratio_ifn_il17_adj_sex_age_L[3], ")", sep=""),
               " ", " ", paste(t2_ratio_il12_il21_adj_sex_age_L[1], " (", t2_ratio_il12_il21_adj_sex_age_L[2], ", ", t2_ratio_il12_il21_adj_sex_age_L[3], ")", sep=""),
               " ", " ", paste(t2_ratio_ifn_il21_adj_sex_age_L[1], " (", t2_ratio_ifn_il21_adj_sex_age_L[2], ", ", t2_ratio_ifn_il21_adj_sex_age_L[3], ")", sep=""),
               " ", " ", paste(t2_ratio_pro_il10_adj_sex_age_L[1], " (", t2_ratio_pro_il10_adj_sex_age_L[2], ", ", t2_ratio_pro_il10_adj_sex_age_L[3], ")", sep=""),
               " ", " ", paste(t2_ratio_th1_il10_adj_sex_age_L[1], " (", t2_ratio_th1_il10_adj_sex_age_L[2], ", ", t2_ratio_th1_il10_adj_sex_age_L[3], ")", sep=""),
               " ", " ", paste(t2_ratio_th2_il10_adj_sex_age_L[1], " (", t2_ratio_th2_il10_adj_sex_age_L[2], ", ", t2_ratio_th2_il10_adj_sex_age_L[3], ")", sep=""),
               " ", " ", paste(t2_ratio_th17_il10_adj_sex_age_L[1], " (", t2_ratio_th17_il10_adj_sex_age_L[2], ", ", t2_ratio_th17_il10_adj_sex_age_L[3], ")", sep=""),
               " ", " ", paste(t2_ratio_th1_th2_adj_sex_age_L[1], " (", t2_ratio_th1_th2_adj_sex_age_L[2], ", ", t2_ratio_th1_th2_adj_sex_age_L[3], ")", sep=""),
               " ", " ", paste(t2_ratio_th1_th17_adj_sex_age_L[1], " (", t2_ratio_th1_th17_adj_sex_age_L[2], ", ", t2_ratio_th1_th17_adj_sex_age_L[3], ")", sep=""))

ipcws3<-c("95% CI", " ", " ", maketblvalue(ratio_il1_il10_t2_adj_ipcw_L$`unlist(ratio_il1_il10_t2_adj_ipcw$estimates$ATE)`),
          " ", " ", maketblvalue(ratio_il6_il10_t2_adj_ipcw_L$`unlist(ratio_il6_il10_t2_adj_ipcw$estimates$ATE)`),
          " ", " ", maketblvalue(ratio_tnf_il10_t2_adj_ipcw_L$`unlist(ratio_tnf_il10_t2_adj_ipcw$estimates$ATE)`),
          " ", " ", maketblvalue(ratio_il12_il10_t2_adj_ipcw_L$`unlist(ratio_il12_il10_t2_adj_ipcw$estimates$ATE)`),
          " ", " ", maketblvalue(ratio_ifn_il10_t2_adj_ipcw_L$`unlist(ratio_ifn_il10_t2_adj_ipcw$estimates$ATE)`),
          " ", " ", maketblvalue(ratio_il4_il10_t2_adj_ipcw_L$`unlist(ratio_il4_il10_t2_adj_ipcw$estimates$ATE)`),
          " ", " ", maketblvalue(ratio_il5_il10_t2_adj_ipcw_L$`unlist(ratio_il5_il10_t2_adj_ipcw$estimates$ATE)`),
          " ", " ", maketblvalue(ratio_il13_il10_t2_adj_ipcw_L$`unlist(ratio_il13_il10_t2_adj_ipcw$estimates$ATE)`),
          " ", " ", maketblvalue(ratio_il17_il10_t2_adj_ipcw_L$`unlist(ratio_il17_il10_t2_adj_ipcw$estimates$ATE)`),
          " ", " ", maketblvalue(ratio_il21_il10_t2_adj_ipcw_L$`unlist(ratio_il21_il10_t2_adj_ipcw$estimates$ATE)`),
          " ", " ", maketblvalue(ratio_il2_il10_t2_adj_ipcw_L$`unlist(ratio_il2_il10_t2_adj_ipcw$estimates$ATE)`),
          " ", " ", maketblvalue(ratio_gmc_il10_t2_adj_ipcw_L$`unlist(ratio_gmc_il10_t2_adj_ipcw$estimates$ATE)`),
          " ", " ", maketblvalue(ratio_il12_il4_t2_adj_ipcw_L$`unlist(ratio_il12_il4_t2_adj_ipcw$estimates$ATE)`),
          " ", " ", maketblvalue(ratio_ifn_il4_t2_adj_ipcw_L$`unlist(ratio_ifn_il4_t2_adj_ipcw$estimates$ATE)`),
          " ", " ", maketblvalue(ratio_il12_il5_t2_adj_ipcw_L$`unlist(ratio_il12_il5_t2_adj_ipcw$estimates$ATE)`),
          " ", " ", maketblvalue(ratio_ifn_il5_t2_adj_ipcw_L$`unlist(ratio_ifn_il5_t2_adj_ipcw$estimates$ATE)`),
          " ", " ", maketblvalue(ratio_il12_il13_t2_adj_ipcw_L$`unlist(ratio_il12_il13_t2_adj_ipcw$estimates$ATE)`),
          " ", " ", maketblvalue(ratio_ifn_il13_t2_adj_ipcw_L$`unlist(ratio_ifn_il13_t2_adj_ipcw$estimates$ATE)`),
          " ", " ", maketblvalue(ratio_il12_il17_t2_adj_ipcw_L$`unlist(ratio_il12_il17_t2_adj_ipcw$estimates$ATE)`),
          " ", " ", maketblvalue(ratio_ifn_il17_t2_adj_ipcw_L$`unlist(ratio_ifn_il17_t2_adj_ipcw$estimates$ATE)`),
          " ", " ", maketblvalue(ratio_il12_il21_t2_adj_ipcw_L$`unlist(ratio_il12_il21_t2_adj_ipcw$estimates$ATE)`),
          " ", " ", maketblvalue(ratio_ifn_il21_t2_adj_ipcw_L$`unlist(ratio_ifn_il21_t2_adj_ipcw$estimates$ATE)`),
          " ", " ", maketblvalue(ratio_pro_il10_t2_adj_ipcw_L$`unlist(ratio_pro_il10_t2_adj_ipcw$estimates$ATE)`),
          " ", " ", maketblvalue(ratio_th1_il10_t2_adj_ipcw_L$`unlist(ratio_th1_il10_t2_adj_ipcw$estimates$ATE)`),
          " ", " ", maketblvalue(ratio_th2_il10_t2_adj_ipcw_L$`unlist(ratio_th2_il10_t2_adj_ipcw$estimates$ATE)`),
          " ", " ", maketblvalue(ratio_th17_il10_t2_adj_ipcw_L$`unlist(ratio_th17_il10_t2_adj_ipcw$estimates$ATE)`),
          " ", " ", maketblvalue(ratio_th1_th2_t2_adj_ipcw_L$`unlist(ratio_th1_th2_t2_adj_ipcw$estimates$ATE)`),
          " ", " ", maketblvalue(ratio_th1_th17_t2_adj_ipcw_L$`unlist(ratio_th1_th17_t2_adj_ipcw$estimates$ATE)`))

rawipcwpvals3 <- c("P-value", " ", " ", rawpval(ratio_il1_il10_t2_adj_ipcw_L[5,1]),
          " ", " ", rawpval(ratio_il6_il10_t2_adj_ipcw_L[5,1]),
          " ", " ", rawpval(ratio_tnf_il10_t2_adj_ipcw_L[5,1]),
          " ", " ", rawpval(ratio_il12_il10_t2_adj_ipcw_L[5,1]),
          " ", " ", rawpval(ratio_ifn_il10_t2_adj_ipcw_L[5,1]),
          " ", " ", rawpval(ratio_il4_il10_t2_adj_ipcw_L[5,1]),
          " ", " ", rawpval(ratio_il5_il10_t2_adj_ipcw_L[5,1]),
          " ", " ", rawpval(ratio_il13_il10_t2_adj_ipcw_L[5,1]),
          " ", " ", rawpval(ratio_il17_il10_t2_adj_ipcw_L[5,1]),
          " ", " ", rawpval(ratio_il21_il10_t2_adj_ipcw_L[5,1]),
          " ", " ", rawpval(ratio_il2_il10_t2_adj_ipcw_L[5,1]),
          " ", " ", rawpval(ratio_gmc_il10_t2_adj_ipcw_L[5,1]),
          " ", " ", rawpval(ratio_il12_il4_t2_adj_ipcw_L[5,1]),
          " ", " ", rawpval(ratio_ifn_il4_t2_adj_ipcw_L[5,1]),
          " ", " ", rawpval(ratio_il12_il5_t2_adj_ipcw_L[5,1]),
          " ", " ", rawpval(ratio_ifn_il5_t2_adj_ipcw_L[5,1]),
          " ", " ", rawpval(ratio_il12_il13_t2_adj_ipcw_L[5,1]),
          " ", " ", rawpval(ratio_ifn_il13_t2_adj_ipcw_L[5,1]),
          " ", " ", rawpval(ratio_il12_il17_t2_adj_ipcw_L[5,1]),
          " ", " ", rawpval(ratio_ifn_il17_t2_adj_ipcw_L[5,1]),
          " ", " ", rawpval(ratio_il12_il21_t2_adj_ipcw_L[5,1]),
          " ", " ", rawpval(ratio_ifn_il21_t2_adj_ipcw_L[5,1]),
          " ", " ", rawpval(ratio_pro_il10_t2_adj_ipcw_L[5,1]),
          " ", " ", rawpval(ratio_th1_il10_t2_adj_ipcw_L[5,1]),
          " ", " ", rawpval(ratio_th2_il10_t2_adj_ipcw_L[5,1]),
          " ", " ", rawpval(ratio_th17_il10_t2_adj_ipcw_L[5,1]),
          " ", " ", rawpval(ratio_th1_th2_t2_adj_ipcw_L[5,1]),
          " ", " ", rawpval(ratio_th1_th17_t2_adj_ipcw_L[5,1]))

ipcwpvals3 <- c("Bonferroni P-value", " ", " ", bonpval(ratio_il1_il10_t2_adj_ipcw_L[5,1]),
          " ", " ", bonpval(ratio_il6_il10_t2_adj_ipcw_L[5,1]),
          " ", " ", bonpval(ratio_tnf_il10_t2_adj_ipcw_L[5,1]),
          " ", " ", bonpval(ratio_il12_il10_t2_adj_ipcw_L[5,1]),
          " ", " ", bonpval(ratio_ifn_il10_t2_adj_ipcw_L[5,1]),
          " ", " ", bonpval(ratio_il4_il10_t2_adj_ipcw_L[5,1]),
          " ", " ", bonpval(ratio_il5_il10_t2_adj_ipcw_L[5,1]),
          " ", " ", bonpval(ratio_il13_il10_t2_adj_ipcw_L[5,1]),
          " ", " ", bonpval(ratio_il17_il10_t2_adj_ipcw_L[5,1]),
          " ", " ", bonpval(ratio_il21_il10_t2_adj_ipcw_L[5,1]),
          " ", " ", bonpval(ratio_il2_il10_t2_adj_ipcw_L[5,1]),
          " ", " ", bonpval(ratio_gmc_il10_t2_adj_ipcw_L[5,1]),
          " ", " ", bonpval(ratio_il12_il4_t2_adj_ipcw_L[5,1]),
          " ", " ", bonpval(ratio_ifn_il4_t2_adj_ipcw_L[5,1]),
          " ", " ", bonpval(ratio_il12_il5_t2_adj_ipcw_L[5,1]),
          " ", " ", bonpval(ratio_ifn_il5_t2_adj_ipcw_L[5,1]),
          " ", " ", bonpval(ratio_il12_il13_t2_adj_ipcw_L[5,1]),
          " ", " ", bonpval(ratio_ifn_il13_t2_adj_ipcw_L[5,1]),
          " ", " ", bonpval(ratio_il12_il17_t2_adj_ipcw_L[5,1]),
          " ", " ", bonpval(ratio_ifn_il17_t2_adj_ipcw_L[5,1]),
          " ", " ", bonpval(ratio_il12_il21_t2_adj_ipcw_L[5,1]),
          " ", " ", bonpval(ratio_ifn_il21_t2_adj_ipcw_L[5,1]),
          " ", " ", bonpval(ratio_pro_il10_t2_adj_ipcw_L[5,1]),
          " ", " ", bonpval(ratio_th1_il10_t2_adj_ipcw_L[5,1]),
          " ", " ", bonpval(ratio_th2_il10_t2_adj_ipcw_L[5,1]),
          " ", " ", bonpval(ratio_th17_il10_t2_adj_ipcw_L[5,1]),
          " ", " ", bonpval(ratio_th1_th2_t2_adj_ipcw_L[5,1]),
          " ", " ", bonpval(ratio_th1_th17_t2_adj_ipcw_L[5,1]))

tbls3<-data.table(
  " " = tbl2[[1]],
  " " = tbl2[[2]],
  " " = tbl2[[3]],
  " " = tbl2[[4]],
  " " = tbl2[[5]],
  " " = tbl2[[6]],
  "Unadjusted difference: Intervention vs. Control" = tbl2[[7]],
  " " = rawpval2,
  " " = as_bonf_col(tbl2[[8]]),
  "Age- and sex- adjusted difference: Intervention vs. Control" = asadjtbls3,
  " " = rawpvals3,
  " " = pvals3,
  "Fully adjusted difference: Intervention vs. Control" = tbl2[[9]],
  " " = rawadjpval2,
  " " = as_bonf_col(tbl2[[10]]),
  "IPCW adjusted difference: Intervention vs. Control" = ipcws3,
  " " = rawipcwpvals3,
  " " = ipcwpvals3)

write.csv(tbls3, file=here('tables/supplementary/immune_supptable3.csv'))
print(xtable(tbls3), type="html", file=here("tables/supplementary/immune_supptable3.html"))



#### TABLE S4 ####

outcomes4 <- c("Outcome, Arm", paste("Ln IL-1", "\u03b2", " (pg/ml)", sep=""), "Control", "Nutrition + WSH", 
             "Ln IL-6 (pg/ml)", "Control", "Nutrition + WSH", 
             paste("Ln TNF-", "\u03b1", " (pg/ml)", sep=""), "Control", "Nutrition + WSH",
             "Ln CRP (mg/L)", "Control", "Nutrition + WSH", 
             "Ln IL-12 (pg/ml)", "Control", "Nutrition + WSH", 
             paste("Ln IFN-", "\u03b3", " (pg/ml)", sep=""), "Control", "Nutrition + WSH", 
             "Ln IL-4 (pg/ml)", "Control", "Nutrition + WSH", 
             "Ln IL-5 (pg/ml)", "Control", "Nutrition + WSH", 
             "Ln IL-13 (pg/ml)", "Control", "Nutrition + WSH", 
             "Ln IL-17A (pg/ml)", "Control", "Nutrition + WSH", 
             "Ln IL-21 (pg/ml)", "Control", "Nutrition + WSH", 
             "Ln IL-10 (pg/ml)", "Control", "Nutrition + WSH", 
             "Ln IL-2 (pg/ml)", "Control", "Nutrition + WSH", 
             "Ln GM-CSF (pg/ml)", "Control", "Nutrition + WSH", 
             "Ln AGP (g/L)", "Control", "Nutrition + WSH", 
             paste("Ln IGF-1 (", "\u03bc", "g/L)", sep=""), "Control", "Nutrition + WSH")

Ns4 <- c("N", " ", as.character(il1_t3_N_tr$t3_ln_il1_N_tr[1]), as.character(il1_t3_N_tr$t3_ln_il1_N_tr[2]),
       " ", as.character(il6_t3_N_tr$t3_ln_il6_N_tr[1]), as.character(il6_t3_N_tr$t3_ln_il6_N_tr[2]),
       " ", as.character(tnf_t3_N_tr$t3_ln_tnf_N_tr[1]), as.character(tnf_t3_N_tr$t3_ln_tnf_N_tr[2]),
       " ", as.character(crp_t3_N_tr$t3_ln_crp_N_tr[1]), as.character(crp_t3_N_tr$t3_ln_crp_N_tr[3]),
       " ", as.character(il12_t3_N_tr$t3_ln_il12_N_tr[1]), as.character(il12_t3_N_tr$t3_ln_il12_N_tr[2]),
       " ", as.character(ifn_t3_N_tr$t3_ln_ifn_N_tr[1]), as.character(ifn_t3_N_tr$t3_ln_ifn_N_tr[2]),
       " ", as.character(il4_t3_N_tr$t3_ln_il4_N_tr[1]), as.character(il4_t3_N_tr$t3_ln_il4_N_tr[2]),
       " ", as.character(il5_t3_N_tr$t3_ln_il5_N_tr[1]), as.character(il5_t3_N_tr$t3_ln_il5_N_tr[2]),
       " ", as.character(il13_t3_N_tr$t3_ln_il13_N_tr[1]), as.character(il13_t3_N_tr$t3_ln_il13_N_tr[2]),
       " ", as.character(il17_t3_N_tr$t3_ln_il17_N_tr[1]), as.character(il17_t3_N_tr$t3_ln_il17_N_tr[2]),
       " ", as.character(il21_t3_N_tr$t3_ln_il21_N_tr[1]), as.character(il21_t3_N_tr$t3_ln_il21_N_tr[2]),
       " ", as.character(il10_t3_N_tr$t3_ln_il10_N_tr[1]), as.character(il10_t3_N_tr$t3_ln_il10_N_tr[2]),
       " ", as.character(il2_t3_N_tr$t3_ln_il2_N_tr[1]), as.character(il2_t3_N_tr$t3_ln_il2_N_tr[2]),
       " ", as.character(gmc_t3_N_tr$t3_ln_gmc_N_tr[1]), as.character(gmc_t3_N_tr$t3_ln_gmc_N_tr[2]),
       " ", as.character(agp_t3_N_tr$t3_ln_agp_N_tr[1]), as.character(agp_t3_N_tr$t3_ln_agp_N_tr[2]),
       " ", as.character(igf_t3_N_tr$t3_ln_igf_N_tr[1]), as.character(igf_t3_N_tr$t3_ln_igf_N_tr[2]))

absmeans4 <- c("Absolute Mean", " ", as.character(round(abs_il1_t3_N_tr$mean[1], 2)), as.character(round(abs_il1_t3_N_tr$mean[2], 2)),
             " ", as.character(round(abs_il6_t3_N_tr$mean[1], 2)), as.character(round(abs_il6_t3_N_tr$mean[2], 2)),
             " ", as.character(round(abs_tnf_t3_N_tr$mean[1], 2)), as.character(round(abs_tnf_t3_N_tr$mean[2], 2)),
             " ", as.character(round(abs_crp_t3_N_tr$mean[1], 2)), as.character(round(abs_crp_t3_N_tr$mean[3], 2)),
             " ", as.character(round(abs_il12_t3_N_tr$mean[1], 2)), as.character(round(abs_il12_t3_N_tr$mean[2], 2)),
             " ", as.character(round(abs_ifn_t3_N_tr$mean[1], 2)), as.character(round(abs_ifn_t3_N_tr$mean[2], 2)),
             " ", as.character(round(abs_il4_t3_N_tr$mean[1], 2)), as.character(round(abs_il4_t3_N_tr$mean[2], 2)),
             " ", as.character(round(abs_il5_t3_N_tr$mean[1], 2)), as.character(round(abs_il5_t3_N_tr$mean[2], 2)),
             " ", as.character(round(abs_il13_t3_N_tr$mean[1], 2)), as.character(round(abs_il13_t3_N_tr$mean[2], 2)),
             " ", as.character(round(abs_il17_t3_N_tr$mean[1], 2)), as.character(round(abs_il17_t3_N_tr$mean[2], 2)),
             " ", as.character(round(abs_il21_t3_N_tr$mean[1], 2)), as.character(round(abs_il21_t3_N_tr$mean[2], 2)),
             " ", as.character(round(abs_il10_t3_N_tr$mean[1], 2)), as.character(round(abs_il10_t3_N_tr$mean[2], 2)),
             " ", as.character(round(abs_il2_t3_N_tr$mean[1], 2)), as.character(round(abs_il2_t3_N_tr$mean[2], 2)),
             " ", as.character(round(abs_gmc_t3_N_tr$mean[1], 2)), as.character(round(abs_gmc_t3_N_tr$mean[2], 2)),
             " ", as.character(round(abs_agp_t3_N_tr$mean[1], 2)), as.character(round(abs_agp_t3_N_tr$mean[3], 2)),
             " ", as.character(round(abs_igf_t3_N_tr$mean[1], 2)), as.character(round(abs_igf_t3_N_tr$mean[2], 2)))

abssds4 <- c("Absolute SD", " ", as.character(round(abs_il1_t3_N_tr$sd[1], 2)), as.character(round(abs_il1_t3_N_tr$sd[2], 2)),
           " ", as.character(round(abs_il6_t3_N_tr$sd[1], 2)), as.character(round(abs_il6_t3_N_tr$sd[2], 2)),
           " ", as.character(round(abs_tnf_t3_N_tr$sd[1], 2)), as.character(round(abs_tnf_t3_N_tr$sd[2], 2)),
           " ", as.character(round(abs_crp_t3_N_tr$sd[1], 2)), as.character(round(abs_crp_t3_N_tr$sd[3], 2)),
           " ", as.character(round(abs_il12_t3_N_tr$sd[1], 2)), as.character(round(abs_il12_t3_N_tr$sd[2], 2)),
           " ", as.character(round(abs_ifn_t3_N_tr$sd[1], 2)), as.character(round(abs_ifn_t3_N_tr$sd[2], 2)),
           " ", as.character(round(abs_il4_t3_N_tr$sd[1], 2)), as.character(round(abs_il4_t3_N_tr$sd[2], 2)),
           " ", as.character(round(abs_il5_t3_N_tr$sd[1], 2)), as.character(round(abs_il5_t3_N_tr$sd[2], 2)),
           " ", as.character(round(abs_il13_t3_N_tr$sd[1], 2)), as.character(round(abs_il13_t3_N_tr$sd[2], 2)),
           " ", as.character(round(abs_il17_t3_N_tr$sd[1], 2)), as.character(round(abs_il17_t3_N_tr$sd[2], 2)),
           " ", as.character(round(abs_il21_t3_N_tr$sd[1], 2)), as.character(round(abs_il21_t3_N_tr$sd[2], 2)),
           " ", as.character(round(abs_il10_t3_N_tr$sd[1], 2)), as.character(round(abs_il10_t3_N_tr$sd[2], 2)),
           " ", as.character(round(abs_il2_t3_N_tr$sd[1], 2)), as.character(round(abs_il2_t3_N_tr$sd[2], 2)),
           " ", as.character(round(abs_gmc_t3_N_tr$sd[1], 2)), as.character(round(abs_gmc_t3_N_tr$sd[2], 2)),
           " ", as.character(round(abs_agp_t3_N_tr$sd[1], 2)), as.character(round(abs_agp_t3_N_tr$sd[3], 2)),
           " ", as.character(round(abs_igf_t3_N_tr$sd[1], 2)), as.character(round(abs_igf_t3_N_tr$sd[2], 2)))

means4 <- c("Mean", " ", as.character(round(il1_t3_N_tr$mean[1], 2)), as.character(round(il1_t3_N_tr$mean[2], 2)),
          " ", as.character(round(il6_t3_N_tr$mean[1], 2)), as.character(round(il6_t3_N_tr$mean[2], 2)),
          " ", as.character(round(tnf_t3_N_tr$mean[1], 2)), as.character(round(tnf_t3_N_tr$mean[2], 2)),
          " ", as.character(round(crp_t3_N_tr$mean[1], 2)), as.character(round(crp_t3_N_tr$mean[3], 2)),
          " ", as.character(round(il12_t3_N_tr$mean[1], 2)), as.character(round(il12_t3_N_tr$mean[2], 2)),
          " ", as.character(round(ifn_t3_N_tr$mean[1], 2)), as.character(round(ifn_t3_N_tr$mean[2], 2)),
          " ", as.character(round(il4_t3_N_tr$mean[1], 2)), as.character(round(il4_t3_N_tr$mean[2], 2)),
          " ", as.character(round(il5_t3_N_tr$mean[1], 2)), as.character(round(il5_t3_N_tr$mean[2], 2)),
          " ", as.character(round(il13_t3_N_tr$mean[1], 2)), as.character(round(il13_t3_N_tr$mean[2], 2)),
          " ", as.character(round(il17_t3_N_tr$mean[1], 2)), as.character(round(il17_t3_N_tr$mean[2], 2)),
          " ", as.character(round(il21_t3_N_tr$mean[1], 2)), as.character(round(il21_t3_N_tr$mean[2], 2)),
          " ", as.character(round(il10_t3_N_tr$mean[1], 2)), as.character(round(il10_t3_N_tr$mean[2], 2)),
          " ", as.character(round(il2_t3_N_tr$mean[1], 2)), as.character(round(il2_t3_N_tr$mean[2], 2)),
          " ", as.character(round(gmc_t3_N_tr$mean[1], 2)), as.character(round(gmc_t3_N_tr$mean[2], 2)),
          " ", as.character(round(agp_t3_N_tr$mean[1], 2)), as.character(round(agp_t3_N_tr$mean[3], 2)),
          " ", as.character(round(igf_t3_N_tr$mean[1], 2)), as.character(round(igf_t3_N_tr$mean[2], 2)))

sds4 <- c("SD", " ", as.character(round(il1_t3_N_tr$sd[1], 2)), as.character(round(il1_t3_N_tr$sd[2], 2)),
        " ", as.character(round(il6_t3_N_tr$sd[1], 2)), as.character(round(il6_t3_N_tr$sd[2], 2)),
        " ", as.character(round(tnf_t3_N_tr$sd[1], 2)), as.character(round(tnf_t3_N_tr$sd[2], 2)),
        " ", as.character(round(crp_t3_N_tr$sd[1], 2)), as.character(round(crp_t3_N_tr$sd[3], 2)),
        " ", as.character(round(il12_t3_N_tr$sd[1], 2)), as.character(round(il12_t3_N_tr$sd[2], 2)),
        " ", as.character(round(ifn_t3_N_tr$sd[1], 2)), as.character(round(ifn_t3_N_tr$sd[2], 2)),
        " ", as.character(round(il4_t3_N_tr$sd[1], 2)), as.character(round(il4_t3_N_tr$sd[2], 2)),
        " ", as.character(round(il5_t3_N_tr$sd[1], 2)), as.character(round(il5_t3_N_tr$sd[2], 2)),
        " ", as.character(round(il13_t3_N_tr$sd[1], 2)), as.character(round(il13_t3_N_tr$sd[2], 2)),
        " ", as.character(round(il17_t3_N_tr$sd[1], 2)), as.character(round(il17_t3_N_tr$sd[2], 2)),
        " ", as.character(round(il21_t3_N_tr$sd[1], 2)), as.character(round(il21_t3_N_tr$sd[2], 2)),
        " ", as.character(round(il10_t3_N_tr$sd[1], 2)), as.character(round(il10_t3_N_tr$sd[2], 2)),
        " ", as.character(round(il2_t3_N_tr$sd[1], 2)), as.character(round(il2_t3_N_tr$sd[2], 2)),
        " ", as.character(round(gmc_t3_N_tr$sd[1], 2)), as.character(round(gmc_t3_N_tr$sd[2], 2)),
        " ", as.character(round(agp_t3_N_tr$sd[1], 2)), as.character(round(agp_t3_N_tr$sd[3], 2)),
        " ", as.character(round(igf_t3_N_tr$sd[1], 2)), as.character(round(igf_t3_N_tr$sd[2], 2)))

rawpvals4 <- c("P-value", " ", " ", rawpval(t3_il1_unadj_L[6]), " ", " ", rawpval(t3_il6_unadj_L[6]), " ", " ", rawpval(t3_tnf_unadj_L[6]),
            " ", " ", rawpval(t3_crp_unadj_L[6]), " ", " ", rawpval(t3_il12_unadj_L[6]), " ", " ", rawpval(t3_ifn_unadj_L[6]),
            " ", " ", rawpval(t3_il4_unadj_L[6]), " ", " ", rawpval(t3_il5_unadj_L[6]), " ", " ", rawpval(t3_il13_unadj_L[6]),
            " ", " ", rawpval(t3_il17_unadj_L[6]), " ", " ", rawpval(t3_il21_unadj_L[6]), " ", " ", rawpval(t3_il10_unadj_L[6]),
            " ", " ", rawpval(t3_il2_unadj_L[6]), " ", " ", rawpval(t3_gmc_unadj_L[6]),
            " ", " ", rawpval(t3_agp_unadj_L[6]), " ", " ", rawpval(t3_igf_unadj_L[6]))

pvals4 <- c("Bonferroni P-value", " ", " ", bonpval(t3_il1_unadj_L[6]), " ", " ", bonpval(t3_il6_unadj_L[6]), " ", " ", bonpval(t3_tnf_unadj_L[6]),
            " ", " ", bonpval(t3_crp_unadj_L[6]), " ", " ", bonpval(t3_il12_unadj_L[6]), " ", " ", bonpval(t3_ifn_unadj_L[6]),
            " ", " ", bonpval(t3_il4_unadj_L[6]), " ", " ", bonpval(t3_il5_unadj_L[6]), " ", " ", bonpval(t3_il13_unadj_L[6]),
            " ", " ", bonpval(t3_il17_unadj_L[6]), " ", " ", bonpval(t3_il21_unadj_L[6]), " ", " ", bonpval(t3_il10_unadj_L[6]),
            " ", " ", bonpval(t3_il2_unadj_L[6]), " ", " ", bonpval(t3_gmc_unadj_L[6]),
            " ", " ", bonpval(t3_agp_unadj_L[6]), " ", " ", bonpval(t3_igf_unadj_L[6]))

t3_agp_unadj_L <- round(t3_agp_unadj_L, 2)
t3_crp_unadj_L <- round(t3_crp_unadj_L, 2)
t3_gmc_unadj_L <- round(t3_gmc_unadj_L, 2)
t3_ifn_unadj_L <- round(t3_ifn_unadj_L, 2)
t3_igf_unadj_L <- round(t3_igf_unadj_L, 2)
t3_il1_unadj_L <- round(t3_il1_unadj_L, 2)
t3_il10_unadj_L <- round(t3_il10_unadj_L, 2)
t3_il12_unadj_L <- round(t3_il12_unadj_L, 2)
t3_il13_unadj_L <- round(t3_il13_unadj_L, 2)
t3_il17_unadj_L <- round(t3_il17_unadj_L, 2)
t3_il2_unadj_L <- round(t3_il2_unadj_L, 2)
t3_il21_unadj_L <- round(t3_il21_unadj_L, 2)
t3_il4_unadj_L <- round(t3_il4_unadj_L, 2)
t3_il5_unadj_L <- round(t3_il5_unadj_L, 2)
t3_il6_unadj_L <- round(t3_il6_unadj_L, 2)
t3_tnf_unadj_L <- round(t3_tnf_unadj_L, 2)

unadjs4 <- c("95% CI", " ", " ", paste(t3_il1_unadj_L[1], " (", t3_il1_unadj_L[2], ", ", t3_il1_unadj_L[3], ")", sep=""),
           " ", " ", paste(t3_il6_unadj_L[1], " (", t3_il6_unadj_L[2], ", ", t3_il6_unadj_L[3], ")", sep=""),
           " ", " ", paste(t3_tnf_unadj_L[1], " (", t3_tnf_unadj_L[2], ", ", t3_tnf_unadj_L[3], ")", sep=""),
           " ", " ", paste(t3_crp_unadj_L[1], " (", t3_crp_unadj_L[2], ", ", t3_crp_unadj_L[3], ")", sep=""),
           " ", " ", paste(t3_il12_unadj_L[1], " (", t3_il12_unadj_L[2], ", ", t3_il12_unadj_L[3], ")", sep=""),
           " ", " ", paste(t3_ifn_unadj_L[1], " (", t3_ifn_unadj_L[2], ", ", t3_ifn_unadj_L[3], ")", sep=""),
           " ", " ", paste(t3_il4_unadj_L[1], " (", t3_il4_unadj_L[2], ", ", t3_il4_unadj_L[3], ")", sep=""),
           " ", " ", paste(t3_il5_unadj_L[1], " (", t3_il5_unadj_L[2], ", ", t3_il5_unadj_L[3], ")", sep=""),
           " ", " ", paste(t3_il13_unadj_L[1], " (", t3_il13_unadj_L[2], ", ", t3_il13_unadj_L[3], ")", sep=""),
           " ", " ", paste(t3_il17_unadj_L[1], " (", t3_il17_unadj_L[2], ", ", t3_il17_unadj_L[3], ")", sep=""),
           " ", " ", paste(t3_il21_unadj_L[1], " (", t3_il21_unadj_L[2], ", ", t3_il21_unadj_L[3], ")", sep=""),
           " ", " ", paste(t3_il10_unadj_L[1], " (", t3_il10_unadj_L[2], ", ", t3_il10_unadj_L[3], ")", sep=""),
           " ", " ", paste(t3_il2_unadj_L[1], " (", t3_il2_unadj_L[2], ", ", t3_il2_unadj_L[3], ")", sep=""),
           " ", " ", paste(t3_gmc_unadj_L[1], " (", t3_gmc_unadj_L[2], ", ", t3_gmc_unadj_L[3], ")", sep=""),
           " ", " ", paste(t3_agp_unadj_L[1], " (", t3_agp_unadj_L[2], ", ", t3_agp_unadj_L[3], ")", sep=""),
           " ", " ", paste(t3_igf_unadj_L[1], " (", t3_igf_unadj_L[2], ", ", t3_igf_unadj_L[3], ")", sep=""))

rawasapvals4 <- c("P-value", " ", " ", rawpval(t3_il1_adj_sex_age_L[6]), " ", " ", rawpval(t3_il6_adj_sex_age_L[6]), " ", " ", rawpval(t3_tnf_adj_sex_age_L[6]),
               " ", " ", rawpval(t3_crp_adj_sex_age_L[6]), " ", " ", rawpval(t3_il12_adj_sex_age_L[6]), " ", " ", rawpval(t3_ifn_adj_sex_age_L[6]),
               " ", " ", rawpval(t3_il4_adj_sex_age_L[6]), " ", " ", rawpval(t3_il5_adj_sex_age_L[6]), " ", " ", rawpval(t3_il13_adj_sex_age_L[6]),
               " ", " ", rawpval(t3_il17_adj_sex_age_L[6]), " ", " ", rawpval(t3_il21_adj_sex_age_L[6]), " ", " ", rawpval(t3_il10_adj_sex_age_L[6]),
               " ", " ", rawpval(t3_il2_adj_sex_age_L[6]), " ", " ", rawpval(t3_gmc_adj_sex_age_L[6]),
               " ", " ", rawpval(t3_agp_adj_sex_age_L[6]), " ", " ", rawpval(t3_igf_adj_sex_age_L[6]))

asapvals4 <- c("Bonferroni P-value", " ", " ", bonpval(t3_il1_adj_sex_age_L[6]), " ", " ", bonpval(t3_il6_adj_sex_age_L[6]), " ", " ", bonpval(t3_tnf_adj_sex_age_L[6]),
               " ", " ", bonpval(t3_crp_adj_sex_age_L[6]), " ", " ", bonpval(t3_il12_adj_sex_age_L[6]), " ", " ", bonpval(t3_ifn_adj_sex_age_L[6]),
               " ", " ", bonpval(t3_il4_adj_sex_age_L[6]), " ", " ", bonpval(t3_il5_adj_sex_age_L[6]), " ", " ", bonpval(t3_il13_adj_sex_age_L[6]),
               " ", " ", bonpval(t3_il17_adj_sex_age_L[6]), " ", " ", bonpval(t3_il21_adj_sex_age_L[6]), " ", " ", bonpval(t3_il10_adj_sex_age_L[6]),
               " ", " ", bonpval(t3_il2_adj_sex_age_L[6]), " ", " ", bonpval(t3_gmc_adj_sex_age_L[6]),
               " ", " ", bonpval(t3_agp_adj_sex_age_L[6]), " ", " ", bonpval(t3_igf_adj_sex_age_L[6]))

t3_crp_adj_sex_age_L <- round(t3_crp_adj_sex_age_L, 2)
t3_agp_adj_sex_age_L <- round(t3_agp_adj_sex_age_L, 2)
t3_gmc_adj_sex_age_L <- round(t3_gmc_adj_sex_age_L, 2)
t3_ifn_adj_sex_age_L <- round(t3_ifn_adj_sex_age_L, 2)
t3_igf_adj_sex_age_L <- round(t3_igf_adj_sex_age_L, 2)
t3_il1_adj_sex_age_L <- round(t3_il1_adj_sex_age_L, 2)
t3_il10_adj_sex_age_L <- round(t3_il10_adj_sex_age_L, 2)
t3_il12_adj_sex_age_L <- round(t3_il12_adj_sex_age_L, 2)
t3_il13_adj_sex_age_L <- round(t3_il13_adj_sex_age_L, 2)
t3_il17_adj_sex_age_L <- round(t3_il17_adj_sex_age_L, 2)
t3_il2_adj_sex_age_L <- round(t3_il2_adj_sex_age_L, 2)
t3_il21_adj_sex_age_L <- round(t3_il21_adj_sex_age_L, 2)
t3_il4_adj_sex_age_L <- round(t3_il4_adj_sex_age_L, 2)
t3_il5_adj_sex_age_L <- round(t3_il5_adj_sex_age_L, 2)
t3_il6_adj_sex_age_L <- round(t3_il6_adj_sex_age_L, 2)
t3_tnf_adj_sex_age_L <- round(t3_tnf_adj_sex_age_L, 2)

asadjs4 <- c("95% CI", " ", " ", paste(t3_il1_adj_sex_age_L[1], " (", t3_il1_adj_sex_age_L[2], ", ", t3_il1_adj_sex_age_L[3], ")", sep=""),
           " ", " ", paste(t3_il6_adj_sex_age_L[1], " (", t3_il6_adj_sex_age_L[2], ", ", t3_il6_adj_sex_age_L[3], ")", sep=""),
           " ", " ", paste(t3_tnf_adj_sex_age_L[1], " (", t3_tnf_adj_sex_age_L[2], ", ", t3_tnf_adj_sex_age_L[3], ")", sep=""),
           " ", " ", paste(t3_crp_adj_sex_age_L[1], " (", t3_crp_adj_sex_age_L[2], ", ", t3_crp_adj_sex_age_L[3], ")", sep=""),
           " ", " ", paste(t3_il12_adj_sex_age_L[1], " (", t3_il12_adj_sex_age_L[2], ", ", t3_il12_adj_sex_age_L[3], ")", sep=""),
           " ", " ", paste(t3_ifn_adj_sex_age_L[1], " (", t3_ifn_adj_sex_age_L[2], ", ", t3_ifn_adj_sex_age_L[3], ")", sep=""),
           " ", " ", paste(t3_il4_adj_sex_age_L[1], " (", t3_il4_adj_sex_age_L[2], ", ", t3_il4_adj_sex_age_L[3], ")", sep=""),
           " ", " ", paste(t3_il5_adj_sex_age_L[1], " (", t3_il5_adj_sex_age_L[2], ", ", t3_il5_adj_sex_age_L[3], ")", sep=""),
           " ", " ", paste(t3_il13_adj_sex_age_L[1], " (", t3_il13_adj_sex_age_L[2], ", ", t3_il13_adj_sex_age_L[3], ")", sep=""),
           " ", " ", paste(t3_il17_adj_sex_age_L[1], " (", t3_il17_adj_sex_age_L[2], ", ", t3_il17_adj_sex_age_L[3], ")", sep=""),
           " ", " ", paste(t3_il21_adj_sex_age_L[1], " (", t3_il21_adj_sex_age_L[2], ", ", t3_il21_adj_sex_age_L[3], ")", sep=""),
           " ", " ", paste(t3_il10_adj_sex_age_L[1], " (", t3_il10_adj_sex_age_L[2], ", ", t3_il10_adj_sex_age_L[3], ")", sep=""),
           " ", " ", paste(t3_il2_adj_sex_age_L[1], " (", t3_il2_adj_sex_age_L[2], ", ", t3_il2_adj_sex_age_L[3], ")", sep=""),
           " ", " ", paste(t3_gmc_adj_sex_age_L[1], " (", t3_gmc_adj_sex_age_L[2], ", ", t3_gmc_adj_sex_age_L[3], ")", sep=""),
           " ", " ", paste(t3_agp_adj_sex_age_L[1], " (", t3_agp_adj_sex_age_L[2], ", ", t3_agp_adj_sex_age_L[3], ")", sep=""),
           " ", " ", paste(t3_igf_adj_sex_age_L[1], " (", t3_igf_adj_sex_age_L[2], ", ", t3_igf_adj_sex_age_L[3], ")", sep=""))

rawadjpvals4 <- c("P-value", " ", " ", rawpval(t3_il1_adj_L[6]), " ", " ", rawpval(t3_il6_adj_L[6]), " ", " ", rawpval(t3_tnf_adj_L[6]),
               " ", " ", rawpval(t3_crp_adj_L[6]), " ", " ", rawpval(t3_il12_adj_L[6]), " ", " ", rawpval(t3_ifn_adj_L[6]),
               " ", " ", rawpval(t3_il4_adj_L[6]), " ", " ", rawpval(t3_il5_adj_L[6]), " ", " ", rawpval(t3_il13_adj_L[6]),
               " ", " ", rawpval(t3_il17_adj_L[6]), " ", " ", rawpval(t3_il21_adj_L[6]), " ", " ", rawpval(t3_il10_adj_L[6]),
               " ", " ", rawpval(t3_il2_adj_L[6]), " ", " ", rawpval(t3_gmc_adj_L[6]),
               " ", " ", rawpval(t3_agp_adj_L[6]), " ", " ", rawpval(t3_igf_adj_L[6]))

adjpvals4 <- c("Bonferroni P-value", " ", " ", bonpval(t3_il1_adj_L[6]), " ", " ", bonpval(t3_il6_adj_L[6]), " ", " ", bonpval(t3_tnf_adj_L[6]),
               " ", " ", bonpval(t3_crp_adj_L[6]), " ", " ", bonpval(t3_il12_adj_L[6]), " ", " ", bonpval(t3_ifn_adj_L[6]),
               " ", " ", bonpval(t3_il4_adj_L[6]), " ", " ", bonpval(t3_il5_adj_L[6]), " ", " ", bonpval(t3_il13_adj_L[6]),
               " ", " ", bonpval(t3_il17_adj_L[6]), " ", " ", bonpval(t3_il21_adj_L[6]), " ", " ", bonpval(t3_il10_adj_L[6]),
               " ", " ", bonpval(t3_il2_adj_L[6]), " ", " ", bonpval(t3_gmc_adj_L[6]),
               " ", " ", bonpval(t3_agp_adj_L[6]), " ", " ", bonpval(t3_igf_adj_L[6]))

t3_crp_adj_L <- round(t3_crp_adj_L, 2)
t3_agp_adj_L <- round(t3_agp_adj_L, 2)
t3_gmc_adj_L <- round(t3_gmc_adj_L, 2)
t3_ifn_adj_L <- round(t3_ifn_adj_L, 2)
t3_igf_adj_L <- round(t3_igf_adj_L, 2)
t3_il1_adj_L <- round(t3_il1_adj_L, 2)
t3_il10_adj_L <- round(t3_il10_adj_L, 2)
t3_il12_adj_L <- round(t3_il12_adj_L, 2)
t3_il13_adj_L <- round(t3_il13_adj_L, 2)
t3_il17_adj_L <- round(t3_il17_adj_L, 2)
t3_il2_adj_L <- round(t3_il2_adj_L, 2)
t3_il21_adj_L <- round(t3_il21_adj_L, 2)
t3_il4_adj_L <- round(t3_il4_adj_L, 2)
t3_il5_adj_L <- round(t3_il5_adj_L, 2)
t3_il6_adj_L <- round(t3_il6_adj_L, 2)
t3_tnf_adj_L <- round(t3_tnf_adj_L, 2)

adjs4 <- c("95% CI", " ", " ", paste(t3_il1_adj_L[1], " (", t3_il1_adj_L[2], ", ", t3_il1_adj_L[3], ")", sep=""),
         " ", " ", paste(t3_il6_adj_L[1], " (", t3_il6_adj_L[2], ", ", t3_il6_adj_L[3], ")", sep=""),
         " ", " ", paste(t3_tnf_adj_L[1], " (", t3_tnf_adj_L[2], ", ", t3_tnf_adj_L[3], ")", sep=""),
         " ", " ", paste(t3_crp_adj_L[1], " (", t3_crp_adj_L[2], ", ", t3_crp_adj_L[3], ")", sep=""),
         " ", " ", paste(t3_il12_adj_L[1], " (", t3_il12_adj_L[2], ", ", t3_il12_adj_L[3], ")", sep=""),
         " ", " ", paste(t3_ifn_adj_L[1], " (", t3_ifn_adj_L[2], ", ", t3_ifn_adj_L[3], ")", sep=""),
         " ", " ", paste(t3_il4_adj_L[1], " (", t3_il4_adj_L[2], ", ", t3_il4_adj_L[3], ")", sep=""),
         " ", " ", paste(t3_il5_adj_L[1], " (", t3_il5_adj_L[2], ", ", t3_il5_adj_L[3], ")", sep=""),
         " ", " ", paste(t3_il13_adj_L[1], " (", t3_il13_adj_L[2], ", ", t3_il13_adj_L[3], ")", sep=""),
         " ", " ", paste(t3_il17_adj_L[1], " (", t3_il17_adj_L[2], ", ", t3_il17_adj_L[3], ")", sep=""),
         " ", " ", paste(t3_il21_adj_L[1], " (", t3_il21_adj_L[2], ", ", t3_il21_adj_L[3], ")", sep=""),
         " ", " ", paste(t3_il10_adj_L[1], " (", t3_il10_adj_L[2], ", ", t3_il10_adj_L[3], ")", sep=""),
         " ", " ", paste(t3_il2_adj_L[1], " (", t3_il2_adj_L[2], ", ", t3_il2_adj_L[3], ")", sep=""),
         " ", " ", paste(t3_gmc_adj_L[1], " (", t3_gmc_adj_L[2], ", ", t3_gmc_adj_L[3], ")", sep=""),
         " ", " ", paste(t3_agp_adj_L[1], " (", t3_agp_adj_L[2], ", ", t3_agp_adj_L[3], ")", sep=""),
         " ", " ", paste(t3_igf_adj_L[1], " (", t3_igf_adj_L[2], ", ", t3_igf_adj_L[3], ")", sep=""))

ipcws4<-c("95% CI", " ", " ", maketblvalue(il1_t3_adj_ipcw_L$`unlist(il1_t3_adj_ipcw$estimates$ATE)`),
          " ", " ", maketblvalue(il6_t3_adj_ipcw_L$`unlist(il6_t3_adj_ipcw$estimates$ATE)`),
          " ", " ", maketblvalue(tnf_t3_adj_ipcw_L$`unlist(tnf_t3_adj_ipcw$estimates$ATE)`),
          " ", " ", maketblvalue(crp_t3_ipcw_L$`unlist(crp_t3_ipcw_L$estimates$ATE)`),
          " ", " ", maketblvalue(il12_t3_adj_ipcw_L$`unlist(il12_t3_adj_ipcw$estimates$ATE)`),
          " ", " ", maketblvalue(ifn_t3_adj_ipcw_L$`unlist(ifn_t3_adj_ipcw$estimates$ATE)`),
          " ", " ", maketblvalue(il4_t3_adj_ipcw_L$`unlist(il4_t3_adj_ipcw$estimates$ATE)`),
          " ", " ", maketblvalue(il5_t3_adj_ipcw_L$`unlist(il5_t3_adj_ipcw$estimates$ATE)`),
          " ", " ", maketblvalue(il13_t3_adj_ipcw_L$`unlist(il13_t3_adj_ipcw$estimates$ATE)`),
          " ", " ", maketblvalue(il17_t3_adj_ipcw_L$`unlist(il17_t3_adj_ipcw$estimates$ATE)`),
          " ", " ", maketblvalue(il21_t3_adj_ipcw_L$`unlist(il21_t3_adj_ipcw$estimates$ATE)`),
          " ", " ", maketblvalue(il10_t3_adj_ipcw_L$`unlist(il10_t3_adj_ipcw$estimates$ATE)`),
          " ", " ", maketblvalue(il2_t3_adj_ipcw_L$`unlist(il2_t3_adj_ipcw$estimates$ATE)`),
          " ", " ", maketblvalue(gmc_t3_adj_ipcw_L$`unlist(gmc_t3_adj_ipcw$estimates$ATE)`), 
          " ", " ", maketblvalue(agp_t3_ipcw_L$`unlist(agp_t3_ipcw_L$estimates$ATE)`), 
          " ", " ", maketblvalue(igf_t3_adj_ipcw_L$`unlist(igf_t3_adj_ipcw$estimates$ATE)`))

rawipcwpvals4<-c("P-value", " ", " ", rawpval(il1_t3_adj_ipcw_L[5,1]),
          " ", " ", rawpval(il6_t3_adj_ipcw_L[5,1]),
          " ", " ", rawpval(tnf_t3_adj_ipcw_L[5,1]),
          " ", " ", rawpval(crp_t3_ipcw_L[5,1]),
          " ", " ", rawpval(il12_t3_adj_ipcw_L[5,1]),
          " ", " ", rawpval(ifn_t3_adj_ipcw_L[5,1]),
          " ", " ", rawpval(il4_t3_adj_ipcw_L[5,1]),
          " ", " ", rawpval(il5_t3_adj_ipcw_L[5,1]),
          " ", " ", rawpval(il13_t3_adj_ipcw_L[5,1]),
          " ", " ", rawpval(il17_t3_adj_ipcw_L[5,1]),
          " ", " ", rawpval(il21_t3_adj_ipcw_L[5,1]),
          " ", " ", rawpval(il10_t3_adj_ipcw_L[5,1]),
          " ", " ", rawpval(il2_t3_adj_ipcw_L[5,1]),
          " ", " ", rawpval(gmc_t3_adj_ipcw_L[5,1]), 
          " ", " ", rawpval(agp_t3_ipcw_L[5,1]), 
          " ", " ", rawpval(igf_t3_adj_ipcw_L[5,1]))

ipcwpvals4<-c("Bonferroni P-value", " ", " ", bonpval(il1_t3_adj_ipcw_L[5,1]),
          " ", " ", bonpval(il6_t3_adj_ipcw_L[5,1]),
          " ", " ", bonpval(tnf_t3_adj_ipcw_L[5,1]),
          " ", " ", bonpval(crp_t3_ipcw_L[5,1]),
          " ", " ", bonpval(il12_t3_adj_ipcw_L[5,1]),
          " ", " ", bonpval(ifn_t3_adj_ipcw_L[5,1]),
          " ", " ", bonpval(il4_t3_adj_ipcw_L[5,1]),
          " ", " ", bonpval(il5_t3_adj_ipcw_L[5,1]),
          " ", " ", bonpval(il13_t3_adj_ipcw_L[5,1]),
          " ", " ", bonpval(il17_t3_adj_ipcw_L[5,1]),
          " ", " ", bonpval(il21_t3_adj_ipcw_L[5,1]),
          " ", " ", bonpval(il10_t3_adj_ipcw_L[5,1]),
          " ", " ", bonpval(il2_t3_adj_ipcw_L[5,1]),
          " ", " ", bonpval(gmc_t3_adj_ipcw_L[5,1]), 
          " ", " ", bonpval(agp_t3_ipcw_L[5,1]), 
          " ", " ", bonpval(igf_t3_adj_ipcw_L[5,1]))

# Table S4: Effect of intervention on individual immune status and growth factor measurements at age 28 months
tbls4 <- data.table(
  " " = outcomes4,
  " " = Ns4, 
  " " = absmeans4,
  " " = abssds4,
  " " = means4, 
  " " = sds4,
  "Unadjusted difference: Intervention vs. Control" = unadjs4,
  " " = rawpvals4,
  " " = pvals4,
  "Age- and sex- adjusted difference: Intervention vs. Control" = asadjs4, 
  " " = rawasapvals4,
  " " = asapvals4,
  "Fully adjusted difference: Intervention vs. Control" = adjs4,
  " " = rawadjpvals4,
  " " = adjpvals4, 
  "IPCW adjusted difference: Intervention vs. Control" = ipcws4,
  " " = rawipcwpvals4,
  " " = ipcwpvals4
)

write.csv(tbls4, file=here('tables/supplementary/immune_supptable4.csv'))
print(xtable(tbls4), type="html", file=here("tables/supplementary/immune_supptable4.html"))


#### TABLE S5 ####

rawpvals5 <- c("P-value", " ", " ", rawpval(t3_ratio_il1_il10_adj_sex_age_L[6]), " ", " ", rawpval(t3_ratio_il6_il10_adj_sex_age_L[6]), " ", " ", rawpval(t3_ratio_tnf_il10_adj_sex_age_L[6]), 
            " ", " ", rawpval(t3_ratio_il12_il10_adj_sex_age_L[6]), " ", " ", rawpval(t3_ratio_ifn_il10_adj_sex_age_L[6]), " ", " ", rawpval(t3_ratio_il4_il10_adj_sex_age_L[6]), 
            " ", " ", rawpval(t3_ratio_il5_il10_adj_sex_age_L[6]), " ", " ", rawpval(t3_ratio_il13_il10_adj_sex_age_L[6]), " ", " ", rawpval(t3_ratio_il17_il10_adj_sex_age_L[6]), 
            " ", " ", rawpval(t3_ratio_il21_il10_adj_sex_age_L[6]), " ", " ", rawpval(t3_ratio_il2_il10_adj_sex_age_L[6]), " ", " ", rawpval(t3_ratio_gmc_il10_adj_sex_age_L[6]), 
            " ", " ", rawpval(t3_ratio_il12_il4_adj_sex_age_L[6]), " ", " ", rawpval(t3_ratio_ifn_il4_adj_sex_age_L[6]), " ", " ", rawpval(t3_ratio_il12_il5_adj_sex_age_L[6]), 
            " ", " ", rawpval(t3_ratio_ifn_il5_adj_sex_age_L[6]), " ", " ", rawpval(t3_ratio_il12_il13_adj_sex_age_L[6]), " ", " ", rawpval(t3_ratio_ifn_il13_adj_sex_age_L[6]), 
            " ", " ", rawpval(t3_ratio_il12_il17_adj_sex_age_L[6]), " ", " ", rawpval(t3_ratio_ifn_il17_adj_sex_age_L[6]), " ", " ", rawpval(t3_ratio_il12_il21_adj_sex_age_L[6]), 
            " ", " ", rawpval(t3_ratio_ifn_il21_adj_sex_age_L[6]), " ", " ", rawpval(t3_ratio_pro_il10_adj_sex_age_L[6]), " ", " ", rawpval(t3_ratio_th1_il10_adj_sex_age_L[6]), 
            " ", " ", rawpval(t3_ratio_th2_il10_adj_sex_age_L[6]), " ", " ", rawpval(t3_ratio_th17_il10_adj_sex_age_L[6]), " ", " ", rawpval(t3_ratio_th1_th2_adj_sex_age_L[6]), " ", " ", rawpval(t3_ratio_th1_th17_adj_sex_age_L[6])
)

pvals5 <- c("Bonferroni P-value", " ", " ", bonpval(t3_ratio_il1_il10_adj_sex_age_L[6]), " ", " ", bonpval(t3_ratio_il6_il10_adj_sex_age_L[6]), " ", " ", bonpval(t3_ratio_tnf_il10_adj_sex_age_L[6]), 
            " ", " ", bonpval(t3_ratio_il12_il10_adj_sex_age_L[6]), " ", " ", bonpval(t3_ratio_ifn_il10_adj_sex_age_L[6]), " ", " ", bonpval(t3_ratio_il4_il10_adj_sex_age_L[6]), 
            " ", " ", bonpval(t3_ratio_il5_il10_adj_sex_age_L[6]), " ", " ", bonpval(t3_ratio_il13_il10_adj_sex_age_L[6]), " ", " ", bonpval(t3_ratio_il17_il10_adj_sex_age_L[6]), 
            " ", " ", bonpval(t3_ratio_il21_il10_adj_sex_age_L[6]), " ", " ", bonpval(t3_ratio_il2_il10_adj_sex_age_L[6]), " ", " ", bonpval(t3_ratio_gmc_il10_adj_sex_age_L[6]), 
            " ", " ", bonpval(t3_ratio_il12_il4_adj_sex_age_L[6]), " ", " ", bonpval(t3_ratio_ifn_il4_adj_sex_age_L[6]), " ", " ", bonpval(t3_ratio_il12_il5_adj_sex_age_L[6]), 
            " ", " ", bonpval(t3_ratio_ifn_il5_adj_sex_age_L[6]), " ", " ", bonpval(t3_ratio_il12_il13_adj_sex_age_L[6]), " ", " ", bonpval(t3_ratio_ifn_il13_adj_sex_age_L[6]), 
            " ", " ", bonpval(t3_ratio_il12_il17_adj_sex_age_L[6]), " ", " ", bonpval(t3_ratio_ifn_il17_adj_sex_age_L[6]), " ", " ", bonpval(t3_ratio_il12_il21_adj_sex_age_L[6]), 
            " ", " ", bonpval(t3_ratio_ifn_il21_adj_sex_age_L[6]), " ", " ", bonpval(t3_ratio_pro_il10_adj_sex_age_L[6]), " ", " ", bonpval(t3_ratio_th1_il10_adj_sex_age_L[6]), 
            " ", " ", bonpval(t3_ratio_th2_il10_adj_sex_age_L[6]), " ", " ", bonpval(t3_ratio_th17_il10_adj_sex_age_L[6]), " ", " ", bonpval(t3_ratio_th1_th2_adj_sex_age_L[6]), " ", " ", bonpval(t3_ratio_th1_th17_adj_sex_age_L[6])
)

t3_ratio_il1_il10_adj_sex_age_L <- round(t3_ratio_il1_il10_adj_sex_age_L, 2)
t3_ratio_il6_il10_adj_sex_age_L <- round(t3_ratio_il6_il10_adj_sex_age_L, 2)
t3_ratio_tnf_il10_adj_sex_age_L <- round(t3_ratio_tnf_il10_adj_sex_age_L, 2)
t3_ratio_il12_il10_adj_sex_age_L <- round(t3_ratio_il12_il10_adj_sex_age_L, 2)
t3_ratio_ifn_il10_adj_sex_age_L <- round(t3_ratio_ifn_il10_adj_sex_age_L, 2)
t3_ratio_il4_il10_adj_sex_age_L <- round(t3_ratio_il4_il10_adj_sex_age_L, 2)
t3_ratio_il5_il10_adj_sex_age_L <- round(t3_ratio_il5_il10_adj_sex_age_L, 2)
t3_ratio_il13_il10_adj_sex_age_L <- round(t3_ratio_il13_il10_adj_sex_age_L, 2)
t3_ratio_il17_il10_adj_sex_age_L <- round(t3_ratio_il17_il10_adj_sex_age_L, 2)
t3_ratio_il21_il10_adj_sex_age_L <- round(t3_ratio_il21_il10_adj_sex_age_L, 2)
t3_ratio_il2_il10_adj_sex_age_L <- round(t3_ratio_il2_il10_adj_sex_age_L, 2)
t3_ratio_gmc_il10_adj_sex_age_L <- round(t3_ratio_gmc_il10_adj_sex_age_L, 2)
t3_ratio_il12_il4_adj_sex_age_L <- round(t3_ratio_il12_il4_adj_sex_age_L, 2)
t3_ratio_ifn_il4_adj_sex_age_L <- round(t3_ratio_ifn_il4_adj_sex_age_L, 2)
t3_ratio_il12_il5_adj_sex_age_L <- round(t3_ratio_il12_il5_adj_sex_age_L, 2)
t3_ratio_ifn_il5_adj_sex_age_L <- round(t3_ratio_ifn_il5_adj_sex_age_L, 2)
t3_ratio_il12_il13_adj_sex_age_L <- round(t3_ratio_il12_il13_adj_sex_age_L, 2)
t3_ratio_ifn_il13_adj_sex_age_L <- round(t3_ratio_ifn_il13_adj_sex_age_L, 2)
t3_ratio_il12_il17_adj_sex_age_L <- round(t3_ratio_il12_il17_adj_sex_age_L, 2)
t3_ratio_ifn_il17_adj_sex_age_L <- round(t3_ratio_ifn_il17_adj_sex_age_L, 2)
t3_ratio_il12_il21_adj_sex_age_L <- round(t3_ratio_il12_il21_adj_sex_age_L, 2)
t3_ratio_ifn_il21_adj_sex_age_L <- round(t3_ratio_ifn_il21_adj_sex_age_L, 2)
t3_ratio_pro_il10_adj_sex_age_L <- round(t3_ratio_pro_il10_adj_sex_age_L, 2)
t3_ratio_th1_il10_adj_sex_age_L <- round(t3_ratio_th1_il10_adj_sex_age_L, 2)
t3_ratio_th2_il10_adj_sex_age_L <- round(t3_ratio_th2_il10_adj_sex_age_L, 2)
t3_ratio_th17_il10_adj_sex_age_L <- round(t3_ratio_th17_il10_adj_sex_age_L, 2)
t3_ratio_th1_th2_adj_sex_age_L <- round(t3_ratio_th1_th2_adj_sex_age_L, 2)
t3_ratio_th1_th17_adj_sex_age_L <- round(t3_ratio_th1_th17_adj_sex_age_L, 2)

asadjtbls5 <- c("95% CI", " ", " ", paste(t3_ratio_il1_il10_adj_sex_age_L[1], " (", t3_ratio_il1_il10_adj_sex_age_L[2], ", ", t3_ratio_il1_il10_adj_sex_age_L[3], ")", sep=""),
               " ", " ", paste(t3_ratio_il6_il10_adj_sex_age_L[1], " (", t3_ratio_il6_il10_adj_sex_age_L[2], ", ", t3_ratio_il6_il10_adj_sex_age_L[3], ")", sep=""),
               " ", " ", paste(t3_ratio_tnf_il10_adj_sex_age_L[1], " (", t3_ratio_tnf_il10_adj_sex_age_L[2], ", ", t3_ratio_tnf_il10_adj_sex_age_L[3], ")", sep=""),
               " ", " ", paste(t3_ratio_il12_il10_adj_sex_age_L[1], " (", t3_ratio_il12_il10_adj_sex_age_L[2], ", ", t3_ratio_il12_il10_adj_sex_age_L[3], ")", sep=""),
               " ", " ", paste(t3_ratio_ifn_il10_adj_sex_age_L[1], " (", t3_ratio_ifn_il10_adj_sex_age_L[2], ", ", t3_ratio_ifn_il10_adj_sex_age_L[3], ")", sep=""),
               " ", " ", paste(t3_ratio_il4_il10_adj_sex_age_L[1], " (", t3_ratio_il4_il10_adj_sex_age_L[2], ", ", t3_ratio_il4_il10_adj_sex_age_L[3], ")", sep=""),
               " ", " ", paste(t3_ratio_il5_il10_adj_sex_age_L[1], " (", t3_ratio_il5_il10_adj_sex_age_L[2], ", ", t3_ratio_il5_il10_adj_sex_age_L[3], ")", sep=""),
               " ", " ", paste(t3_ratio_il13_il10_adj_sex_age_L[1], " (", t3_ratio_il13_il10_adj_sex_age_L[2], ", ", t3_ratio_il13_il10_adj_sex_age_L[3], ")", sep=""),
               " ", " ", paste(t3_ratio_il17_il10_adj_sex_age_L[1], " (", t3_ratio_il17_il10_adj_sex_age_L[2], ", ", t3_ratio_il17_il10_adj_sex_age_L[3], ")", sep=""),
               " ", " ", paste(t3_ratio_il21_il10_adj_sex_age_L[1], " (", t3_ratio_il21_il10_adj_sex_age_L[2], ", ", t3_ratio_il21_il10_adj_sex_age_L[3], ")", sep=""),
               " ", " ", paste(t3_ratio_il2_il10_adj_sex_age_L[1], " (", t3_ratio_il2_il10_adj_sex_age_L[2], ", ", t3_ratio_il2_il10_adj_sex_age_L[3], ")", sep=""),
               " ", " ", paste(t3_ratio_gmc_il10_adj_sex_age_L[1], " (", t3_ratio_gmc_il10_adj_sex_age_L[2], ", ", t3_ratio_gmc_il10_adj_sex_age_L[3], ")", sep=""),
               " ", " ", paste(t3_ratio_il12_il4_adj_sex_age_L[1], " (", t3_ratio_il12_il4_adj_sex_age_L[2], ", ", t3_ratio_il12_il4_adj_sex_age_L[3], ")", sep=""),
               " ", " ", paste(t3_ratio_ifn_il4_adj_sex_age_L[1], " (", t3_ratio_ifn_il4_adj_sex_age_L[2], ", ", t3_ratio_ifn_il4_adj_sex_age_L[3], ")", sep=""),
               " ", " ", paste(t3_ratio_il12_il5_adj_sex_age_L[1], " (", t3_ratio_il12_il5_adj_sex_age_L[2], ", ", t3_ratio_il12_il5_adj_sex_age_L[3], ")", sep=""),
               " ", " ", paste(t3_ratio_ifn_il5_adj_sex_age_L[1], " (", t3_ratio_ifn_il5_adj_sex_age_L[2], ", ", t3_ratio_ifn_il5_adj_sex_age_L[3], ")", sep=""),
               " ", " ", paste(t3_ratio_il12_il13_adj_sex_age_L[1], " (", t3_ratio_il12_il13_adj_sex_age_L[2], ", ", t3_ratio_il12_il13_adj_sex_age_L[3], ")", sep=""),
               " ", " ", paste(t3_ratio_ifn_il13_adj_sex_age_L[1], " (", t3_ratio_ifn_il13_adj_sex_age_L[2], ", ", t3_ratio_ifn_il13_adj_sex_age_L[3], ")", sep=""),
               " ", " ", paste(t3_ratio_il12_il17_adj_sex_age_L[1], " (", t3_ratio_il12_il17_adj_sex_age_L[2], ", ", t3_ratio_il12_il17_adj_sex_age_L[3], ")", sep=""),
               " ", " ", paste(t3_ratio_ifn_il17_adj_sex_age_L[1], " (", t3_ratio_ifn_il17_adj_sex_age_L[2], ", ", t3_ratio_ifn_il17_adj_sex_age_L[3], ")", sep=""),
               " ", " ", paste(t3_ratio_il12_il21_adj_sex_age_L[1], " (", t3_ratio_il12_il21_adj_sex_age_L[2], ", ", t3_ratio_il12_il21_adj_sex_age_L[3], ")", sep=""),
               " ", " ", paste(t3_ratio_ifn_il21_adj_sex_age_L[1], " (", t3_ratio_ifn_il21_adj_sex_age_L[2], ", ", t3_ratio_ifn_il21_adj_sex_age_L[3], ")", sep=""),
               " ", " ", paste(t3_ratio_pro_il10_adj_sex_age_L[1], " (", t3_ratio_pro_il10_adj_sex_age_L[2], ", ", t3_ratio_pro_il10_adj_sex_age_L[3], ")", sep=""),
               " ", " ", paste(t3_ratio_th1_il10_adj_sex_age_L[1], " (", t3_ratio_th1_il10_adj_sex_age_L[2], ", ", t3_ratio_th1_il10_adj_sex_age_L[3], ")", sep=""),
               " ", " ", paste(t3_ratio_th2_il10_adj_sex_age_L[1], " (", t3_ratio_th2_il10_adj_sex_age_L[2], ", ", t3_ratio_th2_il10_adj_sex_age_L[3], ")", sep=""),
               " ", " ", paste(t3_ratio_th17_il10_adj_sex_age_L[1], " (", t3_ratio_th17_il10_adj_sex_age_L[2], ", ", t3_ratio_th17_il10_adj_sex_age_L[3], ")", sep=""),
               " ", " ", paste(t3_ratio_th1_th2_adj_sex_age_L[1], " (", t3_ratio_th1_th2_adj_sex_age_L[2], ", ", t3_ratio_th1_th2_adj_sex_age_L[3], ")", sep=""),
               " ", " ", paste(t3_ratio_th1_th17_adj_sex_age_L[1], " (", t3_ratio_th1_th17_adj_sex_age_L[2], ", ", t3_ratio_th1_th17_adj_sex_age_L[3], ")", sep=""))

ipcws5<-c("95% CI", " ", " ", maketblvalue(ratio_il1_il10_t3_adj_ipcw_L$`unlist(ratio_il1_il10_t3_adj_ipcw$estimates$ATE)`),
          " ", " ", maketblvalue(ratio_il6_il10_t3_adj_ipcw_L$`unlist(ratio_il6_il10_t3_adj_ipcw$estimates$ATE)`),
          " ", " ", maketblvalue(ratio_tnf_il10_t3_adj_ipcw_L$`unlist(ratio_tnf_il10_t3_adj_ipcw$estimates$ATE)`),
          " ", " ", maketblvalue(ratio_il12_il10_t3_adj_ipcw_L$`unlist(ratio_il12_il10_t3_adj_ipcw$estimates$ATE)`),
          " ", " ", maketblvalue(ratio_ifn_il10_t3_adj_ipcw_L$`unlist(ratio_ifn_il10_t3_adj_ipcw$estimates$ATE)`),
          " ", " ", maketblvalue(ratio_il4_il10_t3_adj_ipcw_L$`unlist(ratio_il4_il10_t3_adj_ipcw$estimates$ATE)`),
          " ", " ", maketblvalue(ratio_il5_il10_t3_adj_ipcw_L$`unlist(ratio_il5_il10_t3_adj_ipcw$estimates$ATE)`),
          " ", " ", maketblvalue(ratio_il13_il10_t3_adj_ipcw_L$`unlist(ratio_il13_il10_t3_adj_ipcw$estimates$ATE)`),
          " ", " ", maketblvalue(ratio_il17_il10_t3_adj_ipcw_L$`unlist(ratio_il17_il10_t3_adj_ipcw$estimates$ATE)`),
          " ", " ", maketblvalue(ratio_il21_il10_t3_adj_ipcw_L$`unlist(ratio_il21_il10_t3_adj_ipcw$estimates$ATE)`),
          " ", " ", maketblvalue(ratio_il2_il10_t3_adj_ipcw_L$`unlist(ratio_il2_il10_t3_adj_ipcw$estimates$ATE)`),
          " ", " ", maketblvalue(ratio_gmc_il10_t3_adj_ipcw_L$`unlist(ratio_gmc_il10_t3_adj_ipcw$estimates$ATE)`),
          " ", " ", maketblvalue(ratio_il12_il4_t3_adj_ipcw_L$`unlist(ratio_il12_il4_t3_adj_ipcw$estimates$ATE)`),
          " ", " ", maketblvalue(ratio_ifn_il4_t3_adj_ipcw_L$`unlist(ratio_ifn_il4_t3_adj_ipcw$estimates$ATE)`),
          " ", " ", maketblvalue(ratio_il12_il5_t3_adj_ipcw_L$`unlist(ratio_il12_il5_t3_adj_ipcw$estimates$ATE)`),
          " ", " ", maketblvalue(ratio_ifn_il5_t3_adj_ipcw_L$`unlist(ratio_ifn_il5_t3_adj_ipcw$estimates$ATE)`),
          " ", " ", maketblvalue(ratio_il12_il13_t3_adj_ipcw_L$`unlist(ratio_il12_il13_t3_adj_ipcw$estimates$ATE)`),
          " ", " ", maketblvalue(ratio_ifn_il13_t3_adj_ipcw_L$`unlist(ratio_ifn_il13_t3_adj_ipcw$estimates$ATE)`),
          " ", " ", maketblvalue(ratio_il12_il17_t3_adj_ipcw_L$`unlist(ratio_il12_il17_t3_adj_ipcw$estimates$ATE)`),
          " ", " ", maketblvalue(ratio_ifn_il17_t3_adj_ipcw_L$`unlist(ratio_ifn_il17_t3_adj_ipcw$estimates$ATE)`),
          " ", " ", maketblvalue(ratio_il12_il21_t3_adj_ipcw_L$`unlist(ratio_il12_il21_t3_adj_ipcw$estimates$ATE)`),
          " ", " ", maketblvalue(ratio_ifn_il21_t3_adj_ipcw_L$`unlist(ratio_ifn_il21_t3_adj_ipcw$estimates$ATE)`),
          " ", " ", maketblvalue(ratio_pro_il10_t3_adj_ipcw_L$`unlist(ratio_pro_il10_t3_adj_ipcw$estimates$ATE)`),
          " ", " ", maketblvalue(ratio_th1_il10_t3_adj_ipcw_L$`unlist(ratio_th1_il10_t3_adj_ipcw$estimates$ATE)`),
          " ", " ", maketblvalue(ratio_th2_il10_t3_adj_ipcw_L$`unlist(ratio_th2_il10_t3_adj_ipcw$estimates$ATE)`),
          " ", " ", maketblvalue(ratio_th17_il10_t3_adj_ipcw_L$`unlist(ratio_th17_il10_t3_adj_ipcw$estimates$ATE)`),
          " ", " ", maketblvalue(ratio_th1_th2_t3_adj_ipcw_L$`unlist(ratio_th1_th2_t3_adj_ipcw$estimates$ATE)`),
          " ", " ", maketblvalue(ratio_th1_th17_t3_adj_ipcw_L$`unlist(ratio_th1_th17_t3_adj_ipcw$estimates$ATE)`))

rawipcwpvals5<-c("P-value", " ", " ", rawpval(ratio_il1_il10_t3_adj_ipcw_L[5,1]),
          " ", " ", rawpval(ratio_il6_il10_t3_adj_ipcw_L[5,1]),
          " ", " ", rawpval(ratio_tnf_il10_t3_adj_ipcw_L[5,1]),
          " ", " ", rawpval(ratio_il12_il10_t3_adj_ipcw_L[5,1]),
          " ", " ", rawpval(ratio_ifn_il10_t3_adj_ipcw_L[5,1]),
          " ", " ", rawpval(ratio_il4_il10_t3_adj_ipcw_L[5,1]),
          " ", " ", rawpval(ratio_il5_il10_t3_adj_ipcw_L[5,1]),
          " ", " ", rawpval(ratio_il13_il10_t3_adj_ipcw_L[5,1]),
          " ", " ", rawpval(ratio_il17_il10_t3_adj_ipcw_L[5,1]),
          " ", " ", rawpval(ratio_il21_il10_t3_adj_ipcw_L[5,1]),
          " ", " ", rawpval(ratio_il2_il10_t3_adj_ipcw_L[5,1]),
          " ", " ", rawpval(ratio_gmc_il10_t3_adj_ipcw_L[5,1]),
          " ", " ", rawpval(ratio_il12_il4_t3_adj_ipcw_L[5,1]),
          " ", " ", rawpval(ratio_ifn_il4_t3_adj_ipcw_L[5,1]),
          " ", " ", rawpval(ratio_il12_il5_t3_adj_ipcw_L[5,1]),
          " ", " ", rawpval(ratio_ifn_il5_t3_adj_ipcw_L[5,1]),
          " ", " ", rawpval(ratio_il12_il13_t3_adj_ipcw_L[5,1]),
          " ", " ", rawpval(ratio_ifn_il13_t3_adj_ipcw_L[5,1]),
          " ", " ", rawpval(ratio_il12_il17_t3_adj_ipcw_L[5,1]),
          " ", " ", rawpval(ratio_ifn_il17_t3_adj_ipcw_L[5,1]),
          " ", " ", rawpval(ratio_il12_il21_t3_adj_ipcw_L[5,1]),
          " ", " ", rawpval(ratio_ifn_il21_t3_adj_ipcw_L[5,1]),
          " ", " ", rawpval(ratio_pro_il10_t3_adj_ipcw_L[5,1]),
          " ", " ", rawpval(ratio_th1_il10_t3_adj_ipcw_L[5,1]),
          " ", " ", rawpval(ratio_th2_il10_t3_adj_ipcw_L[5,1]),
          " ", " ", rawpval(ratio_th17_il10_t3_adj_ipcw_L[5,1]),
          " ", " ", rawpval(ratio_th1_th2_t3_adj_ipcw_L[5,1]),
          " ", " ", rawpval(ratio_th1_th17_t3_adj_ipcw_L[5,1]))

ipcwpvals5<-c("Bonferroni P-value", " ", " ", bonpval(ratio_il1_il10_t3_adj_ipcw_L[5,1]),
          " ", " ", bonpval(ratio_il6_il10_t3_adj_ipcw_L[5,1]),
          " ", " ", bonpval(ratio_tnf_il10_t3_adj_ipcw_L[5,1]),
          " ", " ", bonpval(ratio_il12_il10_t3_adj_ipcw_L[5,1]),
          " ", " ", bonpval(ratio_ifn_il10_t3_adj_ipcw_L[5,1]),
          " ", " ", bonpval(ratio_il4_il10_t3_adj_ipcw_L[5,1]),
          " ", " ", bonpval(ratio_il5_il10_t3_adj_ipcw_L[5,1]),
          " ", " ", bonpval(ratio_il13_il10_t3_adj_ipcw_L[5,1]),
          " ", " ", bonpval(ratio_il17_il10_t3_adj_ipcw_L[5,1]),
          " ", " ", bonpval(ratio_il21_il10_t3_adj_ipcw_L[5,1]),
          " ", " ", bonpval(ratio_il2_il10_t3_adj_ipcw_L[5,1]),
          " ", " ", bonpval(ratio_gmc_il10_t3_adj_ipcw_L[5,1]),
          " ", " ", bonpval(ratio_il12_il4_t3_adj_ipcw_L[5,1]),
          " ", " ", bonpval(ratio_ifn_il4_t3_adj_ipcw_L[5,1]),
          " ", " ", bonpval(ratio_il12_il5_t3_adj_ipcw_L[5,1]),
          " ", " ", bonpval(ratio_ifn_il5_t3_adj_ipcw_L[5,1]),
          " ", " ", bonpval(ratio_il12_il13_t3_adj_ipcw_L[5,1]),
          " ", " ", bonpval(ratio_ifn_il13_t3_adj_ipcw_L[5,1]),
          " ", " ", bonpval(ratio_il12_il17_t3_adj_ipcw_L[5,1]),
          " ", " ", bonpval(ratio_ifn_il17_t3_adj_ipcw_L[5,1]),
          " ", " ", bonpval(ratio_il12_il21_t3_adj_ipcw_L[5,1]),
          " ", " ", bonpval(ratio_ifn_il21_t3_adj_ipcw_L[5,1]),
          " ", " ", bonpval(ratio_pro_il10_t3_adj_ipcw_L[5,1]),
          " ", " ", bonpval(ratio_th1_il10_t3_adj_ipcw_L[5,1]),
          " ", " ", bonpval(ratio_th2_il10_t3_adj_ipcw_L[5,1]),
          " ", " ", bonpval(ratio_th17_il10_t3_adj_ipcw_L[5,1]),
          " ", " ", bonpval(ratio_th1_th2_t3_adj_ipcw_L[5,1]),
          " ", " ", bonpval(ratio_th1_th17_t3_adj_ipcw_L[5,1]))


tbls5<-data.table(
  " " = tbl3[[1]],
  " " = tbl3[[2]],
  " " = tbl3[[3]],
  " " = tbl3[[4]],
  " " = tbl3[[5]],
  " " = tbl3[[6]],
  "Unadjusted difference: Intervention vs. Control" = tbl3[[7]],
  " " = rawpval3,
  " " = as_bonf_col(tbl3[[8]]),
  "Age- and sex- adjusted difference: Intervention vs. Control" = asadjtbls5,
  " " = rawpvals5,
  " " = pvals5,
  "Fully adjusted difference: Intervention vs. Control" = tbl3[[9]],
  " " = rawadjpval3,
  " " = as_bonf_col(tbl3[[10]]),
  "IPCW adjusted difference: Intervention vs. Control" = ipcws5,
  " " = rawipcwpvals5,
  " " = ipcwpvals5)

write.csv(tbls5, file=here('tables/supplementary/immune_supptable5.csv'))
print(xtable(tbls5), type="html", file=here("tables/supplementary/immune_supptable5.html"))



#### TABLE S6 ####

outcomes6 <- c("Outcome, Arm", paste("Ln \u0394IL-1", "\u03b2", " (pg/ml)", sep=""), "Control", "Nutrition + WSH", 
                 "Ln \u0394IL-6 (pg/ml)", "Control", "Nutrition + WSH", 
                 paste("Ln \u0394TNF-", "\u03b1", " (pg/ml)", sep=""), "Control", "Nutrition + WSH",
                 "Ln \u0394CRP (mg/L)", "Control", "Nutrition + WSH", 
                 "Ln \u0394IL-12 (pg/ml)", "Control", "Nutrition + WSH", 
                 paste("Ln \u0394IFN-", "\u03b3", " (pg/ml)", sep=""), "Control", "Nutrition + WSH", 
                 "Ln \u0394IL-4 (pg/ml)", "Control", "Nutrition + WSH", 
                 "Ln \u0394IL-5 (pg/ml)", "Control", "Nutrition + WSH", 
                 "Ln \u0394IL-13 (pg/ml)", "Control", "Nutrition + WSH", 
                 "Ln \u0394IL-17A (pg/ml)", "Control", "Nutrition + WSH", 
                 "Ln \u0394IL-21 (pg/ml)", "Control", "Nutrition + WSH", 
                 "Ln \u0394IL-10 (pg/ml)", "Control", "Nutrition + WSH", 
                 "Ln \u0394IL-2 (pg/ml)", "Control", "Nutrition + WSH", 
                 "Ln \u0394GM-CSF (pg/ml)", "Control", "Nutrition + WSH", 
                 "Ln \u0394AGP (g/L)", "Control", "Nutrition + WSH", 
                 paste("Ln \u0394IGF-1 (", "\u03bc", "g/L)", sep=""), "Control", "Nutrition + WSH")

Ns6 <- c("N", " ", as.character(d23_ln_il1_N_tr$d23_ln_il1_N_tr[1]), as.character(d23_ln_il1_N_tr$d23_ln_il1_N_tr[2]), 
           " ", as.character(d23_ln_il6_N_tr$d23_ln_il6_N_tr[1]), as.character(d23_ln_il6_N_tr$d23_ln_il6_N_tr[2]),
         " ", as.character(d23_ln_tnf_N_tr$d23_ln_tnf_N_tr[1]), as.character(d23_ln_tnf_N_tr$d23_ln_tnf_N_tr[2]),
         " ", as.character(d23_ln_crp_N_tr$d23_ln_crp_N_tr[1]), as.character(d23_ln_crp_N_tr$d23_ln_crp_N_tr[2]),
         " ", as.character(d23_ln_il12_N_tr$d23_ln_il12_N_tr[1]), as.character(d23_ln_il12_N_tr$d23_ln_il12_N_tr[2]),
           " ", as.character(d23_ln_ifn_N_tr$d23_ln_ifn_N_tr[1]), as.character(d23_ln_ifn_N_tr$d23_ln_ifn_N_tr[2]),
           " ", as.character(d23_ln_il4_N_tr$d23_ln_il4_N_tr[1]), as.character(d23_ln_il4_N_tr$d23_ln_il4_N_tr[2]),
           " ", as.character(d23_ln_il5_N_tr$d23_ln_il5_N_tr[1]), as.character(d23_ln_il5_N_tr$d23_ln_il5_N_tr[2]),
           " ", as.character(d23_ln_il13_N_tr$d23_ln_il13_N_tr[1]), as.character(d23_ln_il13_N_tr$d23_ln_il13_N_tr[2]),
           " ", as.character(d23_ln_il17_N_tr$d23_ln_il17_N_tr[1]), as.character(d23_ln_il17_N_tr$d23_ln_il17_N_tr[2]),
           " ", as.character(d23_ln_il21_N_tr$d23_ln_il21_N_tr[1]), as.character(d23_ln_il21_N_tr$d23_ln_il21_N_tr[2]),
           " ", as.character(d23_ln_il10_N_tr$d23_ln_il10_N_tr[1]), as.character(d23_ln_il10_N_tr$d23_ln_il10_N_tr[2]),
           " ", as.character(d23_ln_il2_N_tr$d23_ln_il2_N_tr[1]), as.character(d23_ln_il2_N_tr$d23_ln_il2_N_tr[2]),
         " ", as.character(d23_ln_gmc_N_tr$d23_ln_gmc_N_tr[1]), as.character(d23_ln_gmc_N_tr$d23_ln_gmc_N_tr[2]),
         " ", as.character(d23_ln_agp_N_tr$d23_ln_agp_N_tr[1]), as.character(d23_ln_agp_N_tr$d23_ln_agp_N_tr[2]),
         " ", as.character(d23_ln_igf_N_tr$d23_ln_igf_N_tr[1]), as.character(d23_ln_igf_N_tr$d23_ln_igf_N_tr[2]))

absmeans6 <- c("Absolute Mean", " ", as.character(round(abs_d23_il1_N_tr$mean[1], 2)), as.character(round(abs_d23_il1_N_tr$mean[2], 2)), 
                 " ", as.character(round(abs_d23_il6_N_tr$mean[1], 2)), as.character(round(abs_d23_il6_N_tr$mean[2], 2)),
               " ", as.character(round(abs_d23_tnf_N_tr$mean[1], 2)), as.character(round(abs_d23_tnf_N_tr$mean[2], 2)),
               " ", as.character(round(abs_d23_crp_N_tr$mean[1], 2)), as.character(round(abs_d23_crp_N_tr$mean[2], 2)),
               " ", as.character(round(abs_d23_il12_N_tr$mean[1], 2)), as.character(round(abs_d23_il12_N_tr$mean[2], 2)),
                 " ", as.character(round(abs_d23_ifn_N_tr$mean[1], 2)), as.character(round(abs_d23_ifn_N_tr$mean[2], 2)),
                 " ", as.character(round(abs_d23_il4_N_tr$mean[1], 2)), as.character(round(abs_d23_il4_N_tr$mean[2], 2)),
                 " ", as.character(round(abs_d23_il5_N_tr$mean[1], 2)), as.character(round(abs_d23_il5_N_tr$mean[2], 2)),
                 " ", as.character(round(abs_d23_il13_N_tr$mean[1], 2)), as.character(round(abs_d23_il13_N_tr$mean[2], 2)),
                 " ", as.character(round(abs_d23_il17_N_tr$mean[1], 2)), as.character(round(abs_d23_il17_N_tr$mean[2], 2)),
                 " ", as.character(round(abs_d23_il21_N_tr$mean[1], 2)), as.character(round(abs_d23_il21_N_tr$mean[2], 2)),
                 " ", as.character(round(abs_d23_il10_N_tr$mean[1], 2)), as.character(round(abs_d23_il10_N_tr$mean[2], 2)),
                 " ", as.character(round(abs_d23_il2_N_tr$mean[1], 2)), as.character(round(abs_d23_il2_N_tr$mean[2], 2)),
               " ", as.character(round(abs_d23_gmc_N_tr$mean[1], 2)), as.character(round(abs_d23_gmc_N_tr$mean[2], 2)),
               " ", as.character(round(abs_d23_agp_N_tr$mean[1], 2)), as.character(round(abs_d23_agp_N_tr$mean[2], 2)),
               " ", as.character(round(abs_d23_igf_N_tr$mean[1], 2)), as.character(round(abs_d23_igf_N_tr$mean[2], 2)))

abssds6 <- c("Absolute SD", " ", as.character(round(abs_d23_il1_N_tr$sd[1], 2)), as.character(round(abs_d23_il1_N_tr$sd[2], 2)), 
               " ", as.character(round(abs_d23_il6_N_tr$sd[1], 2)), as.character(round(abs_d23_il6_N_tr$sd[2], 2)),
             " ", as.character(round(abs_d23_tnf_N_tr$sd[1], 2)), as.character(round(abs_d23_tnf_N_tr$sd[2], 2)),
             " ", as.character(round(abs_d23_crp_N_tr$sd[1], 2)), as.character(round(abs_d23_crp_N_tr$sd[2], 2)),
             " ", as.character(round(abs_d23_il12_N_tr$sd[1], 2)), as.character(round(abs_d23_il12_N_tr$sd[2], 2)),
               " ", as.character(round(abs_d23_ifn_N_tr$sd[1], 2)), as.character(round(abs_d23_ifn_N_tr$sd[2], 2)),
               " ", as.character(round(abs_d23_il4_N_tr$sd[1], 2)), as.character(round(abs_d23_il4_N_tr$sd[2], 2)),
               " ", as.character(round(abs_d23_il5_N_tr$sd[1], 2)), as.character(round(abs_d23_il5_N_tr$sd[2], 2)),
               " ", as.character(round(abs_d23_il13_N_tr$sd[1], 2)), as.character(round(abs_d23_il13_N_tr$sd[2], 2)),
               " ", as.character(round(abs_d23_il17_N_tr$sd[1], 2)), as.character(round(abs_d23_il17_N_tr$sd[2], 2)),
               " ", as.character(round(abs_d23_il21_N_tr$sd[1], 2)), as.character(round(abs_d23_il21_N_tr$sd[2], 2)),
               " ", as.character(round(abs_d23_il10_N_tr$sd[1], 2)), as.character(round(abs_d23_il10_N_tr$sd[2], 2)),
               " ", as.character(round(abs_d23_il2_N_tr$sd[1], 2)), as.character(round(abs_d23_il2_N_tr$sd[2], 2)),
             " ", as.character(round(abs_d23_gmc_N_tr$sd[1], 2)), as.character(round(abs_d23_gmc_N_tr$sd[2], 2)),
             " ", as.character(round(abs_d23_agp_N_tr$sd[1], 2)), as.character(round(abs_d23_agp_N_tr$sd[2], 2)),
             " ", as.character(round(abs_d23_igf_N_tr$sd[1], 2)), as.character(round(abs_d23_igf_N_tr$sd[2], 2)))

means6 <- c("Mean", " ", as.character(round(d23_ln_il1_N_tr$mean[1], 2)), as.character(round(d23_ln_il1_N_tr$mean[2], 2)), 
              " ", as.character(round(d23_ln_il6_N_tr$mean[1], 2)), as.character(round(d23_ln_il6_N_tr$mean[2], 2)),
            " ", as.character(round(d23_ln_tnf_N_tr$mean[1], 2)), as.character(round(d23_ln_tnf_N_tr$mean[2], 2)),
            " ", as.character(round(d23_ln_crp_N_tr$mean[1], 2)), as.character(round(d23_ln_crp_N_tr$mean[2], 2)),
            " ", as.character(round(d23_ln_il12_N_tr$mean[1], 2)), as.character(round(d23_ln_il12_N_tr$mean[2], 2)),
              " ", as.character(round(d23_ln_ifn_N_tr$mean[1], 2)), as.character(round(d23_ln_ifn_N_tr$mean[2], 2)),
              " ", as.character(round(d23_ln_il4_N_tr$mean[1], 2)), as.character(round(d23_ln_il4_N_tr$mean[2], 2)),
              " ", as.character(round(d23_ln_il5_N_tr$mean[1], 2)), as.character(round(d23_ln_il5_N_tr$mean[2], 2)),
              " ", as.character(round(d23_ln_il13_N_tr$mean[1], 2)), as.character(round(d23_ln_il13_N_tr$mean[2], 2)),
              " ", as.character(round(d23_ln_il17_N_tr$mean[1], 2)), as.character(round(d23_ln_il17_N_tr$mean[2], 2)),
              " ", as.character(round(d23_ln_il21_N_tr$mean[1], 2)), as.character(round(d23_ln_il21_N_tr$mean[2], 2)),
              " ", as.character(round(d23_ln_il10_N_tr$mean[1], 2)), as.character(round(d23_ln_il10_N_tr$mean[2], 2)),
              " ", as.character(round(d23_ln_il2_N_tr$mean[1], 2)), as.character(round(d23_ln_il2_N_tr$mean[2], 2)),
            " ", as.character(round(d23_ln_gmc_N_tr$mean[1], 2)), as.character(round(d23_ln_gmc_N_tr$mean[2], 2)),
            " ", as.character(round(d23_ln_agp_N_tr$mean[1], 2)), as.character(round(d23_ln_agp_N_tr$mean[2], 2)),
            " ", as.character(round(d23_ln_igf_N_tr$mean[1], 2)), as.character(round(d23_ln_igf_N_tr$mean[2], 2)))

sds6 <- c("SD", " ", as.character(round(d23_ln_il1_N_tr$sd[1], 2)), as.character(round(d23_ln_il1_N_tr$sd[2], 2)), 
            " ", as.character(round(d23_ln_il6_N_tr$sd[1], 2)), as.character(round(d23_ln_il6_N_tr$sd[2], 2)),
          " ", as.character(round(d23_ln_tnf_N_tr$sd[1], 2)), as.character(round(d23_ln_tnf_N_tr$sd[2], 2)),
          " ", as.character(round(d23_ln_crp_N_tr$sd[1], 2)), as.character(round(d23_ln_crp_N_tr$sd[2], 2)),
          " ", as.character(round(d23_ln_il12_N_tr$sd[1], 2)), as.character(round(d23_ln_il12_N_tr$sd[2], 2)),
            " ", as.character(round(d23_ln_ifn_N_tr$sd[1], 2)), as.character(round(d23_ln_ifn_N_tr$sd[2], 2)),
            " ", as.character(round(d23_ln_il4_N_tr$sd[1], 2)), as.character(round(d23_ln_il4_N_tr$sd[2], 2)),
            " ", as.character(round(d23_ln_il5_N_tr$sd[1], 2)), as.character(round(d23_ln_il5_N_tr$sd[2], 2)),
            " ", as.character(round(d23_ln_il13_N_tr$sd[1], 2)), as.character(round(d23_ln_il13_N_tr$sd[2], 2)),
            " ", as.character(round(d23_ln_il17_N_tr$sd[1], 2)), as.character(round(d23_ln_il17_N_tr$sd[2], 2)),
            " ", as.character(round(d23_ln_il21_N_tr$sd[1], 2)), as.character(round(d23_ln_il21_N_tr$sd[2], 2)),
            " ", as.character(round(d23_ln_il10_N_tr$sd[1], 2)), as.character(round(d23_ln_il10_N_tr$sd[2], 2)),
            " ", as.character(round(d23_ln_il2_N_tr$sd[1], 2)), as.character(round(d23_ln_il2_N_tr$sd[2], 2)),
          " ", as.character(round(d23_ln_gmc_N_tr$sd[1], 2)), as.character(round(d23_ln_gmc_N_tr$sd[2], 2)),
          " ", as.character(round(d23_ln_agp_N_tr$sd[1], 2)), as.character(round(d23_ln_agp_N_tr$sd[2], 2)),
          " ", as.character(round(d23_ln_igf_N_tr$sd[1], 2)), as.character(round(d23_ln_igf_N_tr$sd[2], 2)))

rawpvals6 <- c("P-value", " ", " ", rawpval(d23_ln_il1_unadj_L[6]), " ", " ", rawpval(d23_ln_il6_unadj_L[6]), " ", " ", rawpval(d23_ln_tnf_unadj_L[6]),
            " ", " ", rawpval(d23_ln_crp_unadj_L[6]), " ", " ", rawpval(d23_ln_il12_unadj_L[6]), " ", " ", rawpval(d23_ln_ifn_unadj_L[6]),
            " ", " ", rawpval(d23_ln_il4_unadj_L[6]), " ", " ", rawpval(d23_ln_il5_unadj_L[6]), " ", " ", rawpval(d23_ln_il13_unadj_L[6]),
            " ", " ", rawpval(d23_ln_il17_unadj_L[6]), " ", " ", rawpval(d23_ln_il21_unadj_L[6]), " ", " ", rawpval(d23_ln_il10_unadj_L[6]),
            " ", " ", rawpval(d23_ln_il2_unadj_L[6]), " ", " ", rawpval(d23_ln_gmc_unadj_L[6]),
            " ", " ", rawpval(d23_ln_agp_unadj_L[6]), " ", " ", rawpval(d23_ln_igf_unadj_L[6]))

pvals6 <- c("Bonferroni P-value", " ", " ", bonpval(d23_ln_il1_unadj_L[6]), " ", " ", bonpval(d23_ln_il6_unadj_L[6]), " ", " ", bonpval(d23_ln_tnf_unadj_L[6]),
            " ", " ", bonpval(d23_ln_crp_unadj_L[6]), " ", " ", bonpval(d23_ln_il12_unadj_L[6]), " ", " ", bonpval(d23_ln_ifn_unadj_L[6]),
            " ", " ", bonpval(d23_ln_il4_unadj_L[6]), " ", " ", bonpval(d23_ln_il5_unadj_L[6]), " ", " ", bonpval(d23_ln_il13_unadj_L[6]),
            " ", " ", bonpval(d23_ln_il17_unadj_L[6]), " ", " ", bonpval(d23_ln_il21_unadj_L[6]), " ", " ", bonpval(d23_ln_il10_unadj_L[6]),
            " ", " ", bonpval(d23_ln_il2_unadj_L[6]), " ", " ", bonpval(d23_ln_gmc_unadj_L[6]),
            " ", " ", bonpval(d23_ln_agp_unadj_L[6]), " ", " ", bonpval(d23_ln_igf_unadj_L[6]))

d23_ln_crp_unadj_L <- round(d23_ln_crp_unadj_L, 2)
d23_ln_agp_unadj_L <- round(d23_ln_agp_unadj_L, 2)
d23_ln_il1_unadj_L <- round(d23_ln_il1_unadj_L, 2)
d23_ln_il6_unadj_L <- round(d23_ln_il6_unadj_L, 2)
d23_ln_tnf_unadj_L <- round(d23_ln_tnf_unadj_L, 2)
d23_ln_il12_unadj_L <- round(d23_ln_il12_unadj_L, 2)
d23_ln_ifn_unadj_L <- round(d23_ln_ifn_unadj_L, 2)
d23_ln_il4_unadj_L <- round(d23_ln_il4_unadj_L, 2)
d23_ln_il5_unadj_L <- round(d23_ln_il5_unadj_L, 2)
d23_ln_il13_unadj_L <- round(d23_ln_il13_unadj_L, 2)
d23_ln_il17_unadj_L <- round(d23_ln_il17_unadj_L, 2)
d23_ln_il21_unadj_L <- round(d23_ln_il21_unadj_L, 2)
d23_ln_il10_unadj_L <- round(d23_ln_il10_unadj_L, 2)
d23_ln_il2_unadj_L <- round(d23_ln_il2_unadj_L, 2)
d23_ln_gmc_unadj_L <- round(d23_ln_gmc_unadj_L, 2)
d23_ln_igf_unadj_L <- round(d23_ln_igf_unadj_L, 2)

unadjs6 <- c("95% CI", " ", " ", paste(d23_ln_il1_unadj_L[1], " (", d23_ln_il1_unadj_L[2], ", ", d23_ln_il1_unadj_L[3], ")", sep=""),
               " ", " ", paste(d23_ln_il6_unadj_L[1], " (", d23_ln_il6_unadj_L[2], ", ", d23_ln_il6_unadj_L[3], ")", sep=""),
             " ", " ", paste(d23_ln_tnf_unadj_L[1], " (", d23_ln_tnf_unadj_L[2], ", ", d23_ln_tnf_unadj_L[3], ")", sep=""),
             " ", " ", paste(d23_ln_crp_unadj_L[1], " (", d23_ln_crp_unadj_L[2], ", ", d23_ln_crp_unadj_L[3], ")", sep=""),
             " ", " ", paste(d23_ln_il12_unadj_L[1], " (", d23_ln_il12_unadj_L[2], ", ", d23_ln_il12_unadj_L[3], ")", sep=""),
               " ", " ", paste(d23_ln_ifn_unadj_L[1], " (", d23_ln_ifn_unadj_L[2], ", ", d23_ln_ifn_unadj_L[3], ")", sep=""),
               " ", " ", paste(d23_ln_il4_unadj_L[1], " (", d23_ln_il4_unadj_L[2], ", ", d23_ln_il4_unadj_L[3], ")", sep=""),
               " ", " ", paste(d23_ln_il5_unadj_L[1], " (", d23_ln_il5_unadj_L[2], ", ", d23_ln_il5_unadj_L[3], ")", sep=""),
               " ", " ", paste(d23_ln_il13_unadj_L[1], " (", d23_ln_il13_unadj_L[2], ", ", d23_ln_il13_unadj_L[3], ")", sep=""),
               " ", " ", paste(d23_ln_il17_unadj_L[1], " (", d23_ln_il17_unadj_L[2], ", ", d23_ln_il17_unadj_L[3], ")", sep=""),
               " ", " ", paste(d23_ln_il21_unadj_L[1], " (", d23_ln_il21_unadj_L[2], ", ", d23_ln_il21_unadj_L[3], ")", sep=""),
               " ", " ", paste(d23_ln_il10_unadj_L[1], " (", d23_ln_il10_unadj_L[2], ", ", d23_ln_il10_unadj_L[3], ")", sep=""),
               " ", " ", paste(d23_ln_il2_unadj_L[1], " (", d23_ln_il2_unadj_L[2], ", ", d23_ln_il2_unadj_L[3], ")", sep=""),
             " ", " ", paste(d23_ln_gmc_unadj_L[1], " (", d23_ln_gmc_unadj_L[2], ", ", d23_ln_gmc_unadj_L[3], ")", sep=""),
             " ", " ", paste(d23_ln_agp_unadj_L[1], " (", d23_ln_agp_unadj_L[2], ", ", d23_ln_agp_unadj_L[3], ")", sep=""),
             " ", " ", paste(d23_ln_igf_unadj_L[1], " (", d23_ln_igf_unadj_L[2], ", ", d23_ln_igf_unadj_L[3], ")", sep=""))

rawasapvals6 <- c("P-value", " ", " ", rawpval(d23_ln_il1_adj_sex_age_L[6]), " ", " ", rawpval(d23_ln_il6_adj_sex_age_L[6]), " ", " ", rawpval(d23_ln_tnf_adj_sex_age_L[6]),
               " ", " ", rawpval(d23_ln_crp_adj_sex_age_L[6])," ", " ", rawpval(d23_ln_il12_adj_sex_age_L[6]),  " ", " ", rawpval(d23_ln_ifn_adj_sex_age_L[6]),
               " ", " ", rawpval(d23_ln_il4_adj_sex_age_L[6]), " ", " ", rawpval(d23_ln_il5_adj_sex_age_L[6]), " ", " ", rawpval(d23_ln_il13_adj_sex_age_L[6]),
               " ", " ", rawpval(d23_ln_il17_adj_sex_age_L[6]), " ", " ", rawpval(d23_ln_il21_adj_sex_age_L[6]), " ", " ", rawpval(d23_ln_il10_adj_sex_age_L[6]),
               " ", " ", rawpval(d23_ln_il2_adj_sex_age_L[6]), " ", " ", rawpval(d23_ln_gmc_adj_sex_age_L[6]),
               " ", " ", rawpval(d23_ln_agp_adj_sex_age_L[6]), " ", " ", rawpval(d23_ln_igf_adj_sex_age_L[6]))

asapvals6 <- c("Bonferroni P-value", " ", " ", bonpval(d23_ln_il1_adj_sex_age_L[6]), " ", " ", bonpval(d23_ln_il6_adj_sex_age_L[6]), " ", " ", bonpval(d23_ln_tnf_adj_sex_age_L[6]),
               " ", " ", bonpval(d23_ln_crp_adj_sex_age_L[6])," ", " ", bonpval(d23_ln_il12_adj_sex_age_L[6]),  " ", " ", bonpval(d23_ln_ifn_adj_sex_age_L[6]),
               " ", " ", bonpval(d23_ln_il4_adj_sex_age_L[6]), " ", " ", bonpval(d23_ln_il5_adj_sex_age_L[6]), " ", " ", bonpval(d23_ln_il13_adj_sex_age_L[6]),
               " ", " ", bonpval(d23_ln_il17_adj_sex_age_L[6]), " ", " ", bonpval(d23_ln_il21_adj_sex_age_L[6]), " ", " ", bonpval(d23_ln_il10_adj_sex_age_L[6]),
               " ", " ", bonpval(d23_ln_il2_adj_sex_age_L[6]), " ", " ", bonpval(d23_ln_gmc_adj_sex_age_L[6]),
               " ", " ", bonpval(d23_ln_agp_adj_sex_age_L[6]), " ", " ", bonpval(d23_ln_igf_adj_sex_age_L[6]))

d23_ln_crp_adj_sex_age_L <- round(d23_ln_crp_adj_sex_age_L, 2)
d23_ln_agp_adj_sex_age_L <- round(d23_ln_agp_adj_sex_age_L, 2)
d23_ln_il1_adj_sex_age_L <- round(d23_ln_il1_adj_sex_age_L, 2)
d23_ln_il6_adj_sex_age_L <- round(d23_ln_il6_adj_sex_age_L, 2)
d23_ln_tnf_adj_sex_age_L <- round(d23_ln_tnf_adj_sex_age_L, 2)
d23_ln_il12_adj_sex_age_L <- round(d23_ln_il12_adj_sex_age_L, 2)
d23_ln_ifn_adj_sex_age_L <- round(d23_ln_ifn_adj_sex_age_L, 2)
d23_ln_il4_adj_sex_age_L <- round(d23_ln_il4_adj_sex_age_L, 2)
d23_ln_il5_adj_sex_age_L <- round(d23_ln_il5_adj_sex_age_L, 2)
d23_ln_il13_adj_sex_age_L <- round(d23_ln_il13_adj_sex_age_L, 2)
d23_ln_il17_adj_sex_age_L <- round(d23_ln_il17_adj_sex_age_L, 2)
d23_ln_il21_adj_sex_age_L <- round(d23_ln_il21_adj_sex_age_L, 2)
d23_ln_il10_adj_sex_age_L <- round(d23_ln_il10_adj_sex_age_L, 2)
d23_ln_il2_adj_sex_age_L <- round(d23_ln_il2_adj_sex_age_L, 2)
d23_ln_gmc_adj_sex_age_L <- round(d23_ln_gmc_adj_sex_age_L, 2)
d23_ln_igf_adj_sex_age_L <- round(d23_ln_igf_adj_sex_age_L, 2)

asadjs6 <- c("95% CI", " ", " ", paste(d23_ln_il1_adj_sex_age_L[1], " (", d23_ln_il1_adj_sex_age_L[2], ", ", d23_ln_il1_adj_sex_age_L[3], ")", sep=""),
               " ", " ", paste(d23_ln_il6_adj_sex_age_L[1], " (", d23_ln_il6_adj_sex_age_L[2], ", ", d23_ln_il6_adj_sex_age_L[3], ")", sep=""),
             " ", " ", paste(d23_ln_tnf_adj_sex_age_L[1], " (", d23_ln_tnf_adj_sex_age_L[2], ", ", d23_ln_tnf_adj_sex_age_L[3], ")", sep=""),
             " ", " ", paste(d23_ln_crp_adj_sex_age_L[1], " (", d23_ln_crp_adj_sex_age_L[2], ", ", d23_ln_crp_adj_sex_age_L[3], ")", sep=""),
             " ", " ", paste(d23_ln_il12_adj_sex_age_L[1], " (", d23_ln_il12_adj_sex_age_L[2], ", ", d23_ln_il12_adj_sex_age_L[3], ")", sep=""),
               " ", " ", paste(d23_ln_ifn_adj_sex_age_L[1], " (", d23_ln_ifn_adj_sex_age_L[2], ", ", d23_ln_ifn_adj_sex_age_L[3], ")", sep=""),
               " ", " ", paste(d23_ln_il4_adj_sex_age_L[1], " (", d23_ln_il4_adj_sex_age_L[2], ", ", d23_ln_il4_adj_sex_age_L[3], ")", sep=""),
               " ", " ", paste(d23_ln_il5_adj_sex_age_L[1], " (", d23_ln_il5_adj_sex_age_L[2], ", ", d23_ln_il5_adj_sex_age_L[3], ")", sep=""),
               " ", " ", paste(d23_ln_il13_adj_sex_age_L[1], " (", d23_ln_il13_adj_sex_age_L[2], ", ", d23_ln_il13_adj_sex_age_L[3], ")", sep=""),
               " ", " ", paste(d23_ln_il17_adj_sex_age_L[1], " (", d23_ln_il17_adj_sex_age_L[2], ", ", d23_ln_il17_adj_sex_age_L[3], ")", sep=""),
               " ", " ", paste(d23_ln_il21_adj_sex_age_L[1], " (", d23_ln_il21_adj_sex_age_L[2], ", ", d23_ln_il21_adj_sex_age_L[3], ")", sep=""),
               " ", " ", paste(d23_ln_il10_adj_sex_age_L[1], " (", d23_ln_il10_adj_sex_age_L[2], ", ", d23_ln_il10_adj_sex_age_L[3], ")", sep=""),
               " ", " ", paste(d23_ln_il2_adj_sex_age_L[1], " (", d23_ln_il2_adj_sex_age_L[2], ", ", d23_ln_il2_adj_sex_age_L[3], ")", sep=""),
             " ", " ", paste(d23_ln_gmc_adj_sex_age_L[1], " (", d23_ln_gmc_adj_sex_age_L[2], ", ", d23_ln_gmc_adj_sex_age_L[3], ")", sep=""),
             " ", " ", paste(d23_ln_agp_adj_sex_age_L[1], " (", d23_ln_agp_adj_sex_age_L[2], ", ", d23_ln_agp_adj_sex_age_L[3], ")", sep=""),
             " ", " ", paste(d23_ln_igf_adj_sex_age_L[1], " (", d23_ln_igf_adj_sex_age_L[2], ", ", d23_ln_igf_adj_sex_age_L[3], ")", sep=""))

rawadjpvals6 <- c("P-value", " ", " ", rawpval(d23_ln_il1_adj_L[6]), " ", " ", rawpval(d23_ln_il6_adj_L[6]), " ", " ", rawpval(d23_ln_tnf_adj_L[6]),
               " ", " ", rawpval(d23_ln_crp_adj_L[6]), " ", " ", rawpval(d23_ln_il12_adj_L[6]), " ", " ", rawpval(d23_ln_ifn_adj_L[6]),
               " ", " ", rawpval(d23_ln_il4_adj_L[6]), " ", " ", rawpval(d23_ln_il5_adj_L[6]), " ", " ", rawpval(d23_ln_il13_adj_L[6]),
               " ", " ", rawpval(d23_ln_il17_adj_L[6]), " ", " ", rawpval(d23_ln_il21_adj_L[6]), " ", " ", rawpval(d23_ln_il10_adj_L[6]),
               " ", " ", rawpval(d23_ln_il2_adj_L[6]), " ", " ", rawpval(d23_ln_gmc_adj_L[6]),
               " ", " ", rawpval(d23_ln_agp_adj_L[6]), " ", " ", rawpval(d23_ln_igf_adj_L[6]))

adjpvals6 <- c("Bonferroni P-value", " ", " ", bonpval(d23_ln_il1_adj_L[6]), " ", " ", bonpval(d23_ln_il6_adj_L[6]), " ", " ", bonpval(d23_ln_tnf_adj_L[6]),
               " ", " ", bonpval(d23_ln_crp_adj_L[6]), " ", " ", bonpval(d23_ln_il12_adj_L[6]), " ", " ", bonpval(d23_ln_ifn_adj_L[6]),
               " ", " ", bonpval(d23_ln_il4_adj_L[6]), " ", " ", bonpval(d23_ln_il5_adj_L[6]), " ", " ", bonpval(d23_ln_il13_adj_L[6]),
               " ", " ", bonpval(d23_ln_il17_adj_L[6]), " ", " ", bonpval(d23_ln_il21_adj_L[6]), " ", " ", bonpval(d23_ln_il10_adj_L[6]),
               " ", " ", bonpval(d23_ln_il2_adj_L[6]), " ", " ", bonpval(d23_ln_gmc_adj_L[6]),
               " ", " ", bonpval(d23_ln_agp_adj_L[6]), " ", " ", bonpval(d23_ln_igf_adj_L[6]))

d23_ln_crp_adj_L <- round(d23_ln_crp_adj_L, 2)
d23_ln_agp_adj_L <- round(d23_ln_agp_adj_L, 2)
d23_ln_il1_adj_L <- round(d23_ln_il1_adj_L, 2)
d23_ln_il6_adj_L <- round(d23_ln_il6_adj_L, 2)
d23_ln_tnf_adj_L <- round(d23_ln_tnf_adj_L, 2)
d23_ln_il12_adj_L <- round(d23_ln_il12_adj_L, 2)
d23_ln_ifn_adj_L <- round(d23_ln_ifn_adj_L, 2)
d23_ln_il4_adj_L <- round(d23_ln_il4_adj_L, 2)
d23_ln_il5_adj_L <- round(d23_ln_il5_adj_L, 2)
d23_ln_il13_adj_L <- round(d23_ln_il13_adj_L, 2)
d23_ln_il17_adj_L <- round(d23_ln_il17_adj_L, 2)
d23_ln_il21_adj_L <- round(d23_ln_il21_adj_L, 2)
d23_ln_il10_adj_L <- round(d23_ln_il10_adj_L, 2)
d23_ln_il2_adj_L <- round(d23_ln_il2_adj_L, 2)
d23_ln_gmc_adj_L <- round(d23_ln_gmc_adj_L, 2)
d23_ln_igf_adj_L <- round(d23_ln_igf_adj_L, 2)

adjs6 <- c("95% CI", " ", " ", paste(d23_ln_il1_adj_L[1], " (", d23_ln_il1_adj_L[2], ", ", d23_ln_il1_adj_L[3], ")", sep=""),
             " ", " ", paste(d23_ln_il6_adj_L[1], " (", d23_ln_il6_adj_L[2], ", ", d23_ln_il6_adj_L[3], ")", sep=""),
           " ", " ", paste(d23_ln_tnf_adj_L[1], " (", d23_ln_tnf_adj_L[2], ", ", d23_ln_tnf_adj_L[3], ")", sep=""),
           " ", " ", paste(d23_ln_crp_adj_L[1], " (", d23_ln_crp_adj_L[2], ", ", d23_ln_crp_adj_L[3], ")", sep=""),
           " ", " ", paste(d23_ln_il12_adj_L[1], " (", d23_ln_il12_adj_L[2], ", ", d23_ln_il12_adj_L[3], ")", sep=""),
             " ", " ", paste(d23_ln_ifn_adj_L[1], " (", d23_ln_ifn_adj_L[2], ", ", d23_ln_ifn_adj_L[3], ")", sep=""),
             " ", " ", paste(d23_ln_il4_adj_L[1], " (", d23_ln_il4_adj_L[2], ", ", d23_ln_il4_adj_L[3], ")", sep=""),
             " ", " ", paste(d23_ln_il5_adj_L[1], " (", d23_ln_il5_adj_L[2], ", ", d23_ln_il5_adj_L[3], ")", sep=""),
             " ", " ", paste(d23_ln_il13_adj_L[1], " (", d23_ln_il13_adj_L[2], ", ", d23_ln_il13_adj_L[3], ")", sep=""),
             " ", " ", paste(d23_ln_il17_adj_L[1], " (", d23_ln_il17_adj_L[2], ", ", d23_ln_il17_adj_L[3], ")", sep=""),
             " ", " ", paste(d23_ln_il21_adj_L[1], " (", d23_ln_il21_adj_L[2], ", ", d23_ln_il21_adj_L[3], ")", sep=""),
             " ", " ", paste(d23_ln_il10_adj_L[1], " (", d23_ln_il10_adj_L[2], ", ", d23_ln_il10_adj_L[3], ")", sep=""),
             " ", " ", paste(d23_ln_il2_adj_L[1], " (", d23_ln_il2_adj_L[2], ", ", d23_ln_il2_adj_L[3], ")", sep=""),
           " ", " ", paste(d23_ln_gmc_adj_L[1], " (", d23_ln_gmc_adj_L[2], ", ", d23_ln_gmc_adj_L[3], ")", sep=""),
           " ", " ", paste(d23_ln_agp_adj_L[1], " (", d23_ln_agp_adj_L[2], ", ", d23_ln_agp_adj_L[3], ")", sep=""),
           " ", " ", paste(d23_ln_igf_adj_L[1], " (", d23_ln_igf_adj_L[2], ", ", d23_ln_igf_adj_L[3], ")", sep=""))

ipcws6<-c("95% CI", " ", " ", maketblvalue(d23_il1_adj_ipcw_L$`unlist(d23_il1_adj_ipcw$estimates$ATE)`),
          " ", " ", maketblvalue(d23_il6_adj_ipcw_L$`unlist(d23_il6_adj_ipcw$estimates$ATE)`),
          " ", " ", maketblvalue(d23_tnf_adj_ipcw_L$`unlist(d23_tnf_adj_ipcw$estimates$ATE)`),
          " ", " ", maketblvalue(d23_crp_adj_ipcw_L$`unlist(d23_crp_adj_ipcw_L$estimates$ATE)`),
          " ", " ", maketblvalue(d23_il12_adj_ipcw_L$`unlist(d23_il12_adj_ipcw$estimates$ATE)`),
          " ", " ", maketblvalue(d23_ifn_adj_ipcw_L$`unlist(d23_ifn_adj_ipcw$estimates$ATE)`),
          " ", " ", maketblvalue(d23_il4_adj_ipcw_L$`unlist(d23_il4_adj_ipcw$estimates$ATE)`),
          " ", " ", maketblvalue(d23_il5_adj_ipcw_L$`unlist(d23_il5_adj_ipcw$estimates$ATE)`),
          " ", " ", maketblvalue(d23_il13_adj_ipcw_L$`unlist(d23_il13_adj_ipcw$estimates$ATE)`),
          " ", " ", maketblvalue(d23_il17_adj_ipcw_L$`unlist(d23_il17_adj_ipcw$estimates$ATE)`),
          " ", " ", maketblvalue(d23_il21_adj_ipcw_L$`unlist(d23_il21_adj_ipcw$estimates$ATE)`),
          " ", " ", maketblvalue(d23_il10_adj_ipcw_L$`unlist(d23_il10_adj_ipcw$estimates$ATE)`),
          " ", " ", maketblvalue(d23_il2_adj_ipcw_L$`unlist(d23_il2_adj_ipcw$estimates$ATE)`),
          " ", " ", maketblvalue(d23_gmc_adj_ipcw_L$`unlist(d23_gmc_adj_ipcw$estimates$ATE)`),
          " ", " ", maketblvalue(d23_agp_adj_ipcw_L$`unlist(d23_agp_adj_ipcw_L$estimates$ATE)`),
          " ", " ", maketblvalue(d23_igf_adj_ipcw_L$`unlist(d23_igf_adj_ipcw$estimates$ATE)`))

rawipcwpvals6<-c("P-value", " ", " ", rawpval(d23_il1_adj_ipcw_L[5, 1]),
          " ", " ", rawpval(d23_il6_adj_ipcw_L[5, 1]),
          " ", " ", rawpval(d23_tnf_adj_ipcw_L[5, 1]),
          " ", " ", rawpval(d23_crp_adj_ipcw_L[5, 1]),
          " ", " ", rawpval(d23_il12_adj_ipcw_L[5, 1]),
          " ", " ", rawpval(d23_ifn_adj_ipcw_L[5, 1]),
          " ", " ", rawpval(d23_il4_adj_ipcw_L[5, 1]),
          " ", " ", rawpval(d23_il5_adj_ipcw_L[5, 1]),
          " ", " ", rawpval(d23_il13_adj_ipcw_L[5, 1]),
          " ", " ", rawpval(d23_il17_adj_ipcw_L[5, 1]),
          " ", " ", rawpval(d23_il21_adj_ipcw_L[5, 1]),
          " ", " ", rawpval(d23_il10_adj_ipcw_L[5, 1]),
          " ", " ", rawpval(d23_il2_adj_ipcw_L[5, 1]),
          " ", " ", rawpval(d23_gmc_adj_ipcw_L[5, 1]),
          " ", " ", rawpval(d23_agp_adj_ipcw_L[5, 1]),
          " ", " ", rawpval(d23_igf_adj_ipcw_L[5, 1]))

ipcwpvals6<-c("Bonferroni P-value", " ", " ", bonpval(d23_il1_adj_ipcw_L[5, 1]),
          " ", " ", bonpval(d23_il6_adj_ipcw_L[5, 1]),
          " ", " ", bonpval(d23_tnf_adj_ipcw_L[5, 1]),
          " ", " ", bonpval(d23_crp_adj_ipcw_L[5, 1]),
          " ", " ", bonpval(d23_il12_adj_ipcw_L[5, 1]),
          " ", " ", bonpval(d23_ifn_adj_ipcw_L[5, 1]),
          " ", " ", bonpval(d23_il4_adj_ipcw_L[5, 1]),
          " ", " ", bonpval(d23_il5_adj_ipcw_L[5, 1]),
          " ", " ", bonpval(d23_il13_adj_ipcw_L[5, 1]),
          " ", " ", bonpval(d23_il17_adj_ipcw_L[5, 1]),
          " ", " ", bonpval(d23_il21_adj_ipcw_L[5, 1]),
          " ", " ", bonpval(d23_il10_adj_ipcw_L[5, 1]),
          " ", " ", bonpval(d23_il2_adj_ipcw_L[5, 1]),
          " ", " ", bonpval(d23_gmc_adj_ipcw_L[5, 1]),
          " ", " ", bonpval(d23_agp_adj_ipcw_L[5, 1]),
          " ", " ", bonpval(d23_igf_adj_ipcw_L[5, 1]))

# Table S6: Effect of intervention on change in individual immune status and growth factor measurements between ages 14 and 28 months
tbls6 <- data.table(
  " " = outcomes6,
  " " = Ns6, 
  " " = absmeans6,
  " " = abssds6,
  " " = means6, 
  " " = sds6,
  "Unadjusted difference: Intervention vs. Control" = unadjs6,
  " " = rawpvals6,
  " " = pvals6,
  "Age- and sex- adjusted difference: Intervention vs. Control" = asadjs6, 
  " " = rawasapvals6,
  " " = asapvals6,
  "Fully adjusted difference: Intervention vs. Control" = adjs6,
  " " = rawadjpvals6,
  " " = adjpvals6,
  "IPCW adjusted difference: Intervention vs. Control" = ipcws6,
  " " = rawipcwpvals6,
  " " = ipcwpvals6
)

write.csv(tbls6, file=here('tables/supplementary/immune_supptable6.csv'))
print(xtable(tbls6), type="html", file=here("tables/supplementary/immune_supptable6.html"))



#### TABLE S7 ####
makeipcw<-function(var){
  rounded<-round(var, 2)
  paste(rounded[1], " (", rounded[3], ", ", rounded[4], ")", sep="")
}

asadjtbls7<-c("95% CI", " ", " ", makecival(d23_ratio_il1_il10_adj_sex_age_L),
                 " ", " ", makecival(d23_ratio_il6_il10_adj_sex_age_L),
                 " ", " ", makecival(d23_ratio_tnf_il10_adj_sex_age_L),
                 " ", " ", makecival(d23_ratio_il12_il10_adj_sex_age_L),
                 " ", " ", makecival(d23_ratio_ifn_il10_adj_sex_age_L),
                 " ", " ", makecival(d23_ratio_il4_il10_adj_sex_age_L),
                 " ", " ", makecival(d23_ratio_il5_il10_adj_sex_age_L),
                 " ", " ", makecival(d23_ratio_il13_il10_adj_sex_age_L),
                 " ", " ", makecival(d23_ratio_il17_il10_adj_sex_age_L),
                 " ", " ", makecival(d23_ratio_il21_il10_adj_sex_age_L),
                 " ", " ", makecival(d23_ratio_il2_il10_adj_sex_age_L),
                 " ", " ", makecival(d23_ratio_gmc_il10_adj_sex_age_L),
                 " ", " ", makecival(d23_ratio_il12_il4_adj_sex_age_L),
                 " ", " ", makecival(d23_ratio_ifn_il4_adj_sex_age_L),
                 " ", " ", makecival(d23_ratio_il12_il5_adj_sex_age_L),
                 " ", " ", makecival(d23_ratio_ifn_il5_adj_sex_age_L),
                 " ", " ", makecival(d23_ratio_il12_il13_adj_sex_age_L),
                 " ", " ", makecival(d23_ratio_ifn_il13_adj_sex_age_L),
                 " ", " ", makecival(d23_ratio_il12_il17_adj_sex_age_L),
                 " ", " ", makecival(d23_ratio_ifn_il17_adj_sex_age_L),
                 " ", " ", makecival(d23_ratio_il12_il21_adj_sex_age_L),
                 " ", " ", makecival(d23_ratio_ifn_il21_adj_sex_age_L),
                 " ", " ", makecival(d23_ratio_pro_il10_adj_sex_age_L),
                 " ", " ", makecival(d23_ratio_th1_il10_adj_sex_age_L),
                 " ", " ", makecival(d23_ratio_th2_il10_adj_sex_age_L),
                 " ", " ", makecival(d23_ratio_th17_il10_adj_sex_age_L),
                 " ", " ", makecival(d23_ratio_th1_th2_adj_sex_age_L),
                 " ", " ", makecival(d23_ratio_th1_th17_adj_sex_age_L))

rawpvals7<-c("P-value", " ", " ", rawpval(d23_ratio_il1_il10_adj_sex_age_L[6]),
             " ", " ", rawpval(d23_ratio_il6_il10_adj_sex_age_L[6]), " ", " ", rawpval(d23_ratio_tnf_il10_adj_sex_age_L[6]), " ", " ", rawpval(d23_ratio_il12_il10_adj_sex_age_L[6]),
             " ", " ", rawpval(d23_ratio_ifn_il10_adj_sex_age_L[6]), " ", " ", rawpval(d23_ratio_il4_il10_adj_sex_age_L[6]), " ", " ", rawpval(d23_ratio_il5_il10_adj_sex_age_L[6]),
             " ", " ", rawpval(d23_ratio_il13_il10_adj_sex_age_L[6]), " ", " ", rawpval(d23_ratio_il17_il10_adj_sex_age_L[6]), " ", " ", rawpval(d23_ratio_il21_il10_adj_sex_age_L[6]),
             " ", " ", rawpval(d23_ratio_il2_il10_adj_sex_age_L[6]), " ", " ", rawpval(d23_ratio_gmc_il10_adj_sex_age_L[6]), " ", " ", rawpval(d23_ratio_il12_il4_adj_sex_age_L[6]),
             " ", " ", rawpval(d23_ratio_ifn_il4_adj_sex_age_L[6]), " ", " ", rawpval(d23_ratio_il12_il5_adj_sex_age_L[6]), " ", " ", rawpval(d23_ratio_ifn_il5_adj_sex_age_L[6]),
             " ", " ", rawpval(d23_ratio_il12_il13_adj_sex_age_L[6]), " ", " ", rawpval(d23_ratio_ifn_il13_adj_sex_age_L[6]), " ", " ", rawpval(d23_ratio_il12_il17_adj_sex_age_L[6]),
             " ", " ", rawpval(d23_ratio_ifn_il17_adj_sex_age_L[6]), " ", " ", rawpval(d23_ratio_il12_il21_adj_sex_age_L[6]), " ", " ", rawpval(d23_ratio_ifn_il21_adj_sex_age_L[6]),
             " ", " ", rawpval(d23_ratio_pro_il10_adj_sex_age_L[6]), " ", " ", rawpval(d23_ratio_th1_il10_adj_sex_age_L[6]), " ", " ", rawpval(d23_ratio_th2_il10_adj_sex_age_L[6]),
             " ", " ", rawpval(d23_ratio_th17_il10_adj_sex_age_L[6]), " ", " ", rawpval(d23_ratio_th1_th2_adj_sex_age_L[6]), " ", " ", rawpval(d23_ratio_th1_th17_adj_sex_age_L[6]))

pvals7<-c("Bonferroni P-value", " ", " ", bonpval(d23_ratio_il1_il10_adj_sex_age_L[6]),
             " ", " ", bonpval(d23_ratio_il6_il10_adj_sex_age_L[6]), " ", " ", bonpval(d23_ratio_tnf_il10_adj_sex_age_L[6]), " ", " ", bonpval(d23_ratio_il12_il10_adj_sex_age_L[6]),
             " ", " ", bonpval(d23_ratio_ifn_il10_adj_sex_age_L[6]), " ", " ", bonpval(d23_ratio_il4_il10_adj_sex_age_L[6]), " ", " ", bonpval(d23_ratio_il5_il10_adj_sex_age_L[6]),
             " ", " ", bonpval(d23_ratio_il13_il10_adj_sex_age_L[6]), " ", " ", bonpval(d23_ratio_il17_il10_adj_sex_age_L[6]), " ", " ", bonpval(d23_ratio_il21_il10_adj_sex_age_L[6]),
             " ", " ", bonpval(d23_ratio_il2_il10_adj_sex_age_L[6]), " ", " ", bonpval(d23_ratio_gmc_il10_adj_sex_age_L[6]), " ", " ", bonpval(d23_ratio_il12_il4_adj_sex_age_L[6]),
             " ", " ", bonpval(d23_ratio_ifn_il4_adj_sex_age_L[6]), " ", " ", bonpval(d23_ratio_il12_il5_adj_sex_age_L[6]), " ", " ", bonpval(d23_ratio_ifn_il5_adj_sex_age_L[6]),
             " ", " ", bonpval(d23_ratio_il12_il13_adj_sex_age_L[6]), " ", " ", bonpval(d23_ratio_ifn_il13_adj_sex_age_L[6]), " ", " ", bonpval(d23_ratio_il12_il17_adj_sex_age_L[6]),
             " ", " ", bonpval(d23_ratio_ifn_il17_adj_sex_age_L[6]), " ", " ", bonpval(d23_ratio_il12_il21_adj_sex_age_L[6]), " ", " ", bonpval(d23_ratio_ifn_il21_adj_sex_age_L[6]),
             " ", " ", bonpval(d23_ratio_pro_il10_adj_sex_age_L[6]), " ", " ", bonpval(d23_ratio_th1_il10_adj_sex_age_L[6]), " ", " ", bonpval(d23_ratio_th2_il10_adj_sex_age_L[6]),
             " ", " ", bonpval(d23_ratio_th17_il10_adj_sex_age_L[6]), " ", " ", bonpval(d23_ratio_th1_th2_adj_sex_age_L[6]), " ", " ", bonpval(d23_ratio_th1_th17_adj_sex_age_L[6]))

ipcws7 <- c("95% CI", " ", " ", makeipcw(d23_ratio_il1_il10_adj_ipcw_L$`unlist(d23_ratio_il1_il10_adj_ipcw$estimates$ATE)`),                                      " ", " ", makeipcw(d23_ratio_il6_il10_adj_ipcw_L$`unlist(d23_ratio_il6_il10_adj_ipcw$estimates$ATE)`),
            " ", " ", makeipcw(d23_ratio_tnf_il10_adj_ipcw_L$`unlist(d23_ratio_tnf_il10_adj_ipcw$estimates$ATE)`),
            " ", " ", makeipcw(d23_ratio_il12_il10_adj_ipcw_L$`unlist(d23_ratio_il12_il10_adj_ipcw$estimates$ATE)`),
            " ", " ", makeipcw(d23_ratio_ifn_il10_adj_ipcw_L$`unlist(d23_ratio_ifn_il10_adj_ipcw$estimates$ATE)`),
            " ", " ", makeipcw(d23_ratio_il4_il10_adj_ipcw_L$`unlist(d23_ratio_il4_il10_adj_ipcw$estimates$ATE)`), 
            " ", " ", makeipcw(d23_ratio_il5_il10_adj_ipcw_L$`unlist(d23_ratio_il5_il10_adj_ipcw$estimates$ATE)`),
            " ", " ", makeipcw(d23_ratio_il13_il10_adj_ipcw_L$`unlist(d23_ratio_il13_il10_adj_ipcw$estimates$ATE)`),
            " ", " ", makeipcw(d23_ratio_il17_il10_adj_ipcw_L$`unlist(d23_ratio_il17_il10_adj_ipcw$estimates$ATE)`),
            " ", " ", makeipcw(d23_ratio_il21_il10_adj_ipcw_L$`unlist(d23_ratio_il21_il10_adj_ipcw$estimates$ATE)`),
            " ", " ", makeipcw(d23_ratio_il2_il10_adj_ipcw_L$`unlist(d23_ratio_il2_il10_adj_ipcw$estimates$ATE)`),
            " ", " ", makeipcw(d23_ratio_gmc_il10_adj_ipcw_L$`unlist(d23_ratio_gmc_il10_adj_ipcw$estimates$ATE)`),
            " ", " ", makeipcw(d23_ratio_il12_il4_adj_ipcw_L$`unlist(d23_ratio_il12_il4_adj_ipcw$estimates$ATE)`),
            " ", " ", makeipcw(d23_ratio_ifn_il4_adj_ipcw_L$`unlist(d23_ratio_ifn_il4_adj_ipcw$estimates$ATE)`),
            " ", " ", makeipcw(d23_ratio_il12_il5_adj_ipcw_L$`unlist(d23_ratio_il12_il5_adj_ipcw$estimates$ATE)`),
            " ", " ", makeipcw(d23_ratio_ifn_il5_adj_ipcw_L$`unlist(d23_ratio_ifn_il5_adj_ipcw$estimates$ATE)`),
            " ", " ", makeipcw(d23_ratio_il12_il13_adj_ipcw_L$`unlist(d23_ratio_il12_il13_adj_ipcw$estimates$ATE)`),
            " ", " ", makeipcw(d23_ratio_ifn_il13_adj_ipcw_L$`unlist(d23_ratio_ifn_il13_adj_ipcw$estimates$ATE)`),
            " ", " ", makeipcw(d23_ratio_il12_il17_adj_ipcw_L$`unlist(d23_ratio_il12_il17_adj_ipcw$estimates$ATE)`),
            " ", " ", makeipcw(d23_ratio_ifn_il17_adj_ipcw_L$`unlist(d23_ratio_ifn_il17_adj_ipcw$estimates$ATE)`),
            " ", " ", makeipcw(d23_ratio_il12_il21_adj_ipcw_L$`unlist(d23_ratio_il12_il21_adj_ipcw$estimates$ATE)`),
            " ", " ", makeipcw(d23_ratio_ifn_il21_adj_ipcw_L$`unlist(d23_ratio_ifn_il21_adj_ipcw$estimates$ATE)`),
            " ", " ", makeipcw(d23_ratio_pro_il10_adj_ipcw_L$`unlist(d23_ratio_pro_il10_adj_ipcw$estimates$ATE)`),
            " ", " ", makeipcw(d23_ratio_th1_il10_adj_ipcw_L$`unlist(d23_ratio_th1_il10_adj_ipcw$estimates$ATE)`),
            " ", " ", makeipcw(d23_ratio_th2_il10_adj_ipcw_L$`unlist(d23_ratio_th2_il10_adj_ipcw$estimates$ATE)`),
            " ", " ", makeipcw(d23_ratio_th17_il10_adj_ipcw_L$`unlist(d23_ratio_th17_il10_adj_ipcw$estimates$ATE)`),
            " ", " ", makeipcw(d23_ratio_th1_th2_adj_ipcw_L$`unlist(d23_ratio_th1_th2_adj_ipcw$estimates$ATE)`),
            " ", " ", makeipcw(d23_ratio_th1_th17_adj_ipcw_L$`unlist(d23_ratio_th1_th17_adj_ipcw$estimates$ATE)`))

rawipcwpvals7 <- c("P-value", " ", " ", rawpval(d23_ratio_il1_il10_adj_ipcw_L[5,1]),                                      
                " ", " ", rawpval(d23_ratio_il6_il10_adj_ipcw_L[5,1]),
            " ", " ", rawpval(d23_ratio_tnf_il10_adj_ipcw_L[5,1]),
            " ", " ", rawpval(d23_ratio_il12_il10_adj_ipcw_L[5,1]),
            " ", " ", rawpval(d23_ratio_ifn_il10_adj_ipcw_L[5,1]),
            " ", " ", rawpval(d23_ratio_il4_il10_adj_ipcw_L[5,1]), 
            " ", " ", rawpval(d23_ratio_il5_il10_adj_ipcw_L[5,1]),
            " ", " ", rawpval(d23_ratio_il13_il10_adj_ipcw_L[5,1]),
            " ", " ", rawpval(d23_ratio_il17_il10_adj_ipcw_L[5,1]),
            " ", " ", rawpval(d23_ratio_il21_il10_adj_ipcw_L[5,1]),
            " ", " ", rawpval(d23_ratio_il2_il10_adj_ipcw_L[5,1]),
            " ", " ", rawpval(d23_ratio_gmc_il10_adj_ipcw_L[5,1]),
            " ", " ", rawpval(d23_ratio_il12_il4_adj_ipcw_L[5,1]),
            " ", " ", rawpval(d23_ratio_ifn_il4_adj_ipcw_L[5,1]),
            " ", " ", rawpval(d23_ratio_il12_il5_adj_ipcw_L[5,1]),
            " ", " ", rawpval(d23_ratio_ifn_il5_adj_ipcw_L[5,1]),
            " ", " ", rawpval(d23_ratio_il12_il13_adj_ipcw_L[5,1]),
            " ", " ", rawpval(d23_ratio_ifn_il13_adj_ipcw_L[5,1]),
            " ", " ", rawpval(d23_ratio_il12_il17_adj_ipcw_L[5,1]),
            " ", " ", rawpval(d23_ratio_ifn_il17_adj_ipcw_L[5,1]),
            " ", " ", rawpval(d23_ratio_il12_il21_adj_ipcw_L[5,1]),
            " ", " ", rawpval(d23_ratio_ifn_il21_adj_ipcw_L[5,1]),
            " ", " ", rawpval(d23_ratio_pro_il10_adj_ipcw_L[5,1]),
            " ", " ", rawpval(d23_ratio_th1_il10_adj_ipcw_L[5,1]),
            " ", " ", rawpval(d23_ratio_th2_il10_adj_ipcw_L[5,1]),
            " ", " ", rawpval(d23_ratio_th17_il10_adj_ipcw_L[5,1]),
            " ", " ", rawpval(d23_ratio_th1_th2_adj_ipcw_L[5,1]),
            " ", " ", rawpval(d23_ratio_th1_th17_adj_ipcw_L[5,1]))

ipcwpvals7 <- c("Bonferroni P-value", " ", " ", bonpval(d23_ratio_il1_il10_adj_ipcw_L[5,1]),                                      
                " ", " ", bonpval(d23_ratio_il6_il10_adj_ipcw_L[5,1]),
            " ", " ", bonpval(d23_ratio_tnf_il10_adj_ipcw_L[5,1]),
            " ", " ", bonpval(d23_ratio_il12_il10_adj_ipcw_L[5,1]),
            " ", " ", bonpval(d23_ratio_ifn_il10_adj_ipcw_L[5,1]),
            " ", " ", bonpval(d23_ratio_il4_il10_adj_ipcw_L[5,1]), 
            " ", " ", bonpval(d23_ratio_il5_il10_adj_ipcw_L[5,1]),
            " ", " ", bonpval(d23_ratio_il13_il10_adj_ipcw_L[5,1]),
            " ", " ", bonpval(d23_ratio_il17_il10_adj_ipcw_L[5,1]),
            " ", " ", bonpval(d23_ratio_il21_il10_adj_ipcw_L[5,1]),
            " ", " ", bonpval(d23_ratio_il2_il10_adj_ipcw_L[5,1]),
            " ", " ", bonpval(d23_ratio_gmc_il10_adj_ipcw_L[5,1]),
            " ", " ", bonpval(d23_ratio_il12_il4_adj_ipcw_L[5,1]),
            " ", " ", bonpval(d23_ratio_ifn_il4_adj_ipcw_L[5,1]),
            " ", " ", bonpval(d23_ratio_il12_il5_adj_ipcw_L[5,1]),
            " ", " ", bonpval(d23_ratio_ifn_il5_adj_ipcw_L[5,1]),
            " ", " ", bonpval(d23_ratio_il12_il13_adj_ipcw_L[5,1]),
            " ", " ", bonpval(d23_ratio_ifn_il13_adj_ipcw_L[5,1]),
            " ", " ", bonpval(d23_ratio_il12_il17_adj_ipcw_L[5,1]),
            " ", " ", bonpval(d23_ratio_ifn_il17_adj_ipcw_L[5,1]),
            " ", " ", bonpval(d23_ratio_il12_il21_adj_ipcw_L[5,1]),
            " ", " ", bonpval(d23_ratio_ifn_il21_adj_ipcw_L[5,1]),
            " ", " ", bonpval(d23_ratio_pro_il10_adj_ipcw_L[5,1]),
            " ", " ", bonpval(d23_ratio_th1_il10_adj_ipcw_L[5,1]),
            " ", " ", bonpval(d23_ratio_th2_il10_adj_ipcw_L[5,1]),
            " ", " ", bonpval(d23_ratio_th17_il10_adj_ipcw_L[5,1]),
            " ", " ", bonpval(d23_ratio_th1_th2_adj_ipcw_L[5,1]),
            " ", " ", bonpval(d23_ratio_th1_th17_adj_ipcw_L[5,1]))

tbls7<-data.table(
  " " = tbl4[[1]],
  " " = tbl4[[2]],
  " " = tbl4[[3]],
  " " = tbl4[[4]],
  " " = tbl4[[5]],
  " " = tbl4[[6]],
  "Unadjusted difference: Intervention vs. Control" = tbl4[[7]],
  " " = rawpval4,
  " " = as_bonf_col(tbl4[[8]]),
  "Age- and sex- adjusted difference: Intervention vs. Control" = asadjtbls7,
  " " = rawpvals7,
  " " = pvals7,
  "Fully adjusted difference: Intervention vs. Control" = tbl4[[9]],
  " " = rawadjpval4,
  " " = as_bonf_col(tbl4[[10]]),
  "IPCW adjusted difference: Intervention vs. Control" = ipcws7,
  " " = rawipcwpvals7,
  " " = ipcwpvals7)

write.csv(tbls7, file=here('tables/supplementary/immune_supptable7.csv'))
print(xtable(tbls7), type="html", file=here("tables/supplementary/immune_supptable7.html"))




## Machine-readable copy for the Quarto supplement (manuscript/supplement.qmd).
## The CSV/HTML above are kept for backwards compatibility.
saveRDS(list(S1 = tbls1, S2 = tbls2, S3 = tbls3, S4 = tbls4,
             S5 = tbls5, S6 = tbls6, S7 = tbls7),
        here("results/supplement_tables_1_7.rds"))
