rm(list=ls())
library("xtable")
source(here::here("0-config.R"))

source(here::here("table scripts/immune-main-tables.R"))
load(here::here("results/immune_ipcw.RData"))

setwd(paste0(dropboxDir,"Data/Cleaned/Audrie/"))
lab<-read.csv("bangladesh-dm-ee-anthro-diar-ee-med-plasma-blind-tr-enrol-covariates-lab.csv", stringsAsFactors = TRUE)

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

outcomes2 <- c(paste("Ln IL-1", "Î²", " (pg/ml)", sep=""), "Control", "Nutrition + WSH", 
             "Ln IL-6 (pg/ml)", "Control", "Nutrition + WSH", 
             paste("Ln TNF-", "Î±", " (pg/ml)", sep=""), "Control", "Nutrition + WSH",
             "Ln CRP (mg/L)", "Control", "Nutrition + WSH", 
             "Ln IL-12 (pg/ml)", "Control", "Nutrition + WSH", 
             paste("Ln IFN-", "Î³", " (pg/ml)", sep=""), "Control", "Nutrition + WSH", 
             "Ln IL-4 (pg/ml)", "Control", "Nutrition + WSH", 
             "Ln IL-5 (pg/ml)", "Control", "Nutrition + WSH", 
             "Ln IL-13 (pg/ml)", "Control", "Nutrition + WSH", 
             "Ln IL-17A (pg/ml)", "Control", "Nutrition + WSH", 
             "Ln IL-21 (pg/ml)", "Control", "Nutrition + WSH", 
             "Ln IL-10 (pg/ml)", "Control", "Nutrition + WSH", 
             "Ln IL-2 (pg/ml)", "Control", "Nutrition + WSH", 
             "Ln GM-CSF (pg/ml)", "Control", "Nutrition + WSH", 
             "Ln AGP (g/L)", "Control", "Nutrition + WSH",
             paste("Ln IGF-1 (", "Î¼", "g/L)", sep=""), "Control", "Nutrition + WSH")

Ns2 <- c(" ", as.character(il1_t2_N_tr$t2_ln_il1_N_tr[1]), as.character(il1_t2_N_tr$t2_ln_il1_N_tr[2]),
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

abssds2 <-c(" ", as.character(round(abs_il1_t2_N_tr$sd[1], 2)), as.character(round(abs_il1_t2_N_tr$sd[2], 2)),
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

absmeans2 <- c(" ", as.character(round(abs_il1_t2_N_tr$mean[1], 2)), as.character(round(abs_il1_t2_N_tr$mean[2], 2)),
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

means2 <- c(" ", as.character(round(il1_t2_N_tr$mean[1], 2)), as.character(round(il1_t2_N_tr$mean[2], 2)),
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

sds2 <- c(" ", as.character(round(il1_t2_N_tr$sd[1], 2)), as.character(round(il1_t2_N_tr$sd[2], 2)),
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

unadjs2 <- c(" ", " ", paste(t2_il1_unadj_L[1], " (", t2_il1_unadj_L[2], ", ", t2_il1_unadj_L[3], ")", sep=""),
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

asadjs2 <- c(" ", " ", paste(t2_il1_adj_sex_age_L[1], " (", t2_il1_adj_sex_age_L[2], ", ", t2_il1_adj_sex_age_L[3], ")", sep=""),
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

adjs2 <- c(" ", " ", paste(t2_il1_adj_L[1], " (", t2_il1_adj_L[2], ", ", t2_il1_adj_L[3], ")", sep=""),
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

ipcws2 <-c(" ", " ", maketblvalue(il12_t2_adj_ipcw_L$`unlist(il12_t2_adj_ipcw$estimates$ATE)`),
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

# Table S2: Effect of intervention on individual immune status and growth factor measurements at age 14 months
tbls2 <- data.table(
  "Outcome, Arm" = outcomes2,
  "N" = Ns2, 
  "Absolute Mean" = absmeans2,
  "Absolute SD" = abssds2,
  "Mean" = means2, 
  "SD" = sds2,
  "Unadjusted difference: Intervention vs. Control (95% CI)" = unadjs2,
  "Age- and sex- adjusted difference: Intervention vs. Control (95% CI)" = asadjs2, 
  "Fully adjusted difference: Intervention vs. Control (95% CI)" = adjs2,
  "IPCW adjusted difference: Intervention vs. Control (95% CI)" = ipcws2)

write.csv(tbls2, file=here('tables/supplementary/immune_supptable2.csv'))
print(xtable(tbls2), type="html", file=here("tables/supplementary/immune_supptable2.html"))



#### TABLE S3 ####

ipcws3<-c(" ", " ", maketblvalue(ratio_il1_il10_t2_adj_ipcw_L$`unlist(ratio_il1_il10_t2_adj_ipcw$estimates$ATE)`),
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

tbls3<-cbind(tbl2, ipcws3)
names(tbls3)[10]<-"IPCW adjusted difference: Intervention vs. Control (95% CI)"

write.csv(tbls3, file=here('tables/supplementary/immune_supptable3.csv'))
print(xtable(tbls3), type="html", file=here("tables/supplementary/immune_supptable3.html"))

#### TABLE S4 ####

outcomes4 <- c(paste("Ln IL-1", "Î²", " (pg/ml)", sep=""), "Control", "Nutrition + WSH", 
             "Ln IL-6 (pg/ml)", "Control", "Nutrition + WSH", 
             paste("Ln TNF-", "Î±", " (pg/ml)", sep=""), "Control", "Nutrition + WSH",
             "Ln IL-12 (pg/ml)", "Control", "Nutrition + WSH", 
             paste("Ln IFN-", "Î³", " (pg/ml)", sep=""), "Control", "Nutrition + WSH", 
             "Ln IL-4 (pg/ml)", "Control", "Nutrition + WSH", 
             "Ln IL-5 (pg/ml)", "Control", "Nutrition + WSH", 
             "Ln IL-13 (pg/ml)", "Control", "Nutrition + WSH", 
             "Ln IL-17A (pg/ml)", "Control", "Nutrition + WSH", 
             "Ln IL-21 (pg/ml)", "Control", "Nutrition + WSH", 
             "Ln IL-10 (pg/ml)", "Control", "Nutrition + WSH", 
             "Ln IL-2 (pg/ml)", "Control", "Nutrition + WSH", 
             "Ln GM-CSF (pg/ml)", "Control", "Nutrition + WSH", 
             paste("Ln IGF-1 (", "Î¼", "g/L)", sep=""), "Control", "Nutrition + WSH")

Ns4 <- c(" ", as.character(il1_t3_N_tr$t3_ln_il1_N_tr[1]), as.character(il1_t3_N_tr$t3_ln_il1_N_tr[2]),
       " ", as.character(il6_t3_N_tr$t3_ln_il6_N_tr[1]), as.character(il6_t3_N_tr$t3_ln_il6_N_tr[2]),
       " ", as.character(tnf_t3_N_tr$t3_ln_tnf_N_tr[1]), as.character(tnf_t3_N_tr$t3_ln_tnf_N_tr[2]),
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
       " ", as.character(igf_t3_N_tr$t3_ln_igf_N_tr[1]), as.character(igf_t3_N_tr$t3_ln_igf_N_tr[2]))

absmeans4 <- c(" ", as.character(round(abs_il1_t3_N_tr$mean[1], 2)), as.character(round(abs_il1_t3_N_tr$mean[2], 2)),
             " ", as.character(round(abs_il6_t3_N_tr$mean[1], 2)), as.character(round(abs_il6_t3_N_tr$mean[2], 2)),
             " ", as.character(round(abs_tnf_t3_N_tr$mean[1], 2)), as.character(round(abs_tnf_t3_N_tr$mean[2], 2)),
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
             " ", as.character(round(abs_igf_t3_N_tr$mean[1], 2)), as.character(round(abs_igf_t3_N_tr$mean[2], 2)))

abssds4 <- c(" ", as.character(round(abs_il1_t3_N_tr$sd[1], 2)), as.character(round(abs_il1_t3_N_tr$sd[2], 2)),
           " ", as.character(round(abs_il6_t3_N_tr$sd[1], 2)), as.character(round(abs_il6_t3_N_tr$sd[2], 2)),
           " ", as.character(round(abs_tnf_t3_N_tr$sd[1], 2)), as.character(round(abs_tnf_t3_N_tr$sd[2], 2)),
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
           " ", as.character(round(abs_igf_t3_N_tr$sd[1], 2)), as.character(round(abs_igf_t3_N_tr$sd[2], 2)))

means4 <- c(" ", as.character(round(il1_t3_N_tr$mean[1], 2)), as.character(round(il1_t3_N_tr$mean[2], 2)),
          " ", as.character(round(il6_t3_N_tr$mean[1], 2)), as.character(round(il6_t3_N_tr$mean[2], 2)),
          " ", as.character(round(tnf_t3_N_tr$mean[1], 2)), as.character(round(tnf_t3_N_tr$mean[2], 2)),
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
          " ", as.character(round(igf_t3_N_tr$mean[1], 2)), as.character(round(igf_t3_N_tr$mean[2], 2)))

sds4 <- c(" ", as.character(round(il1_t3_N_tr$sd[1], 2)), as.character(round(il1_t3_N_tr$sd[2], 2)),
        " ", as.character(round(il6_t3_N_tr$sd[1], 2)), as.character(round(il6_t3_N_tr$sd[2], 2)),
        " ", as.character(round(tnf_t3_N_tr$sd[1], 2)), as.character(round(tnf_t3_N_tr$sd[2], 2)),
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
        " ", as.character(round(igf_t3_N_tr$sd[1], 2)), as.character(round(igf_t3_N_tr$sd[2], 2)))

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

unadjs4 <- c(" ", " ", paste(t3_il1_unadj_L[1], " (", t3_il1_unadj_L[2], ", ", t3_il1_unadj_L[3], ")", sep=""),
           " ", " ", paste(t3_il6_unadj_L[1], " (", t3_il6_unadj_L[2], ", ", t3_il6_unadj_L[3], ")", sep=""),
           " ", " ", paste(t3_tnf_unadj_L[1], " (", t3_tnf_unadj_L[2], ", ", t3_tnf_unadj_L[3], ")", sep=""),
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
           " ", " ", paste(t3_igf_unadj_L[1], " (", t3_igf_unadj_L[2], ", ", t3_igf_unadj_L[3], ")", sep=""))

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

asadjs4 <- c(" ", " ", paste(t3_il1_adj_sex_age_L[1], " (", t3_il1_adj_sex_age_L[2], ", ", t3_il1_adj_sex_age_L[3], ")", sep=""),
           " ", " ", paste(t3_il6_adj_sex_age_L[1], " (", t3_il6_adj_sex_age_L[2], ", ", t3_il6_adj_sex_age_L[3], ")", sep=""),
           " ", " ", paste(t3_tnf_adj_sex_age_L[1], " (", t3_tnf_adj_sex_age_L[2], ", ", t3_tnf_adj_sex_age_L[3], ")", sep=""),
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
           " ", " ", paste(t3_igf_adj_sex_age_L[1], " (", t3_igf_adj_sex_age_L[2], ", ", t3_igf_adj_sex_age_L[3], ")", sep=""))

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

adjs4 <- c(" ", " ", paste(t3_il1_adj_L[1], " (", t3_il1_adj_L[2], ", ", t3_il1_adj_L[3], ")", sep=""),
         " ", " ", paste(t3_il6_adj_L[1], " (", t3_il6_adj_L[2], ", ", t3_il6_adj_L[3], ")", sep=""),
         " ", " ", paste(t3_tnf_adj_L[1], " (", t3_tnf_adj_L[2], ", ", t3_tnf_adj_L[3], ")", sep=""),
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
         " ", " ", paste(t3_igf_adj_L[1], " (", t3_igf_adj_L[2], ", ", t3_igf_adj_L[3], ")", sep=""))

ipcws4<-c(" ", " ", maketblvalue(il1_t3_adj_ipcw_L$`unlist(il1_t3_adj_ipcw$estimates$ATE)`),
          " ", " ", maketblvalue(il6_t3_adj_ipcw_L$`unlist(il6_t3_adj_ipcw$estimates$ATE)`),
          " ", " ", maketblvalue(tnf_t3_adj_ipcw_L$`unlist(tnf_t3_adj_ipcw$estimates$ATE)`),
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
          " ", " ", maketblvalue(igf_t3_adj_ipcw_L$`unlist(igf_t3_adj_ipcw$estimates$ATE)`))


# Table S4: Effect of intervention on individual immune status and growth factor measurements at age 28 months
tbls4 <- data.table(
  "Outcome, Arm" = outcomes4,
  "N" = Ns4, 
  "Absolute Mean" = absmeans4,
  "Absolute SD" = abssds4,
  "Mean" = means4, 
  "SD" = sds4,
  "Unadjusted difference: Intervention vs. Control (95% CI)" = unadjs4,
  "Age- and sex- adjusted difference: Intervention vs. Control (95% CI)" = asadjs4, 
  "Fully adjusted difference: Intervention vs. Control (95% CI)" = adjs4,
  "IPCW adjusted difference: Intervention vs. Control (95% CI)" = ipcws4
)

write.csv(tbls4, file=here('tables/supplementary/immune_supptable4.csv'))
print(xtable(tbls4), type="html", file=here("tables/supplementary/immune_supptable4.html"))


#### TABLE S5 ####
ipcws5<-c(" ", " ", maketblvalue(ratio_il1_il10_t3_adj_ipcw_L$`unlist(ratio_il1_il10_t3_adj_ipcw$estimates$ATE)`),
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

tbls5<-cbind(tbl3, ipcws5)
names(tbls5)[10]<-"IPCW adjusted difference: Intervention vs. Control (95% CI)"

write.csv(tbls5, file=here('tables/supplementary/immune_supptable5.csv'))
print(xtable(tbls5), type="html", file=here("tables/supplementary/immune_supptable5.html"))



#### TABLE S6 ####

outcomes6 <- c(paste("Ln ΔIL-1", "Î²", " (pg/ml)", sep=""), "Control", "Nutrition + WSH", 
                 "Ln ΔIL-6 (pg/ml)", "Control", "Nutrition + WSH", 
                 paste("Ln ΔTNF-", "Î±", " (pg/ml)", sep=""), "Control", "Nutrition + WSH",
                 "Ln ΔIL-12 (pg/ml)", "Control", "Nutrition + WSH", 
                 paste("Ln ΔIFN-", "Î³", " (pg/ml)", sep=""), "Control", "Nutrition + WSH", 
                 "Ln ΔIL-4 (pg/ml)", "Control", "Nutrition + WSH", 
                 "Ln ΔIL-5 (pg/ml)", "Control", "Nutrition + WSH", 
                 "Ln ΔIL-13 (pg/ml)", "Control", "Nutrition + WSH", 
                 "Ln ΔIL-17A (pg/ml)", "Control", "Nutrition + WSH", 
                 "Ln ΔIL-21 (pg/ml)", "Control", "Nutrition + WSH", 
                 "Ln ΔIL-10 (pg/ml)", "Control", "Nutrition + WSH", 
                 "Ln ΔIL-2 (pg/ml)", "Control", "Nutrition + WSH", 
                 "Ln ΔGM-CSF (pg/ml)", "Control", "Nutrition + WSH", 
                 paste("Ln ΔIGF-1 (", "Î¼", "g/L)", sep=""), "Control", "Nutrition + WSH")

Ns6 <- c(" ", as.character(d23_ln_il1_N_tr$d23_ln_il1_N_tr[1]), as.character(d23_ln_il1_N_tr$d23_ln_il1_N_tr[2]), 
           " ", as.character(d23_ln_il6_N_tr$d23_ln_il6_N_tr[1]), as.character(d23_ln_il6_N_tr$d23_ln_il6_N_tr[2]),
           " ", as.character(d23_ln_tnf_N_tr$d23_ln_tnf_N_tr[1]), as.character(d23_ln_tnf_N_tr$d23_ln_tnf_N_tr[2]),
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
           " ", as.character(d23_ln_igf_N_tr$d23_ln_igf_N_tr[1]), as.character(d23_ln_igf_N_tr$d23_ln_igf_N_tr[2]))

absmeans6 <- c(" ", as.character(round(abs_d23_il1_N_tr$mean[1], 2)), as.character(round(abs_d23_il1_N_tr$mean[2], 2)), 
                 " ", as.character(round(abs_d23_il6_N_tr$mean[1], 2)), as.character(round(abs_d23_il6_N_tr$mean[2], 2)),
                 " ", as.character(round(abs_d23_tnf_N_tr$mean[1], 2)), as.character(round(abs_d23_tnf_N_tr$mean[2], 2)),
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
                 " ", as.character(round(abs_d23_igf_N_tr$mean[1], 2)), as.character(round(abs_d23_igf_N_tr$mean[2], 2)))

abssds6 <- c(" ", as.character(round(abs_d23_il1_N_tr$sd[1], 2)), as.character(round(abs_d23_il1_N_tr$sd[2], 2)), 
               " ", as.character(round(abs_d23_il6_N_tr$sd[1], 2)), as.character(round(abs_d23_il6_N_tr$sd[2], 2)),
               " ", as.character(round(abs_d23_tnf_N_tr$sd[1], 2)), as.character(round(abs_d23_tnf_N_tr$sd[2], 2)),
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
               " ", as.character(round(abs_d23_igf_N_tr$sd[1], 2)), as.character(round(abs_d23_igf_N_tr$sd[2], 2)))

means6 <- c(" ", as.character(round(d23_ln_il1_N_tr$mean[1], 2)), as.character(round(d23_ln_il1_N_tr$mean[2], 2)), 
              " ", as.character(round(d23_ln_il6_N_tr$mean[1], 2)), as.character(round(d23_ln_il6_N_tr$mean[2], 2)),
              " ", as.character(round(d23_ln_tnf_N_tr$mean[1], 2)), as.character(round(d23_ln_tnf_N_tr$mean[2], 2)),
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
              " ", as.character(round(d23_ln_igf_N_tr$mean[1], 2)), as.character(round(d23_ln_igf_N_tr$mean[2], 2)))

sds6 <- c(" ", as.character(round(d23_ln_il1_N_tr$sd[1], 2)), as.character(round(d23_ln_il1_N_tr$sd[2], 2)), 
            " ", as.character(round(d23_ln_il6_N_tr$sd[1], 2)), as.character(round(d23_ln_il6_N_tr$sd[2], 2)),
            " ", as.character(round(d23_ln_tnf_N_tr$sd[1], 2)), as.character(round(d23_ln_tnf_N_tr$sd[2], 2)),
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
            " ", as.character(round(d23_ln_igf_N_tr$sd[1], 2)), as.character(round(d23_ln_igf_N_tr$sd[2], 2)))

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

unadjs6 <- c(" ", " ", paste(d23_ln_il1_unadj_L[1], " (", d23_ln_il1_unadj_L[2], ", ", d23_ln_il1_unadj_L[3], ")", sep=""),
               " ", " ", paste(d23_ln_il6_unadj_L[1], " (", d23_ln_il6_unadj_L[2], ", ", d23_ln_il6_unadj_L[3], ")", sep=""),
               " ", " ", paste(d23_ln_tnf_unadj_L[1], " (", d23_ln_tnf_unadj_L[2], ", ", d23_ln_tnf_unadj_L[3], ")", sep=""),
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
               " ", " ", paste(d23_ln_igf_unadj_L[1], " (", d23_ln_igf_unadj_L[2], ", ", d23_ln_igf_unadj_L[3], ")", sep=""))

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

asadjs6 <- c(" ", " ", paste(d23_ln_il1_adj_sex_age_L[1], " (", d23_ln_il1_adj_sex_age_L[2], ", ", d23_ln_il1_adj_sex_age_L[3], ")", sep=""),
               " ", " ", paste(d23_ln_il6_adj_sex_age_L[1], " (", d23_ln_il6_adj_sex_age_L[2], ", ", d23_ln_il6_adj_sex_age_L[3], ")", sep=""),
               " ", " ", paste(d23_ln_tnf_adj_sex_age_L[1], " (", d23_ln_tnf_adj_sex_age_L[2], ", ", d23_ln_tnf_adj_sex_age_L[3], ")", sep=""),
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
               " ", " ", paste(d23_ln_igf_adj_sex_age_L[1], " (", d23_ln_igf_adj_sex_age_L[2], ", ", d23_ln_igf_adj_sex_age_L[3], ")", sep=""))

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

adjs6 <- c(" ", " ", paste(d23_ln_il1_adj_L[1], " (", d23_ln_il1_adj_L[2], ", ", d23_ln_il1_adj_L[3], ")", sep=""),
             " ", " ", paste(d23_ln_il6_adj_L[1], " (", d23_ln_il6_adj_L[2], ", ", d23_ln_il6_adj_L[3], ")", sep=""),
             " ", " ", paste(d23_ln_tnf_adj_L[1], " (", d23_ln_tnf_adj_L[2], ", ", d23_ln_tnf_adj_L[3], ")", sep=""),
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
             " ", " ", paste(d23_ln_igf_adj_L[1], " (", d23_ln_igf_adj_L[2], ", ", d23_ln_igf_adj_L[3], ")", sep=""))

ipcws6<-c(" ", " ", maketblvalue(d23_il1_adj_ipcw_L$`unlist(d23_il1_adj_ipcw$estimates$ATE)`),
          " ", " ", maketblvalue(d23_il6_adj_ipcw_L$`unlist(d23_il6_adj_ipcw$estimates$ATE)`),
          " ", " ", maketblvalue(d23_tnf_adj_ipcw_L$`unlist(d23_tnf_adj_ipcw$estimates$ATE)`),
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
          " ", " ", maketblvalue(d23_igf_adj_ipcw_L$`unlist(d23_igf_adj_ipcw$estimates$ATE)`))

# Table S6: Effect of intervention on change in individual immune status and growth factor measurements between ages 14 and 28 months
tbls6 <- data.table(
  "Outcome, Arm" = outcomes6,
  "N" = Ns6, 
  "Absolute Mean" = absmeans6,
  "Absolute SD" = abssds6,
  "Mean" = means6, 
  "SD" = sds6,
  "Unadjusted difference: Intervention vs. Control (95% CI)" = unadjs6,
  "Age- and sex- adjusted difference: Intervention vs. Control (95% CI)" = asadjs6, 
  "Fully adjusted difference: Intervention vs. Control (95% CI)" = adjs6,
  "IPCW adjusted difference: Intervention vs. Control (95% CI)" = ipcws6
)

write.csv(tbls6, file=here('tables/supplementary/immune_supptable6.csv'))
print(xtable(tbls6), type="html", file=here("tables/supplementary/immune_supptable6.html"))



#### TABLE S7 ####
makeipcw<-function(var){
  rounded<-round(var, 2)
  paste(rounded[1], " (", rounded[3], ", ", rounded[4], ")", sep="")
}

ipcws7 <- c(" ", " ", makeipcw(d23_ratio_il1_il10_adj_ipcw_L$`unlist(d23_ratio_il1_il10_adj_ipcw$estimates$ATE)`),                                      " ", " ", makeipcw(d23_ratio_il6_il10_adj_ipcw_L$`unlist(d23_ratio_il6_il10_adj_ipcw$estimates$ATE)`),
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

tbls7<-cbind(tbl4, ipcws7)
names(tbls7)[10]<-"IPCW adjusted difference: Intervention vs. Control (95% CI)"

write.csv(tbls7, file=here('tables/supplementary/immune_supptable7.csv'))
print(xtable(tbls7), type="html", file=here("tables/supplementary/immune_supptable7.html"))


