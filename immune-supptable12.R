rm(list=ls())
library("xtable")
source(here::here("0-config.R"))

load(here("results/immune_subgroup.RData"))


bonpval <- function(pval){
  bon = round(pval * 2, 2)
  if (pval >= .5)
    bon = 1
  bon 
}

#Calculations
filtering <- function(row){
  any(!is.na(row))
}

y1<-ages[apply(select(ages, grep("t2_ln", names(ages), ignore.case=T)), 1, filtering),]

y2<-ages[apply(select(ages, grep("t3_ln", names(ages), ignore.case=T)), 1, filtering),]


#calculating overall N by arm
y1Nctrl<-sum(!is.na(y1$tr[y1$tr=="Control"]))
y1Nwsh<-sum(!is.na(y1$tr[y1$tr=="Nutrition + WSH"]))
y2Nctrl<-sum(!is.na(y2$tr[y2$tr=="Control"]))
y2Nwsh<-sum(!is.na(y2$tr[y2$tr=="Nutrition + WSH"]))

#functions for calculating %/mean for all variables in table based on arm
meansdfunc <- function(tbl, variable) {
  ctrlmean<-round(mean(variable[tbl$tr=="Control"], na.rm=TRUE))
  ctrlsd<-round(sd(variable[tbl$tr=="Control"], na.rm=TRUE))
  wshmean<-round(mean(variable[tbl$tr=="Nutrition + WSH"], na.rm=TRUE))
  wshsd<-round(sd(variable[tbl$tr=="Nutrition + WSH"], na.rm=TRUE))
  c(ctrlmean, ctrlsd, wshmean, wshsd)
}

npercfunc <- function(tbl, variable) {
  ctrln<-sum(variable[tbl$tr=="Control"], na.rm=TRUE)
  ctrlperc<-round(mean(variable[tbl$tr=="Control"], na.rm=TRUE)*100)
  wshn<-sum(variable[tbl$tr=="Nutrition + WSH"], na.rm=TRUE)
  wshperc<-round(mean(variable[tbl$tr=="Nutrition + WSH"], na.rm=TRUE)*100)
  c(ctrln, ctrlperc, wshn, wshperc)
}

y1momage<-meansdfunc(y1, y1$momage)
y1momeduy<-meansdfunc(y1, y1$momeduy)
y1dadeduy<-meansdfunc(y1, y1$dadeduy)
y1dadagri<-npercfunc(y1, y1$dadagri)
y1Nhh<-meansdfunc(y1, y1$Nhh)
y1elec<-npercfunc(y1, y1$elec)
y1cement<-npercfunc(y1, y1$cement)

y2momage<-meansdfunc(y2, y2$momage)
y2momeduy<-meansdfunc(y2, y2$momeduy)
y2dadeduy<-meansdfunc(y2, y2$dadeduy)
y2dadagri<-npercfunc(y2, y2$dadagri)
y2Nhh<-meansdfunc(y2, y2$Nhh)
y2elec<-npercfunc(y2, y2$elec)
y2cement<-npercfunc(y2, y2$cement)

y1acresctrlm<-round(mean(y1$landacre[y1$tr=="Control"], na.rm=TRUE), 2)
y1acresctrlsd<-round(sd(y1$landacre[y1$tr=="Control"], na.rm=TRUE), 2)
y1acreswshm<-round(mean(y1$landacre[y1$tr=="Nutrition + WSH"], na.rm=TRUE), 2)
y1acreswshsd<-round(mean(y1$landacre[y1$tr=="Nutrition + WSH"], na.rm=TRUE), 2)
y1acres<-c(y1acresctrlm, y1acresctrlsd, y1acreswshm, y1acreswshsd)

y2acresctrlm<-round(mean(y2$landacre[y2$tr=="Control"], na.rm=TRUE), 2)
y2acresctrlsd<-round(sd(y2$landacre[y2$tr=="Control"], na.rm=TRUE), 2)
y2acreswshm<-round(mean(y2$landacre[y2$tr=="Nutrition + WSH"], na.rm=TRUE), 2)
y2acreswshsd<-round(mean(y2$landacre[y2$tr=="Nutrition + WSH"], na.rm=TRUE), 2)
y2acres<-c(y2acresctrlm, y2acresctrlsd, y2acreswshm, y2acreswshsd)

y1tubewell<-npercfunc(y1, y1$tubewell)
y1storewater<-npercfunc(y1, y1$storewat)
y1treatwater<-npercfunc(y1, y1$treatwat)
y1waterdis<-meansdfunc(y1, y1$watmin)
y1odmen<-npercfunc(y1, y1$odmen)
y1odwomen<-npercfunc(y1, y1$odwom)
y1odchild815<-npercfunc(y1, y1$odch815)
y1odchild38<-npercfunc(y1, y1$odch38)
y1odchild03<-npercfunc(y1, y1$odchu3)
y1latowned<-npercfunc(y1, y1$latown)
y1latslab<-npercfunc(y1, y1$latslab)
y1latseal<-npercfunc(y1, y1$latseal)

y1latfctrln<-sum(y1$latfeces[y1$tr=="Control"]==0, na.rm=T)
y1latfctrlperc<-round(y1latfctrln/sum(!is.na(y1$latfeces[y1$tr=="Control"]), na.rm=T)*100)
y1latfwshn<-sum(y1$latfeces[y1$tr=="Nutrition + WSH"]==0, na.rm=T)
y1latfwshperc<-round(y1latfwshn/sum(!is.na(y1$latfeces[y1$tr=="Nutrition + WSH"]), na.rm=T)*100)
y1latfeces<-c(y1latfctrln, y1latfctrlperc, y1latfwshn, y1latfwshperc)

y1potty<-npercfunc(y1, y1$potty)
y1feceshouse<-npercfunc(y1, y1$humfeces)
y1feceschildarea<-npercfunc(y1, y1$humfecesch)
y1handlatwater<-npercfunc(y1, y1$hwlatwat)
y1handlatsoap<-npercfunc(y1, y1$hwlatsoap)
y1handkitwater<-npercfunc(y1, y1$hwkitwat)
y1handkitsoap<-npercfunc(y1, y1$hwkitsoap)

y2tubewell<-npercfunc(y2, y2$tubewell)
y2storewater<-npercfunc(y2, y2$storewat)
y2treatwater<-npercfunc(y2, y2$treatwat)
y2waterdis<-meansdfunc(y2, y2$watmin)
y2odmen<-npercfunc(y2, y2$odmen)
y2odwomen<-npercfunc(y2, y2$odwom)
y2odchild815<-npercfunc(y2, y2$odch815)
y2odchild38<-npercfunc(y2, y2$odch38)
y2odchild03<-npercfunc(y2, y2$odchu3)
y2latowned<-npercfunc(y2, y2$latown)
y2latslab<-npercfunc(y2, y2$latslab)
y2latseal<-npercfunc(y2, y2$latseal)

y2latfctrln<-sum(y2$latfeces[y2$tr=="Control"]==0, na.rm=T)
y2latfctrlperc<-round(y2latfctrln/sum(!is.na(y2$latfeces[y2$tr=="Control"]), na.rm=T)*100)
y2latfwshn<-sum(y2$latfeces[y2$tr=="Nutrition + WSH"]==0, na.rm=T)
y2latfwshperc<-round(y2latfwshn/sum(!is.na(y2$latfeces[y2$tr=="Nutrition + WSH"]), na.rm=T)*100)
y2latfeces<-c(y2latfctrln, y2latfctrlperc, y2latfwshn, y2latfwshperc)

y2potty<-npercfunc(y2, y2$potty)
y2feceshouse<-npercfunc(y2, y2$humfeces)
y2feceschildarea<-npercfunc(y2, y2$humfecesch)
y2handlatwater<-npercfunc(y2, y2$hwlatwat)
y2handlatsoap<-npercfunc(y2, y2$hwlatsoap)
y2handkitwater<-npercfunc(y2, y2$hwkitwat)
y2handkitsoap<-npercfunc(y2, y2$hwkitsoap)

y1fsctrln<-sum(y1$hfiacat[y1$tr=="Control"]=="Food Secure", na.rm=T)
y1fsctrlperc<-round(y1fsctrln/sum(!is.na(y1$hfiacat[y1$tr=="Control"]))*100)
y1fswshn<-sum(y1$hfiacat[y1$tr=="Nutrition + WSH"]=="Food Secure", na.rm=T)
y1fswshperc<-round(y1fswshn/sum(!is.na(y1$hfiacat[y1$tr=="Nutrition + WSH"]))*100)
y1foodsecure<-c(y1fsctrln, y1fsctrlperc, y1fswshn, y1fswshperc)

y2fsctrln<-sum(y2$hfiacat[y2$tr=="Control"]=="Food Secure", na.rm=T)
y2fsctrlperc<-round(y2fsctrln/sum(!is.na(y2$hfiacat[y2$tr=="Control"]))*100)
y2fswshn<-sum(y2$hfiacat[y2$tr=="Nutrition + WSH"]=="Food Secure", na.rm=T)
y2fswshperc<-round(y2fswshn/sum(!is.na(y2$hfiacat[y2$tr=="Nutrition + WSH"]))*100)
y2foodsecure<-c(y2fsctrln, y2fsctrlperc, y2fswshn, y2fswshperc)


#make vectors to put in table
#function combines n and percent or mean and sd for vectors created from npercfunc or meansdfunc
#num is 1 if ctrl group, 3 if wsh
charobject<-function(variable, num) {
  paste(variable[num], " (", variable[num+1], ")", sep="")
}

charobjectperc<-function(variable, num) {
  paste(variable[num], " (", variable[num+1], "%)", sep="")
}

ctrly1<-c(paste("Control (N=", y1Nctrl, ")", sep=""), " ", charobject(y1momage, 1),charobject(y1momeduy, 1), " ", charobject(y1dadeduy, 1), charobjectperc(y1dadagri, 1),
          " ", charobject(y1Nhh, 1), charobjectperc(y1elec, 1), charobjectperc(y1cement, 1), charobject(y1acres, 1),
          " ", charobjectperc(y1tubewell, 1), charobjectperc(y1storewater, 1), charobjectperc(y1treatwater, 1), charobject(y1waterdis, 1), 
          " ", " ", charobjectperc(y1odmen, 1), charobjectperc(y1odwomen, 1), charobjectperc(y1odchild815, 1), charobjectperc(y1odchild38, 1), charobjectperc(y1odchild03, 1), 
          " ", charobjectperc(y1latowned, 1), charobjectperc(y1latslab, 1), charobjectperc(y1latseal, 1), charobjectperc(y1latfeces, 1),
          charobjectperc(y1potty, 1), 
          " ", charobjectperc(y1feceshouse, 1), charobjectperc(y1feceschildarea, 1), 
          " ", " ", charobjectperc(y1handlatwater, 1), charobjectperc(y1handlatsoap, 1), 
          " ", charobjectperc(y1handkitwater, 1), charobjectperc(y1handkitsoap, 1), 
          " ", charobjectperc(y1foodsecure, 1))

wshy1<-c(paste("N+WSH Intervention (N=", y1Nwsh, ")", sep=""), " ", charobject(y1momage, 3),charobject(y1momeduy, 3), " ", charobject(y1dadeduy, 3), charobjectperc(y1dadagri, 3),
         " ", charobject(y1Nhh, 3), charobjectperc(y1elec, 3), charobjectperc(y1cement, 3), charobject(y1acres, 3),
         " ", charobjectperc(y1tubewell, 3), charobjectperc(y1storewater, 3), charobjectperc(y1treatwater, 3), charobject(y1waterdis, 3), 
         " ", " ", charobjectperc(y1odmen, 3), charobjectperc(y1odwomen, 3), charobjectperc(y1odchild815, 3), charobjectperc(y1odchild38, 3), charobjectperc(y1odchild03, 3), 
         " ", charobjectperc(y1latowned, 3), charobjectperc(y1latslab, 3), charobjectperc(y1latseal, 3), charobjectperc(y1latfeces, 3),
         charobjectperc(y1potty, 3), 
         " ", charobjectperc(y1feceshouse, 3), charobjectperc(y1feceschildarea, 3), 
         " ", " ", charobjectperc(y1handlatwater, 3), charobjectperc(y1handlatsoap, 3), 
         " ", charobjectperc(y1handkitwater, 3), charobjectperc(y1handkitsoap, 3), 
         " ", charobjectperc(y1foodsecure, 3))

ctrly2<-c(paste("Control (N=", y2Nctrl, ")", sep=""), " ", charobject(y2momage, 1),charobject(y2momeduy, 1), " ", charobject(y2dadeduy, 1), charobjectperc(y2dadagri, 1),
          " ", charobject(y2Nhh, 1), charobjectperc(y2elec, 1), charobjectperc(y2cement, 1), charobject(y2acres, 1),
          " ", charobjectperc(y2tubewell, 1), charobjectperc(y2storewater, 1), charobjectperc(y2treatwater, 1), charobject(y2waterdis, 1), 
          " ", " ", charobjectperc(y2odmen, 1), charobjectperc(y2odwomen, 1), charobjectperc(y2odchild815, 1), charobjectperc(y2odchild38, 1), charobjectperc(y2odchild03, 1), 
          " ", charobjectperc(y2latowned, 1), charobjectperc(y2latslab, 1), charobjectperc(y2latseal, 1), charobjectperc(y2latfeces, 1),
          charobjectperc(y2potty, 1), 
          " ", charobjectperc(y2feceshouse, 1), charobjectperc(y2feceschildarea, 1), 
          " ", " ", charobjectperc(y2handlatwater, 1), charobjectperc(y2handlatsoap, 1), 
          " ", charobjectperc(y2handkitwater, 1), charobjectperc(y2handkitsoap, 1), 
          " ", charobjectperc(y2foodsecure, 1))

wshy2<-c(paste("N+WSH Intervention (N=", y2Nwsh, ")", sep=""), " ", charobject(y2momage, 3),charobject(y2momeduy, 3), " ", charobject(y2dadeduy, 3), charobjectperc(y2dadagri, 3),
         " ", charobject(y2Nhh, 3), charobjectperc(y2elec, 3), charobjectperc(y2cement, 3), charobject(y2acres, 3),
         " ", charobjectperc(y2tubewell, 3), charobjectperc(y2storewater, 3), charobjectperc(y2treatwater, 3), charobject(y2waterdis, 3), 
         " ", " ", charobjectperc(y2odmen, 3), charobjectperc(y2odwomen, 3), charobjectperc(y2odchild815, 3), charobjectperc(y2odchild38, 3), charobjectperc(y2odchild03, 3), 
         " ", charobjectperc(y2latowned, 3), charobjectperc(y2latslab, 3), charobjectperc(y2latseal, 3), charobjectperc(y2latfeces, 3),
         charobjectperc(y2potty, 3), 
         " ", charobjectperc(y2feceshouse, 3), charobjectperc(y2feceschildarea, 3), 
         " ", " ", charobjectperc(y2handlatwater, 3), charobjectperc(y2handlatsoap, 3), 
         " ", charobjectperc(y2handkitwater, 3), charobjectperc(y2handkitsoap, 3), 
         " ", charobjectperc(y2foodsecure, 3))

####Table 12 ####

#Headings

outcometbl12 <- c(paste("Ln IL-1", "β", "/IL-10", sep=""), "Control", "Nutrition + WSH", 
                 "Ln IL-6/IL-10", "Control", "Nutrition + WSH", 
                 paste("Ln TNF-", "α", "/IL-10", sep=""), "Control", "Nutrition + WSH",
                 "Ln IL-12/IL-10", "Control", "Nutrition + WSH", 
                 paste("Ln IFN-", "γ", "/IL-10", sep=""), "Control", "Nutrition + WSH", 
                 "Ln IL-4/IL-10", "Control", "Nutrition + WSH", 
                 "Ln IL-5/IL-10", "Control", "Nutrition + WSH", 
                 "Ln IL-13/IL-10", "Control", "Nutrition + WSH", 
                 "Ln IL-17A/IL-10", "Control", "Nutrition + WSH", 
                 "Ln IL-21/IL-10", "Control", "Nutrition + WSH", 
                 "Ln IL-2/IL-10", "Control", "Nutrition + WSH", 
                 "Ln GM-CSF/IL-10", "Control", "Nutrition + WSH",
                 "Ln IL-12/IL-4", "Control", "Nutrition + WSH", 
                 paste("Ln IFN-", "γ", "/IL-4", sep=""), "Control", "Nutrition + WSH",
                 "Ln IL-12/IL-5", "Control", "Nutrition + WSH", 
                 paste("Ln IFN-", "γ", "/IL-5", sep=""), "Control", "Nutrition + WSH",
                 "Ln IL-12/IL-13", "Control", "Nutrition + WSH", 
                 paste("Ln IFN-", "γ", "/IL-13", sep=""), "Control", "Nutrition + WSH",
                 "Ln IL-12/IL-17A", "Control", "Nutrition + WSH", 
                 paste("Ln IFN-", "γ", "/IL-17A", sep=""), "Control", "Nutrition + WSH",
                 "Ln IL-12/IL-21", "Control", "Nutrition + WSH", 
                 paste("Ln IFN-", "γ", "/IL-21", sep=""), "Control", "Nutrition + WSH",
                 "Ln Pro-inflammatory cytokines/IL-10", "Control", "Nutrition + WSH",
                 "Ln Th1/IL-10", "Control", "Nutrition + WSH", 
                 "Ln Th2/IL-10", "Control", "Nutrition + WSH", 
                 "Ln Th17/IL-10", "Control", "Nutrition + WSH", 
                 "Ln Th1/Th2", "Control", "Nutrition + WSH", 
                 "Ln Th1/Th17", "Control", "Nutrition + WSH")

#Calculate N

#Female
#X$subgroup==1

Ntbl12f <- c("N", "" , as.character(t2_ratio_il1_il10_subgroup_L$t2_ratio_il1_il10_subgroup_L[1]), as.character(t2_ratio_il1_il10_subgroup_L$t2_ratio_il1_il10_subgroup_L[2]), 
           " ", as.character(t2_ratio_il6_il10_subgroup_L$t2_ratio_il6_il10_subgroup_L[1]), as.character(t2_ratio_il6_il10_subgroup_L$t2_ratio_il6_il10_subgroup_L[2]),  
           " ", as.character(t2_ratio_tnf_il10_subgroup_L$t2_ratio_tnf_il10_subgroup_L[1]), as.character(t2_ratio_tnf_il10_subgroup_L$t2_ratio_tnf_il10_subgroup_L[2]), 
           " ", as.character(t2_ratio_il12_il10_subgroup_L$t2_ratio_il12_il10_subgroup_L[1]), as.character(t2_ratio_il12_il10_subgroup_L$t2_ratio_il12_il10_subgroup_L[2]), 
           " ", as.character(t2_ratio_ifn_il10_subgroup_L$t2_ratio_ifn_il10_subgroup_L[1]), as.character(t2_ratio_ifn_il10_subgroup_L$t2_ratio_ifn_il10_subgroup_L[2]), 
           " ", as.character(t2_ratio_il4_il10_subgroup_L$t2_ratio_il4_il10_subgroup_L[1]), as.character(t2_ratio_il4_il10_subgroup_L$t2_ratio_il4_il10_subgroup_L[2]), 
           " ", as.character(t2_ratio_il5_il10_subgroup_L$t2_ratio_il5_il10_subgroup_L[1]), as.character(t2_ratio_il5_il10_subgroup_L$t2_ratio_il5_il10_subgroup_L[2]),  
           " ", as.character(t2_ratio_il13_il10_subgroup_L$t2_ratio_il13_il10_subgroup_L[1]), as.character(t2_ratio_il13_il10_subgroup_L$t2_ratio_il13_il10_subgroup_L[2]),    
           " ", as.character(t2_ratio_il17_il10_subgroup_L$t2_ratio_il17_il10_subgroup_L[1]), as.character(t2_ratio_il17_il10_subgroup_L$t2_ratio_il17_il10_subgroup_L[2]),  
           " ", as.character(t2_ratio_il21_il10_subgroup_L$t2_ratio_il21_il10_subgroup_L[1]), as.character(t2_ratio_il21_il10_subgroup_L$t2_ratio_il21_il10_subgroup_L[2]),  
           " ", as.character(t2_ratio_il2_il10_subgroup_L$t2_ratio_il2_il10_subgroup_L[1]), as.character(t2_ratio_il2_il10_subgroup_L$t2_ratio_il2_il10_subgroup_L[2]),  
           " ", as.character(t2_ratio_gmc_il10_subgroup_L$t2_ratio_gmc_il10_subgroup_L[1]), as.character(t2_ratio_gmc_il10_subgroup_L$t2_ratio_gmc_il10_subgroup_L[2]),  
           " ", as.character(t2_ratio_il12_il4_subgroup_L$t2_ratio_il12_il4_subgroup_L[1]), as.character(t2_ratio_il12_il4_subgroup_L$t2_ratio_il12_il4_subgroup_L[2]),  
           " ", as.character(t2_ratio_ifn_il4_subgroup_L$t2_ratio_ifn_il4_subgroup_L[1]), as.character(t2_ratio_ifn_il4_subgroup_L$t2_ratio_ifn_il4_subgroup_L[2]),  
           " ", as.character(t2_ratio_il12_il5_subgroup_L$t2_ratio_il12_il5_subgroup_L[1]), as.character(t2_ratio_il12_il5_subgroup_L$t2_ratio_il12_il5_subgroup_L[2]),  
           " ", as.character(t2_ratio_ifn_il5_subgroup_L$t2_ratio_ifn_il5_subgroup_L[1]), as.character(t2_ratio_ifn_il5_subgroup_L$t2_ratio_ifn_il5_subgroup_L[2]), 
           " ", as.character(t2_ratio_il12_il13_subgroup_L$t2_ratio_il12_il13_subgroup_L[1]), as.character(t2_ratio_il12_il13_subgroup_L$t2_ratio_il12_il13_subgroup_L[2]), 
           " ", as.character(t2_ratio_ifn_il13_subgroup_L$t2_ratio_ifn_il13_subgroup_L[1]), as.character(t2_ratio_ifn_il13_subgroup_L$t2_ratio_ifn_il13_subgroup_L[2]), 
           " ", as.character(t2_ratio_il12_il17_subgroup_L$t2_ratio_il12_il17_subgroup_L[1]), as.character(t2_ratio_il12_il17_subgroup_L$t2_ratio_il12_il17_subgroup_L[2]),  
           " ", as.character(t2_ratio_ifn_il17_subgroup_L$t2_ratio_ifn_il17_subgroup_L[1]), as.character(t2_ratio_ifn_il17_subgroup_L$t2_ratio_ifn_il17_subgroup_L[2]),  
           " ", as.character(t2_ratio_il12_il21_subgroup_L$t2_ratio_il12_il21_subgroup_L[1]), as.character(t2_ratio_il12_il21_subgroup_L$t2_ratio_il12_il21_subgroup_L[2]), 
           " ", as.character(t2_ratio_ifn_il21_subgroup_L$t2_ratio_ifn_il21_subgroup_L[1]), as.character(t2_ratio_ifn_il21_subgroup_L$t2_ratio_ifn_il21_subgroup_L[2]),
           " ", as.character(t2_ratio_pro_il10_subgroup_L$t2_ratio_pro_il10_subgroup_L[1]), as.character(t2_ratio_pro_il10_subgroup_L$t2_ratio_pro_il10_subgroup_L[2]),  
           " ", as.character(t2_ratio_th1_il10_subgroup_L$t2_ratio_th1_il10_subgroup_L[1]), as.character(t2_ratio_th1_il10_subgroup_L$t2_ratio_th1_il10_subgroup_L[2]),  
           " ", as.character(t2_ratio_th2_il10_subgroup_L$t2_ratio_th2_il10_subgroup_L[1]), as.character(t2_ratio_th2_il10_subgroup_L$t2_ratio_th2_il10_subgroup_L[2]),    
           " ", as.character(t2_ratio_th17_il10_subgroup_L$t2_ratio_th17_il10_subgroup_L[1]), as.character(t2_ratio_th17_il10_subgroup_L$t2_ratio_th17_il10_subgroup_L[2]), 
           " ", as.character(t2_ratio_th1_th2_subgroup_L$t2_ratio_th1_th2_subgroup_L[1]), as.character(t2_ratio_th1_th2_subgroup_L$t2_ratio_th1_th2_subgroup_L[2]),  
           " ", as.character(t2_ratio_th1_th17_subgroup_L$t2_ratio_th1_th17_subgroup_L[1]), as.character(t2_ratio_th1_th17_subgroup_L$t2_ratio_th1_th17_subgroup_L[2]))

#Male
#X$subgroup==2
Ntbl12m <- c("N", "" , as.character(t2_ratio_il1_il10_subgroup_L$t2_ratio_il1_il10_subgroup_L[1]), as.character(t2_ratio_il1_il10_subgroup_L$t2_ratio_il1_il10_subgroup_L[2]), 
            " ", as.character(t2_ratio_il6_il10_subgroup_L$t2_ratio_il6_il10_subgroup_L[1]), as.character(t2_ratio_il6_il10_subgroup_L$t2_ratio_il6_il10_subgroup_L[2]),  
            " ", as.character(t2_ratio_tnf_il10_subgroup_L$t2_ratio_tnf_il10_subgroup_L[1]), as.character(t2_ratio_tnf_il10_subgroup_L$t2_ratio_tnf_il10_subgroup_L[2]), 
            " ", as.character(t2_ratio_il12_il10_subgroup_L$t2_ratio_il12_il10_subgroup_L[1]), as.character(t2_ratio_il12_il10_subgroup_L$t2_ratio_il12_il10_subgroup_L[2]), 
            " ", as.character(t2_ratio_ifn_il10_subgroup_L$t2_ratio_ifn_il10_subgroup_L[1]), as.character(t2_ratio_ifn_il10_subgroup_L$t2_ratio_ifn_il10_subgroup_L[2]), 
            " ", as.character(t2_ratio_il4_il10_subgroup_L$t2_ratio_il4_il10_subgroup_L[1]), as.character(t2_ratio_il4_il10_subgroup_L$t2_ratio_il4_il10_subgroup_L[2]), 
            " ", as.character(t2_ratio_il5_il10_subgroup_L$t2_ratio_il5_il10_subgroup_L[1]), as.character(t2_ratio_il5_il10_subgroup_L$t2_ratio_il5_il10_subgroup_L[2]),  
            " ", as.character(t2_ratio_il13_il10_subgroup_L$t2_ratio_il13_il10_subgroup_L[1]), as.character(t2_ratio_il13_il10_subgroup_L$t2_ratio_il13_il10_subgroup_L[2]),    
            " ", as.character(t2_ratio_il17_il10_subgroup_L$t2_ratio_il17_il10_subgroup_L[1]), as.character(t2_ratio_il17_il10_subgroup_L$t2_ratio_il17_il10_subgroup_L[2]),  
            " ", as.character(t2_ratio_il21_il10_subgroup_L$t2_ratio_il21_il10_subgroup_L[1]), as.character(t2_ratio_il21_il10_subgroup_L$t2_ratio_il21_il10_subgroup_L[2]),  
            " ", as.character(t2_ratio_il2_il10_subgroup_L$t2_ratio_il2_il10_subgroup_L[1]), as.character(t2_ratio_il2_il10_subgroup_L$t2_ratio_il2_il10_subgroup_L[2]),  
            " ", as.character(t2_ratio_gmc_il10_subgroup_L$t2_ratio_gmc_il10_subgroup_L[1]), as.character(t2_ratio_gmc_il10_subgroup_L$t2_ratio_gmc_il10_subgroup_L[2]),  
            " ", as.character(t2_ratio_il12_il4_subgroup_L$t2_ratio_il12_il4_subgroup_L[1]), as.character(t2_ratio_il12_il4_subgroup_L$t2_ratio_il12_il4_subgroup_L[2]),  
            " ", as.character(t2_ratio_ifn_il4_subgroup_L$t2_ratio_ifn_il4_subgroup_L[1]), as.character(t2_ratio_ifn_il4_subgroup_L$t2_ratio_ifn_il4_subgroup_L[2]),  
            " ", as.character(t2_ratio_il12_il5_subgroup_L$t2_ratio_il12_il5_subgroup_L[1]), as.character(t2_ratio_il12_il5_subgroup_L$t2_ratio_il12_il5_subgroup_L[2]),  
            " ", as.character(t2_ratio_ifn_il5_subgroup_L$t2_ratio_ifn_il5_subgroup_L[1]), as.character(t2_ratio_ifn_il5_subgroup_L$t2_ratio_ifn_il5_subgroup_L[2]), 
            " ", as.character(t2_ratio_il12_il13_subgroup_L$t2_ratio_il12_il13_subgroup_L[1]), as.character(t2_ratio_il12_il13_subgroup_L$t2_ratio_il12_il13_subgroup_L[2]), 
            " ", as.character(t2_ratio_ifn_il13_subgroup_L$t2_ratio_ifn_il13_subgroup_L[1]), as.character(t2_ratio_ifn_il13_subgroup_L$t2_ratio_ifn_il13_subgroup_L[2]), 
            " ", as.character(t2_ratio_il12_il17_subgroup_L$t2_ratio_il12_il17_subgroup_L[1]), as.character(t2_ratio_il12_il17_subgroup_L$t2_ratio_il12_il17_subgroup_L[2]),  
            " ", as.character(t2_ratio_ifn_il17_subgroup_L$t2_ratio_ifn_il17_subgroup_L[1]), as.character(t2_ratio_ifn_il17_subgroup_L$t2_ratio_ifn_il17_subgroup_L[2]),  
            " ", as.character(t2_ratio_il12_il21_subgroup_L$t2_ratio_il12_il21_subgroup_L[1]), as.character(t2_ratio_il12_il21_subgroup_L$t2_ratio_il12_il21_subgroup_L[2]), 
            " ", as.character(t2_ratio_ifn_il21_subgroup_L$t2_ratio_ifn_il21_subgroup_L[1]), as.character(t2_ratio_ifn_il21_subgroup_L$t2_ratio_ifn_il21_subgroup_L[2]),
            " ", as.character(t2_ratio_pro_il10_subgroup_L$t2_ratio_pro_il10_subgroup_L[1]), as.character(t2_ratio_pro_il10_subgroup_L$t2_ratio_pro_il10_subgroup_L[2]),  
            " ", as.character(t2_ratio_th1_il10_subgroup_L$t2_ratio_th1_il10_subgroup_L[1]), as.character(t2_ratio_th1_il10_subgroup_L$t2_ratio_th1_il10_subgroup_L[2]),  
            " ", as.character(t2_ratio_th2_il10_subgroup_L$t2_ratio_th2_il10_subgroup_L[1]), as.character(t2_ratio_th2_il10_subgroup_L$t2_ratio_th2_il10_subgroup_L[2]),    
            " ", as.character(t2_ratio_th17_il10_subgroup_L$t2_ratio_th17_il10_subgroup_L[1]), as.character(t2_ratio_th17_il10_subgroup_L$t2_ratio_th17_il10_subgroup_L[2]), 
            " ", as.character(t2_ratio_th1_th2_subgroup_L$t2_ratio_th1_th2_subgroup_L[1]), as.character(t2_ratio_th1_th2_subgroup_L$t2_ratio_th1_th2_subgroup_L[2]),  
            " ", as.character(t2_ratio_th1_th17_subgroup_L$t2_ratio_th1_th17_subgroup_L[1]), as.character(t2_ratio_th1_th17_subgroup_L$t2_ratio_th1_th17_subgroup_L[2]))

#Risk difference
#Subset to male or female with below code
#female
summary(t2_ratio_th1_th17_subgroup_L$subgroup[1])

rdtbl12f <- c("Risk difference" , " ", as.character(round(t2_ratio_il1_il10_subgroup_L$RD[1], 2)), as.character(round(t2_ratio_il1_il10_subgroup_L$RD[2], 2)), 
              " ", as.character(round(t2_ratio_il6_il10_subgroup_L$RD[1], 2)), as.character(round(t2_ratio_il6_il10_subgroup_L$RD[2], 2)),  
              " ", as.character(round(t2_ratio_tnf_il10_subgroup_L$RD[1], 2)), as.character(round(t2_ratio_tnf_il10_subgroup_L$RD[2], 2)), 
              " ", as.character(round(t2_ratio_il12_il10_subgroup_L$RD[1], 2)), as.character(round(t2_ratio_il12_il10_subgroup_L$RD[2], 2)), 
              " ", as.character(round(t2_ratio_ifn_il10_subgroup_L$RD[1], 2)), as.character(round(t2_ratio_ifn_il10_subgroup_L$RD[2], 2)), 
              " ", as.character(round(t2_ratio_il4_il10_subgroup_L$RD[1], 2)), as.character(round(t2_ratio_il4_il10_subgroup_L$RD[2], 2)), 
              " ", as.character(round(t2_ratio_il5_il10_subgroup_L$RD[1], 2)), as.character(round(t2_ratio_il5_il10_subgroup_L$RD[2], 2)),  
              " ", as.character(round(t2_ratio_il13_il10_subgroup_L$RD[1], 2)), as.character(round(t2_ratio_il13_il10_subgroup_L$RD[2], 2)),    
              " ", as.character(round(t2_ratio_il17_il10_subgroup_L$RD[1], 2)), as.character(round(t2_ratio_il17_il10_subgroup_L$RD[2], 2)),  
              " ", as.character(round(t2_ratio_il21_il10_subgroup_L$RD[1], 2)), as.character(round(t2_ratio_il21_il10_subgroup_L$RD[2], 2)),  
              " ", as.character(round(t2_ratio_il2_il10_subgroup_L$RD[1], 2)), as.character(round(t2_ratio_il2_il10_subgroup_L$RD[2], 2)),  
              " ", as.character(round(t2_ratio_gmc_il10_subgroup_L$RD[1], 2)), as.character(round(t2_ratio_gmc_il10_subgroup_L$RD[2], 2)),  
              " ", as.character(round(t2_ratio_il12_il4_subgroup_L$RD[1], 2)), as.character(round(t2_ratio_il12_il4_subgroup_L$RD[2], 2)),  
              " ", as.character(round(t2_ratio_ifn_il4_subgroup_L$RD[1], 2)), as.character(round(t2_ratio_ifn_il4_subgroup_L$RD[2], 2)),  
              " ", as.character(round(t2_ratio_il12_il5_subgroup_L$RD[1], 2)), as.character(round(t2_ratio_il12_il5_subgroup_L$RD[2], 2)),  
              " ", as.character(round(t2_ratio_ifn_il5_subgroup_L$RD[1], 2)), as.character(round(t2_ratio_ifn_il5_subgroup_L$RD[2], 2)), 
              " ", as.character(round(t2_ratio_il12_il13_subgroup_L$RD[1], 2)), as.character(round(t2_ratio_il12_il13_subgroup_L$RD[2], 2)), 
              " ", as.character(round(t2_ratio_ifn_il13_subgroup_L$RD[1], 2)), as.character(round(t2_ratio_ifn_il13_subgroup_L$RD[2], 2)), 
              " ", as.character(round(t2_ratio_il12_il17_subgroup_L$RD[1], 2)), as.character(round(t2_ratio_il12_il17_subgroup_L$RD[2], 2)),  
              " ", as.character(round(t2_ratio_ifn_il17_subgroup_L$RD[1], 2)), as.character(round(t2_ratio_ifn_il17_subgroup_L$RD[2], 2)),  
              " ", as.character(round(t2_ratio_il12_il21_subgroup_L$RD[1], 2)), as.character(round(t2_ratio_il12_il21_subgroup_L$RD[2], 2)), 
              " ", as.character(round(t2_ratio_ifn_il21_subgroup_L$RD[1], 2)), as.character(round(t2_ratio_ifn_il21_subgroup_L$RD[2], 2)),
              " ", as.character(round(t2_ratio_pro_il10_subgroup_L$RD[1], 2)), as.character(round(t2_ratio_pro_il10_subgroup_L$RD[2], 2)),  
              " ", as.character(round(t2_ratio_th1_il10_subgroup_L$RD[1], 2)), as.character(round(t2_ratio_th1_il10_subgroup_L$RD[2], 2)),  
              " ", as.character(round(t2_ratio_th2_il10_subgroup_L$RD[1], 2)), as.character(round(t2_ratio_th2_il10_subgroup_L$RD[2], 2)),    
              " ", as.character(round(t2_ratio_th17_il10_subgroup_L$RD[1], 2)), as.character(round(t2_ratio_th17_il10_subgroup_L$RD[2], 2)), 
              " ", as.character(round(t2_ratio_th1_th2_subgroup_L$RD[1], 2)), as.character(round(t2_ratio_th1_th2_subgroup_L$RD[2], 2)),  
              " ", as.character(round(t2_ratio_th1_th17_subgroup_L$RD[1], 2)),  as.character(round(t2_ratio_th1_th17_subgroup_L$RD[2], 2)))

#male
summary(t2_ratio_th1_th17_subgroup_L$subgroup[2])

rdtbl12m <- c("Risk difference" , " ", as.character(round(t2_ratio_il1_il10_subgroup_L$RD[1], 2)), as.character(round(t2_ratio_il1_il10_subgroup_L$RD[2], 2)), 
              " ", as.character(round(t2_ratio_il6_il10_subgroup_L$RD[1], 2)), as.character(round(t2_ratio_il6_il10_subgroup_L$RD[2], 2)),  
              " ", as.character(round(t2_ratio_tnf_il10_subgroup_L$RD[1], 2)), as.character(round(t2_ratio_tnf_il10_subgroup_L$RD[2], 2)), 
              " ", as.character(round(t2_ratio_il12_il10_subgroup_L$RD[1], 2)), as.character(round(t2_ratio_il12_il10_subgroup_L$RD[2], 2)), 
              " ", as.character(round(t2_ratio_ifn_il10_subgroup_L$RD[1], 2)), as.character(round(t2_ratio_ifn_il10_subgroup_L$RD[2], 2)), 
              " ", as.character(round(t2_ratio_il4_il10_subgroup_L$RD[1], 2)), as.character(round(t2_ratio_il4_il10_subgroup_L$RD[2], 2)), 
              " ", as.character(round(t2_ratio_il5_il10_subgroup_L$RD[1], 2)), as.character(round(t2_ratio_il5_il10_subgroup_L$RD[2], 2)),  
              " ", as.character(round(t2_ratio_il13_il10_subgroup_L$RD[1], 2)), as.character(round(t2_ratio_il13_il10_subgroup_L$RD[2], 2)),    
              " ", as.character(round(t2_ratio_il17_il10_subgroup_L$RD[1], 2)), as.character(round(t2_ratio_il17_il10_subgroup_L$RD[2], 2)),  
              " ", as.character(round(t2_ratio_il21_il10_subgroup_L$RD[1], 2)), as.character(round(t2_ratio_il21_il10_subgroup_L$RD[2], 2)),  
              " ", as.character(round(t2_ratio_il2_il10_subgroup_L$RD[1], 2)), as.character(round(t2_ratio_il2_il10_subgroup_L$RD[2], 2)),  
              " ", as.character(round(t2_ratio_gmc_il10_subgroup_L$RD[1], 2)), as.character(round(t2_ratio_gmc_il10_subgroup_L$RD[2], 2)),  
              " ", as.character(round(t2_ratio_il12_il4_subgroup_L$RD[1], 2)), as.character(round(t2_ratio_il12_il4_subgroup_L$RD[2], 2)),  
              " ", as.character(round(t2_ratio_ifn_il4_subgroup_L$RD[1], 2)), as.character(round(t2_ratio_ifn_il4_subgroup_L$RD[2], 2)),  
              " ", as.character(round(t2_ratio_il12_il5_subgroup_L$RD[1], 2)), as.character(round(t2_ratio_il12_il5_subgroup_L$RD[2], 2)),  
              " ", as.character(round(t2_ratio_ifn_il5_subgroup_L$RD[1], 2)), as.character(round(t2_ratio_ifn_il5_subgroup_L$RD[2], 2)), 
              " ", as.character(round(t2_ratio_il12_il13_subgroup_L$RD[1], 2)), as.character(round(t2_ratio_il12_il13_subgroup_L$RD[2], 2)), 
              " ", as.character(round(t2_ratio_ifn_il13_subgroup_L$RD[1], 2)), as.character(round(t2_ratio_ifn_il13_subgroup_L$RD[2], 2)), 
              " ", as.character(round(t2_ratio_il12_il17_subgroup_L$RD[1], 2)), as.character(round(t2_ratio_il12_il17_subgroup_L$RD[2], 2)),  
              " ", as.character(round(t2_ratio_ifn_il17_subgroup_L$RD[1], 2)), as.character(round(t2_ratio_ifn_il17_subgroup_L$RD[2], 2)),  
              " ", as.character(round(t2_ratio_il12_il21_subgroup_L$RD[1], 2)), as.character(round(t2_ratio_il12_il21_subgroup_L$RD[2], 2)), 
              " ", as.character(round(t2_ratio_ifn_il21_subgroup_L$RD[1], 2)), as.character(round(t2_ratio_ifn_il21_subgroup_L$RD[2], 2)),
              " ", as.character(round(t2_ratio_pro_il10_subgroup_L$RD[1], 2)), as.character(round(t2_ratio_pro_il10_subgroup_L$RD[2], 2)),  
              " ", as.character(round(t2_ratio_th1_il10_subgroup_L$RD[1], 2)), as.character(round(t2_ratio_th1_il10_subgroup_L$RD[2], 2)),  
              " ", as.character(round(t2_ratio_th2_il10_subgroup_L$RD[1], 2)), as.character(round(t2_ratio_th2_il10_subgroup_L$RD[2], 2)),    
              " ", as.character(round(t2_ratio_th17_il10_subgroup_L$RD[1], 2)), as.character(round(t2_ratio_th17_il10_subgroup_L$RD[2], 2)), 
              " ", as.character(round(t2_ratio_th1_th2_subgroup_L$RD[1], 2)), as.character(round(t2_ratio_th1_th2_subgroup_L$RD[2], 2)),  
              " ", as.character(round(t2_ratio_th1_th17_subgroup_L$RD[1], 2)),  as.character(round(t2_ratio_th1_th17_subgroup_L$RD[2], 2)))

#Confidence Interval

#lower bound
#female

cilbtbl12f <- c( "Lower bound" , " ", as.character(round(t2_ratio_il1_il10_subgroup_L$ci.lb[1], 2)), as.character(round(t2_ratio_il1_il10_subgroup_L$ci.lb[2], 2)), 
            " ", as.character(round(t2_ratio_il6_il10_subgroup_L$ci.lb[1], 2)), as.character(round(t2_ratio_il6_il10_subgroup_L$ci.lb[2], 2)),  
            " ", as.character(round(t2_ratio_tnf_il10_subgroup_L$ci.lb[1], 2)), as.character(round(t2_ratio_tnf_il10_subgroup_L$ci.lb[2], 2)), 
            " ", as.character(round(t2_ratio_il12_il10_subgroup_L$ci.lb[1], 2)), as.character(round(t2_ratio_il12_il10_subgroup_L$ci.lb[2], 2)), 
            " ", as.character(round(t2_ratio_ifn_il10_subgroup_L$ci.lb[1], 2)), as.character(round(t2_ratio_ifn_il10_subgroup_L$ci.lb[2], 2)), 
            " ", as.character(round(t2_ratio_il4_il10_subgroup_L$ci.lb[1], 2)), as.character(round(t2_ratio_il4_il10_subgroup_L$ci.lb[2], 2)), 
            " ", as.character(round(t2_ratio_il5_il10_subgroup_L$ci.lb[1], 2)), as.character(round(t2_ratio_il5_il10_subgroup_L$ci.lb[2], 2)),  
            " ", as.character(round(t2_ratio_il13_il10_subgroup_L$ci.lb[1], 2)), as.character(round(t2_ratio_il13_il10_subgroup_L$ci.lb[2], 2)), 
            " ", as.character(round(t2_ratio_il17_il10_subgroup_L$ci.lb[1], 2)), as.character(round(t2_ratio_il17_il10_subgroup_L$ci.lb[2], 2)),  
            " ", as.character(round(t2_ratio_il21_il10_subgroup_L$ci.lb[1], 2)), as.character(round(t2_ratio_il21_il10_subgroup_L$ci.lb[2], 2)),  
            " ", as.character(round(t2_ratio_il2_il10_subgroup_L$ci.lb[1], 2)), as.character(round(t2_ratio_il2_il10_subgroup_L$ci.lb[2], 2)),  
            " ", as.character(round(t2_ratio_gmc_il10_subgroup_L$ci.lb[1], 2)), as.character(round(t2_ratio_gmc_il10_subgroup_L$ci.lb[2], 2)),  
            " ", as.character(round(t2_ratio_il12_il4_subgroup_L$ci.lb[1], 2)), as.character(round(t2_ratio_il12_il4_subgroup_L$ci.lb[2], 2)),  
            " ", as.character(round(t2_ratio_ifn_il4_subgroup_L$ci.lb[1], 2)), as.character(round(t2_ratio_ifn_il4_subgroup_L$ci.lb[2], 2)),  
            " ", as.character(round(t2_ratio_il12_il5_subgroup_L$ci.lb[1], 2)), as.character(round(t2_ratio_il12_il5_subgroup_L$ci.lb[2], 2)),  
            " ", as.character(round(t2_ratio_ifn_il5_subgroup_L$ci.lb[1], 2)), as.character(round(t2_ratio_ifn_il5_subgroup_L$ci.lb[2], 2)), 
            " ", as.character(round(t2_ratio_il12_il13_subgroup_L$ci.lb[1], 2)), as.character(round(t2_ratio_il12_il13_subgroup_L$ci.lb[2], 2)), 
            " ", as.character(round(t2_ratio_ifn_il13_subgroup_L$ci.lb[1], 2)), as.character(round(t2_ratio_ifn_il13_subgroup_L$ci.lb[2], 2)), 
            " ", as.character(round(t2_ratio_il12_il17_subgroup_L$ci.lb[1], 2)), as.character(round(t2_ratio_il12_il17_subgroup_L$ci.lb[2], 2)),  
            " ", as.character(round(t2_ratio_ifn_il17_subgroup_L$ci.lb[1], 2)), as.character(round(t2_ratio_ifn_il17_subgroup_L$ci.lb[2], 2)),  
            " ", as.character(round(t2_ratio_il12_il21_subgroup_L$ci.lb[1], 2)), as.character(round(t2_ratio_il12_il21_subgroup_L$ci.lb[2], 2)), 
            " ", as.character(round(t2_ratio_ifn_il21_subgroup_L$ci.lb[1], 2)), as.character(round(t2_ratio_ifn_il21_subgroup_L$ci.lb[2], 2)),
            " ", as.character(round(t2_ratio_pro_il10_subgroup_L$ci.lb[1], 2)), as.character(round(t2_ratio_pro_il10_subgroup_L$ci.lb[2], 2)),  
            " ", as.character(round(t2_ratio_th1_il10_subgroup_L$ci.lb[1], 2)), as.character(round(t2_ratio_th1_il10_subgroup_L$ci.lb[2], 2)),  
            " ", as.character(round(t2_ratio_th2_il10_subgroup_L$ci.lb[1], 2)), as.character(round(t2_ratio_th2_il10_subgroup_L$ci.lb[2], 2)),    
            " ", as.character(round(t2_ratio_th17_il10_subgroup_L$ci.lb[1], 2)), as.character(round(t2_ratio_th17_il10_subgroup_L$ci.lb[2], 2)), 
            " ", as.character(round(t2_ratio_th1_th2_subgroup_L$ci.lb[1], 2)), as.character(round(t2_ratio_th1_th2_subgroup_L$ci.lb[2], 2)),  
            " ", as.character(round(t2_ratio_th1_th17_subgroup_L$ci.lb[1], 2)),  as.character(round(t2_ratio_th1_th17_subgroup_L$ci.lb[2], 2)))

#male

cilbtbl12m <- c( "Lower bound" , " ", as.character(round(t2_ratio_il1_il10_subgroup_L$ci.lb[1], 2)), as.character(round(t2_ratio_il1_il10_subgroup_L$ci.lb[2], 2)), 
                 " ", as.character(round(t2_ratio_il6_il10_subgroup_L$ci.lb[1], 2)), as.character(round(t2_ratio_il6_il10_subgroup_L$ci.lb[2], 2)),  
                 " ", as.character(round(t2_ratio_tnf_il10_subgroup_L$ci.lb[1], 2)), as.character(round(t2_ratio_tnf_il10_subgroup_L$ci.lb[2], 2)), 
                 " ", as.character(round(t2_ratio_il12_il10_subgroup_L$ci.lb[1], 2)), as.character(round(t2_ratio_il12_il10_subgroup_L$ci.lb[2], 2)), 
                 " ", as.character(round(t2_ratio_ifn_il10_subgroup_L$ci.lb[1], 2)), as.character(round(t2_ratio_ifn_il10_subgroup_L$ci.lb[2], 2)), 
                 " ", as.character(round(t2_ratio_il4_il10_subgroup_L$ci.lb[1], 2)), as.character(round(t2_ratio_il4_il10_subgroup_L$ci.lb[2], 2)), 
                 " ", as.character(round(t2_ratio_il5_il10_subgroup_L$ci.lb[1], 2)), as.character(round(t2_ratio_il5_il10_subgroup_L$ci.lb[2], 2)),  
                 " ", as.character(round(t2_ratio_il13_il10_subgroup_L$ci.lb[1], 2)), as.character(round(t2_ratio_il13_il10_subgroup_L$ci.lb[2], 2)), 
                 " ", as.character(round(t2_ratio_il17_il10_subgroup_L$ci.lb[1], 2)), as.character(round(t2_ratio_il17_il10_subgroup_L$ci.lb[2], 2)),  
                 " ", as.character(round(t2_ratio_il21_il10_subgroup_L$ci.lb[1], 2)), as.character(round(t2_ratio_il21_il10_subgroup_L$ci.lb[2], 2)),  
                 " ", as.character(round(t2_ratio_il2_il10_subgroup_L$ci.lb[1], 2)), as.character(round(t2_ratio_il2_il10_subgroup_L$ci.lb[2], 2)),  
                 " ", as.character(round(t2_ratio_gmc_il10_subgroup_L$ci.lb[1], 2)), as.character(round(t2_ratio_gmc_il10_subgroup_L$ci.lb[2], 2)),  
                 " ", as.character(round(t2_ratio_il12_il4_subgroup_L$ci.lb[1], 2)), as.character(round(t2_ratio_il12_il4_subgroup_L$ci.lb[2], 2)),  
                 " ", as.character(round(t2_ratio_ifn_il4_subgroup_L$ci.lb[1], 2)), as.character(round(t2_ratio_ifn_il4_subgroup_L$ci.lb[2], 2)),  
                 " ", as.character(round(t2_ratio_il12_il5_subgroup_L$ci.lb[1], 2)), as.character(round(t2_ratio_il12_il5_subgroup_L$ci.lb[2], 2)),  
                 " ", as.character(round(t2_ratio_ifn_il5_subgroup_L$ci.lb[1], 2)), as.character(round(t2_ratio_ifn_il5_subgroup_L$ci.lb[2], 2)), 
                 " ", as.character(round(t2_ratio_il12_il13_subgroup_L$ci.lb[1], 2)), as.character(round(t2_ratio_il12_il13_subgroup_L$ci.lb[2], 2)), 
                 " ", as.character(round(t2_ratio_ifn_il13_subgroup_L$ci.lb[1], 2)), as.character(round(t2_ratio_ifn_il13_subgroup_L$ci.lb[2], 2)), 
                 " ", as.character(round(t2_ratio_il12_il17_subgroup_L$ci.lb[1], 2)), as.character(round(t2_ratio_il12_il17_subgroup_L$ci.lb[2], 2)),  
                 " ", as.character(round(t2_ratio_ifn_il17_subgroup_L$ci.lb[1], 2)), as.character(round(t2_ratio_ifn_il17_subgroup_L$ci.lb[2], 2)),  
                 " ", as.character(round(t2_ratio_il12_il21_subgroup_L$ci.lb[1], 2)), as.character(round(t2_ratio_il12_il21_subgroup_L$ci.lb[2], 2)), 
                 " ", as.character(round(t2_ratio_ifn_il21_subgroup_L$ci.lb[1], 2)), as.character(round(t2_ratio_ifn_il21_subgroup_L$ci.lb[2], 2)),
                 " ", as.character(round(t2_ratio_pro_il10_subgroup_L$ci.lb[1], 2)), as.character(round(t2_ratio_pro_il10_subgroup_L$ci.lb[2], 2)),  
                 " ", as.character(round(t2_ratio_th1_il10_subgroup_L$ci.lb[1], 2)), as.character(round(t2_ratio_th1_il10_subgroup_L$ci.lb[2], 2)),  
                 " ", as.character(round(t2_ratio_th2_il10_subgroup_L$ci.lb[1], 2)), as.character(round(t2_ratio_th2_il10_subgroup_L$ci.lb[2], 2)),    
                 " ", as.character(round(t2_ratio_th17_il10_subgroup_L$ci.lb[1], 2)), as.character(round(t2_ratio_th17_il10_subgroup_L$ci.lb[2], 2)), 
                 " ", as.character(round(t2_ratio_th1_th2_subgroup_L$ci.lb[1], 2)), as.character(round(t2_ratio_th1_th2_subgroup_L$ci.lb[2], 2)),  
                 " ", as.character(round(t2_ratio_th1_th17_subgroup_L$ci.lb[1], 2)),  as.character(round(t2_ratio_th1_th17_subgroup_L$ci.lb[2], 2)))

#Upper bound

#Female
ciubtbl12f <- c( "Upper bound" , " ", as.character(round(t2_ratio_il1_il10_subgroup_L$ci.ub[1], 2)), as.character(round(t2_ratio_il1_il10_subgroup_L$ci.ub[2], 2)), 
             " ", as.character(round(t2_ratio_il6_il10_subgroup_L$ci.ub[1], 2)), as.character(round(t2_ratio_il6_il10_subgroup_L$ci.ub[2], 2)),  
             " ", as.character(round(t2_ratio_tnf_il10_subgroup_L$ci.ub[1], 2)), as.character(round(t2_ratio_tnf_il10_subgroup_L$ci.ub[2], 2)), 
             " ", as.character(round(t2_ratio_il12_il10_subgroup_L$ci.ub[1], 2)), as.character(round(t2_ratio_il12_il10_subgroup_L$ci.ub[2], 2)), 
             " ", as.character(round(t2_ratio_ifn_il10_subgroup_L$ci.ub[1], 2)), as.character(round(t2_ratio_ifn_il10_subgroup_L$ci.ub[2], 2)), 
             " ", as.character(round(t2_ratio_il4_il10_subgroup_L$ci.ub[1], 2)), as.character(round(t2_ratio_il4_il10_subgroup_L$ci.ub[2], 2)), 
             " ", as.character(round(t2_ratio_il5_il10_subgroup_L$ci.ub[1], 2)), as.character(round(t2_ratio_il5_il10_subgroup_L$ci.ub[2], 2)),  
             " ", as.character(round(t2_ratio_il13_il10_subgroup_L$ci.ub[1], 2)), as.character(round(t2_ratio_il13_il10_subgroup_L$ci.ub[2], 2)), 
             " ", as.character(round(t2_ratio_il17_il10_subgroup_L$ci.ub[1], 2)), as.character(round(t2_ratio_il17_il10_subgroup_L$ci.ub[2], 2)),  
             " ", as.character(round(t2_ratio_il21_il10_subgroup_L$ci.ub[1], 2)), as.character(round(t2_ratio_il21_il10_subgroup_L$ci.ub[2], 2)),  
             " ", as.character(round(t2_ratio_il2_il10_subgroup_L$ci.ub[1], 2)), as.character(round(t2_ratio_il2_il10_subgroup_L$ci.ub[2], 2)),  
             " ", as.character(round(t2_ratio_gmc_il10_subgroup_L$ci.ub[1], 2)), as.character(round(t2_ratio_gmc_il10_subgroup_L$ci.ub[2], 2)),  
             " ", as.character(round(t2_ratio_il12_il4_subgroup_L$ci.ub[1], 2)), as.character(round(t2_ratio_il12_il4_subgroup_L$ci.ub[2], 2)),  
             " ", as.character(round(t2_ratio_ifn_il4_subgroup_L$ci.ub[1], 2)), as.character(round(t2_ratio_ifn_il4_subgroup_L$ci.ub[2], 2)),  
             " ", as.character(round(t2_ratio_il12_il5_subgroup_L$ci.ub[1], 2)), as.character(round(t2_ratio_il12_il5_subgroup_L$ci.ub[2], 2)),  
             " ", as.character(round(t2_ratio_ifn_il5_subgroup_L$ci.ub[1], 2)), as.character(round(t2_ratio_ifn_il5_subgroup_L$ci.ub[2], 2)), 
             " ", as.character(round(t2_ratio_il12_il13_subgroup_L$ci.ub[1], 2)), as.character(round(t2_ratio_il12_il13_subgroup_L$ci.ub[2], 2)), 
             " ", as.character(round(t2_ratio_ifn_il13_subgroup_L$ci.ub[1], 2)), as.character(round(t2_ratio_ifn_il13_subgroup_L$ci.ub[2], 2)), 
             " ", as.character(round(t2_ratio_il12_il17_subgroup_L$ci.ub[1], 2)), as.character(round(t2_ratio_il12_il17_subgroup_L$ci.ub[2], 2)),  
             " ", as.character(round(t2_ratio_ifn_il17_subgroup_L$ci.ub[1], 2)), as.character(round(t2_ratio_ifn_il17_subgroup_L$ci.ub[2], 2)),  
             " ", as.character(round(t2_ratio_il12_il21_subgroup_L$ci.ub[1], 2)), as.character(round(t2_ratio_il12_il21_subgroup_L$ci.ub[2], 2)), 
             " ", as.character(round(t2_ratio_ifn_il21_subgroup_L$ci.ub[1], 2)), as.character(round(t2_ratio_ifn_il21_subgroup_L$ci.ub[2], 2)),
             " ", as.character(round(t2_ratio_pro_il10_subgroup_L$ci.ub[1], 2)), as.character(round(t2_ratio_pro_il10_subgroup_L$ci.ub[2], 2)),  
             " ", as.character(round(t2_ratio_th1_il10_subgroup_L$ci.ub[1], 2)), as.character(round(t2_ratio_th1_il10_subgroup_L$ci.ub[2], 2)),  
             " ", as.character(round(t2_ratio_th2_il10_subgroup_L$ci.ub[1], 2)), as.character(round(t2_ratio_th2_il10_subgroup_L$ci.ub[2], 2)),    
             " ", as.character(round(t2_ratio_th17_il10_subgroup_L$ci.ub[1], 2)), as.character(round(t2_ratio_th17_il10_subgroup_L$ci.ub[2], 2)), 
             " ", as.character(round(t2_ratio_th1_th2_subgroup_L$ci.ub[1], 2)), as.character(round(t2_ratio_th1_th2_subgroup_L$ci.ub[2], 2)),  
             " ", as.character(round(t2_ratio_th1_th17_subgroup_L$ci.ub[1], 2)),  as.character(round(t2_ratio_th1_th17_subgroup_L$ci.ub[2], 2)))

#Upper bound

#Male
ciubtbl12m <- c( "Upper bound" , " ", as.character(round(t2_ratio_il1_il10_subgroup_L$ci.ub[1], 2)), as.character(round(t2_ratio_il1_il10_subgroup_L$ci.ub[2], 2)), 
                 " ", as.character(round(t2_ratio_il6_il10_subgroup_L$ci.ub[1], 2)), as.character(round(t2_ratio_il6_il10_subgroup_L$ci.ub[2], 2)),  
                 " ", as.character(round(t2_ratio_tnf_il10_subgroup_L$ci.ub[1], 2)), as.character(round(t2_ratio_tnf_il10_subgroup_L$ci.ub[2], 2)), 
                 " ", as.character(round(t2_ratio_il12_il10_subgroup_L$ci.ub[1], 2)), as.character(round(t2_ratio_il12_il10_subgroup_L$ci.ub[2], 2)), 
                 " ", as.character(round(t2_ratio_ifn_il10_subgroup_L$ci.ub[1], 2)), as.character(round(t2_ratio_ifn_il10_subgroup_L$ci.ub[2], 2)), 
                 " ", as.character(round(t2_ratio_il4_il10_subgroup_L$ci.ub[1], 2)), as.character(round(t2_ratio_il4_il10_subgroup_L$ci.ub[2], 2)), 
                 " ", as.character(round(t2_ratio_il5_il10_subgroup_L$ci.ub[1], 2)), as.character(round(t2_ratio_il5_il10_subgroup_L$ci.ub[2], 2)),  
                 " ", as.character(round(t2_ratio_il13_il10_subgroup_L$ci.ub[1], 2)), as.character(round(t2_ratio_il13_il10_subgroup_L$ci.ub[2], 2)), 
                 " ", as.character(round(t2_ratio_il17_il10_subgroup_L$ci.ub[1], 2)), as.character(round(t2_ratio_il17_il10_subgroup_L$ci.ub[2], 2)),  
                 " ", as.character(round(t2_ratio_il21_il10_subgroup_L$ci.ub[1], 2)), as.character(round(t2_ratio_il21_il10_subgroup_L$ci.ub[2], 2)),  
                 " ", as.character(round(t2_ratio_il2_il10_subgroup_L$ci.ub[1], 2)), as.character(round(t2_ratio_il2_il10_subgroup_L$ci.ub[2], 2)),  
                 " ", as.character(round(t2_ratio_gmc_il10_subgroup_L$ci.ub[1], 2)), as.character(round(t2_ratio_gmc_il10_subgroup_L$ci.ub[2], 2)),  
                 " ", as.character(round(t2_ratio_il12_il4_subgroup_L$ci.ub[1], 2)), as.character(round(t2_ratio_il12_il4_subgroup_L$ci.ub[2], 2)),  
                 " ", as.character(round(t2_ratio_ifn_il4_subgroup_L$ci.ub[1], 2)), as.character(round(t2_ratio_ifn_il4_subgroup_L$ci.ub[2], 2)),  
                 " ", as.character(round(t2_ratio_il12_il5_subgroup_L$ci.ub[1], 2)), as.character(round(t2_ratio_il12_il5_subgroup_L$ci.ub[2], 2)),  
                 " ", as.character(round(t2_ratio_ifn_il5_subgroup_L$ci.ub[1], 2)), as.character(round(t2_ratio_ifn_il5_subgroup_L$ci.ub[2], 2)), 
                 " ", as.character(round(t2_ratio_il12_il13_subgroup_L$ci.ub[1], 2)), as.character(round(t2_ratio_il12_il13_subgroup_L$ci.ub[2], 2)), 
                 " ", as.character(round(t2_ratio_ifn_il13_subgroup_L$ci.ub[1], 2)), as.character(round(t2_ratio_ifn_il13_subgroup_L$ci.ub[2], 2)), 
                 " ", as.character(round(t2_ratio_il12_il17_subgroup_L$ci.ub[1], 2)), as.character(round(t2_ratio_il12_il17_subgroup_L$ci.ub[2], 2)),  
                 " ", as.character(round(t2_ratio_ifn_il17_subgroup_L$ci.ub[1], 2)), as.character(round(t2_ratio_ifn_il17_subgroup_L$ci.ub[2], 2)),  
                 " ", as.character(round(t2_ratio_il12_il21_subgroup_L$ci.ub[1], 2)), as.character(round(t2_ratio_il12_il21_subgroup_L$ci.ub[2], 2)), 
                 " ", as.character(round(t2_ratio_ifn_il21_subgroup_L$ci.ub[1], 2)), as.character(round(t2_ratio_ifn_il21_subgroup_L$ci.ub[2], 2)),
                 " ", as.character(round(t2_ratio_pro_il10_subgroup_L$ci.ub[1], 2)), as.character(round(t2_ratio_pro_il10_subgroup_L$ci.ub[2], 2)),  
                 " ", as.character(round(t2_ratio_th1_il10_subgroup_L$ci.ub[1], 2)), as.character(round(t2_ratio_th1_il10_subgroup_L$ci.ub[2], 2)),  
                 " ", as.character(round(t2_ratio_th2_il10_subgroup_L$ci.ub[1], 2)), as.character(round(t2_ratio_th2_il10_subgroup_L$ci.ub[2], 2)),    
                 " ", as.character(round(t2_ratio_th17_il10_subgroup_L$ci.ub[1], 2)), as.character(round(t2_ratio_th17_il10_subgroup_L$ci.ub[2], 2)), 
                 " ", as.character(round(t2_ratio_th1_th2_subgroup_L$ci.ub[1], 2)), as.character(round(t2_ratio_th1_th2_subgroup_L$ci.ub[2], 2)),  
                 " ", as.character(round(t2_ratio_th1_th17_subgroup_L$ci.ub[1], 2)),  as.character(round(t2_ratio_th1_th17_subgroup_L$ci.ub[2], 2)))


# Table 12:
# supptbl12 <- data.table(
#   "Outcome, Arm" = outcometbl12,
#   "N" = Ntbl12, 
#   "Risk difference" = rdtbl12,
#   "Lower bound" = cilbtbl12,
#   "Upper bound" = ciubtbl12
# )

#
supptbl12 <- data.table(
  " " = outcometbl12,
  "Female" =  Ntbl12f, 
  " " = rdtbl12f,
  " " = cilbtbl12f, 
  " " = ciubtbl12f,
  "Male" = Ntbl12m,
  " " = rdtbl12m,
  " " = cilbtbl12m, 
  " " = ciubtbl12m
)


write.csv(supptbl12, file=here("tables/supplementary/immune_supptable12.html"))
print(xtable(supptbl12), type="html", file=here("tables/supplementary/immune_supptable12.html"))



