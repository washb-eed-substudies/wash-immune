writeqntle<-function(vector) {
  quantiles<-round(quantile(vector, na.rm=TRUE), 2)
  paste(quantiles[3], " (", quantiles[2], ", ", quantiles[4], ")", sep="")
}

outcome<-c("Outcome", "IL-1β (pg/ml)", "Il-6 (pg/ml)", "TNF-α (pg/ml)", "CRP (mg/L)", "IL-12 (pg/ml)",
             "IFN-γ (pg/ml)", "IL-4 (pg/ml)", "IL-5 (pg/ml)", "IL-13 (pg/ml)", "IL-17A (pg/ml)", 
             "IL-21 (pg/ml)", "IL-10 (pg/ml)", "IL-2 (pg/ml)", "GM-CSF (pg/ml)", "AGP (g/L)", "IGF-1 (μg/L)")

t2<-c("Median (25th, 75th percentile)", writeqntle(lab$il1_t2), writeqntle(lab$il6_t2), writeqntle(lab$tnfa_t2), writeqntle(lab$crp_t2), writeqntle(lab$il12_t2), writeqntle(lab$ifng_t2), 
      writeqntle(lab$il4_t2), writeqntle(lab$il5_t2), writeqntle(lab$il13_t2), writeqntle(lab$il17_t2), writeqntle(lab$il21_t2), writeqntle(lab$il10_t2), writeqntle(lab$il2_t2), 
      writeqntle(lab$gmcsf_t2), writeqntle(lab$agp_t2), writeqntle(lab$igf_t2))

t3<-c("Median (25th, 75th percentile)", writeqntle(lab$il1_t3), writeqntle(lab$il6_t3), writeqntle(lab$tnfa_t3), " ", writeqntle(lab$il12_t3), writeqntle(lab$ifng_t3), 
      writeqntle(lab$il4_t3), writeqntle(lab$il5_t3), writeqntle(lab$il13_t3), writeqntle(lab$il17_t3), writeqntle(lab$il21_t3), writeqntle(lab$il10_t3), writeqntle(lab$il2_t3), 
      writeqntle(lab$gmcsf_t3), " ", writeqntle(lab$igf_t3))

tbl<-data.table(" "=outcomes3,
                  "Child Age 14 Months"=t2,
                  "Child Age 28 Months"=t3)

write.csv(tbl, file=here('tables/archived/cytokines_median.csv'))
print(xtable(tbl), type="html", file=here("tables/archived/cytokines_median.html"))



