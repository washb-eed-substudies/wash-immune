rm(list=ls())
source(here::here("0-config.R"))

load(here("results/crp_agp_n_wsh/unadj_crp_agp_t3.RData"))
unadj_agp <- temp_metric_agp
unadj_crp <- temp_metric_crp

load(here("results/crp_agp_n_wsh/age_sex_adj_crp_agp_t3.RData"))
age_sex_agp <- temp_metric_agp
age_sex_crp <- temp_metric_crp

load(here("results/crp_agp_n_wsh/adj_crp_agp_t3.RData"))
adj_agp <- temp_metric_agp
adj_crp <- temp_metric_crp

new_agp <- cbind(unadj_agp, age_sex_agp, adj_agp)
colnames(new_agp) <- c("RD unadj", "ci.lb unadj", "ci.ub unadj", "SE unadj", "z unadj", "Pval unadj",
                      "RD age-sex adj", "ci.lb age-sex adj", "ci.ub age-sex adj", "SE age-sex adj", "z age-sex adj", "Pval age-sex adj",
                      "RD adj", "ci.lb adj", "ci.ub adj", "SE adj", "z adj", "Pval adj")
new_agp <- as.data.frame(new_agp)
agp <- new_agp %>% select(!grep("SE|z", names(new_agp))) %>% round(3)

new_crp <- cbind(unadj_crp, age_sex_crp, adj_crp)
colnames(new_crp) <- c("RD unadj", "ci.lb unadj", "ci.ub unadj", "SE unadj", "z unadj", "Pval unadj",
                       "RD age-sex adj", "ci.lb age-sex adj", "ci.ub age-sex adj", "SE age-sex adj", "z age-sex adj", "Pval age-sex adj",
                       "RD adj", "ci.lb adj", "ci.ub adj", "SE adj", "z adj", "Pval adj")
new_crp <- as.data.frame(new_crp)
crp <- new_crp %>% select(!grep("SE|z", names(new_crp))) %>% round(3)

write.csv(agp, here('tables/agp_t3_allarms.csv'))
write.csv(crp, here('tables/crp_t3_allarms.csv'))
