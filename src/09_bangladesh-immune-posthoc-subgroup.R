#---------------------------------------
# 09_bangladesh-immune-posthoc-subgroup.R
#
# andrew mertens (amertens@berkeley.edu)
#
# Firdaus suggested that we conduct an EMM analysis by pathogen status,
# diarrheal infection, and fever to strengthen the paper and address some
# of the reviewers' points.
#
#---------------------------------------
# TIMING OF THE PATHOGEN MODIFIERS
#
# There are TWO distinct stool collections in this cohort, and they are
# NOT interchangeable:
#
#   (1) 14-month TAC card (multiplex qPCR, 48 targets).
#       Collected essentially the SAME DAY as the T2 blood draw
#       (median gap 0.04 months). Covers protozoa AND STH.
#       Source: eed-substudy-data/sth data/14mo qPCR Data/
#
#   (2) ~30-month year-2 parasite specimen. ONE specimen assayed BOTH by
#       Kato-Katz (Ascaris/Trichuris/hookworm egg counts) AND by
#       qPCR (Giardia/Cryptosporidium/E. histolytica).
#       Kato-Katz and protozoan qPCR are the SAME sample at the SAME
#       visit - they are not two timepoints.
#       Median age 30.3 mo (IQR 29.2-31.6) in the Control/N+WSH sample.
#
# Relative to the immune blood draws (Control + N+WSH, N=588):
#       parasite specimen - T2 blood = +16.0 months (IQR 14.8-17.2)
#       parasite specimen - T3 blood =  +2.0 months (IQR  1.2- 2.7)
# i.e. the 30-month specimen post-dates BOTH draws in 100% of children.
#
# CONSEQUENCE:
#   T2    -> use the 14-month TAC modifiers (contemporaneous).
#            The 30-month specimen is 16 months POST-outcome and must
#            not be used as a T2 effect modifier.
#   T3    -> use the 30-month specimen (near-contemporaneous, +2 mo).
#   Delta -> use the 14-month TAC modifiers: they precede the T3
#            endpoint, so they modify subsequent change. The 30-month
#            specimen post-dates both endpoints of the change score.
#
# WARNING about clean_wbb_sth_pathogens.rds: its agedays/svyweek/svyyear
# columns are NOT the parasite-visit age - they are inherited from the
# diarrhoea-pathogen analysis and run ~7.4 months too young (median 22.3
# vs 30.2 mo; svyyear 2012-13 vs the 2015 collection). The pathogen
# VALUES are correct (verified identical to
# washb-bangladesh-sth-public.csv), but the age must be taken from the
# public file. We attach it below as agem_parasite and never use the
# stale field.
#---------------------------------------

#---------------------------------------
# preamble
#---------------------------------------

rm(list=ls())
source(here::here("0-config.R"))

## Minimum cell count required in EVERY arm x modifier-level cell for a
## modifier to be reported at a given timepoint. Modifiers failing this
## are dropped and LOGGED (never silently). Set to 0 to disable.
MIN_CELL <- 10

eedDir <- "C:/Users/andre/OneDrive/Documents/washb_substudies/eed-substudy-data/"

dfull <- read.csv(paste0(eedDir, "bangladesh-merged-eed-mn-data.csv"))

Wvars<-c("monsoon_bt2","ageday_bt2", "sex","birthord", "momage","momheight","momedu", "hfiacat", "Nlt18","Ncomp", "watmin", "walls", "floor", "elec", "asset_wardrobe", "asset_table", "asset_chair", "asset_clock","asset_khat", "asset_chouki", "asset_radio", "asset_tv", "asset_refrig", "asset_bike", "asset_moto", "asset_sewmach", "asset_mobile", "n_cattle", "n_goat", "n_chicken")

d_illness <- dfull %>%
  select(childid, diar7d_t2, diar7d_t3, fever7d_t2, fever7d_t3) %>%
  mutate(illness_data = 1)

## Public <-> real ID crosswalk (1:1 over 5551 compounds)
IDs <- read.csv(paste0(eedDir, "public-ids.csv"))
id_map <- IDs %>% distinct(dataid_r, dataid)


#---------------------------------------
# (A) ~30-month year-2 parasite specimen: Kato-Katz + protozoan qPCR
#---------------------------------------

wbb_sth <- readRDS(paste0(eedDir, "sth data/clean_wbb_sth_pathogens.rds")) %>%
  rename(block_r = block, clusterid_r = clusterid, dataid_r = dataid) %>%
  left_join(IDs, by = c("block_r", "clusterid_r", "dataid_r"))

## TRUE age at the parasite stool collection, taken from the public file
## (the age carried in clean_wbb_sth_pathogens.rds is the wrong visit).
sth_pub_age <- read.csv(paste0(eedDir, "sth data/washb-bangladesh-sth-public.csv")) %>%
  filter(personid %in% c("T1", "T2")) %>%
  select(dataid_r = dataid, childid = personid, agem_parasite = agem)

ch <- wbb_sth %>%
  left_join(sth_pub_age, by = c("dataid_r", "childid")) %>%
  select(dataid, childid, agem_parasite,
         logalepg, logttepg,          # Kato-Katz egg counts
         posgi, poseh, poscr) %>%     # protozoan qPCR, same specimen
  filter(!is.na(logalepg) | !is.na(logttepg) |
         !is.na(posgi) | !is.na(poseh) | !is.na(poscr)) %>%
  rename(ch_pos_giardia   = posgi,
         ch_pos_entamoeba = poseh,
         ch_pos_crypto    = poscr) %>%
  mutate(ch_pos_ascaris   = 1 * (logalepg > 0),
         ch_pos_trichuris = 1 * (logttepg > 0)) %>%
  select(-logalepg, -logttepg)

ch$childid <- ch$dataid * 10 + as.numeric(gsub("T", "", ch$childid))
ch <- ch %>% select(-dataid) %>% mutate(pathogen_data = 1)


#---------------------------------------
# (B) 14-month TAC multiplex qPCR (contemporaneous with the T2 draw)
#---------------------------------------
## df.cov.rds and pathogen_log10_quantity.rds are keyed on PUBLIC ids
## (dataid_r * 10 + child number); the immune dataset uses real ids, so
## the crosswalk is required. log10 quantity: 0 = negative, >0 = positive,
## NA = not assayed.

tac_cov <- readRDS(paste0(eedDir, "sth data/14mo qPCR Data/df.cov.rds")) %>%
  select(childid_pub = childid, dataid_pub = dataid, st_agem2) %>%
  mutate(childid_pub = as.character(childid_pub),
         dataid_r    = as.numeric(as.character(dataid_pub)),
         childno     = as.numeric(substr(childid_pub, nchar(childid_pub), nchar(childid_pub)))) %>%
  left_join(id_map, by = "dataid_r") %>%
  mutate(childid = dataid * 10 + childno) %>%
  select(childid_pub, childid, agem_tac = st_agem2)

tac_q <- readRDS(paste0(eedDir, "sth data/14mo qPCR Data/pathogen_log10_quantity.rds")) %>%
  mutate(childid_pub = as.character(childid)) %>%
  select(childid_pub,
         tac_giardia   = Giardia,
         tac_crypto    = Cryptosporidium,
         tac_ehist     = `E.histolytica`,
         tac_entamoeba = `pan Entamoeba`,
         tac_ascaris   = Ascaris,
         tac_trichuris = Trichuris,
         tac_ancy      = Ancyclostoma,
         tac_necator   = Necator)

tac <- tac_cov %>%
  inner_join(tac_q, by = "childid_pub") %>%
  mutate(tac14_pos_giardia   = 1 * (tac_giardia   > 0),
         tac14_pos_crypto    = 1 * (tac_crypto    > 0),
         tac14_pos_entamoeba = 1 * (tac_ehist > 0 | tac_entamoeba > 0),
         tac14_pos_ascaris   = 1 * (tac_ascaris   > 0),
         tac14_pos_trichuris = 1 * (tac_trichuris > 0),
         tac14_pos_hookworm  = 1 * (tac_ancy > 0 | tac_necator > 0),
         tac_data = 1) %>%
  select(childid, agem_tac, tac_data, starts_with("tac14_pos_"))


#---------------------------------------
# Load the analysis dataset and assemble
#---------------------------------------
## NOTE (change from the previous version): the base sample is now the
## full immune analysis dataset. Previously every EMM - including the
## diarrhoea and fever ones - was silently restricted to children with
## 30-month parasite data, which dropped ~30% of the cohort for no
## reason. Each EMM is now subset only on its own modifier + outcome.

EE_df <- readRDS(paste0(dropboxDir,"Data/Cleaned/Audrie/bangladesh-immune-analysis-dataset.rds"))

d <- EE_df %>%
  left_join(ch,        by = "childid") %>%
  left_join(tac,       by = "childid") %>%
  left_join(d_illness, by = "childid")

cat("\n--- assembled analysis dataset ---\n")
cat("immune analysis dataset N :", nrow(EE_df), "\n")
cat("with 30-mo parasite data  :", sum(!is.na(d$pathogen_data)), "\n")
cat("with 14-mo TAC data       :", sum(!is.na(d$tac_data)), "\n")
cat("with illness (diar/fever) :", sum(!is.na(d$illness_data)), "\n")

## Timing audit, restricted to the two arms actually contrasted.
aud <- d %>% filter(tr %in% c("Control", "Nutrition + WSH"))
qi <- function(x) sprintf("median %.1f (IQR %.1f-%.1f)",
                          median(x, na.rm = TRUE),
                          quantile(x, .25, na.rm = TRUE),
                          quantile(x, .75, na.rm = TRUE))
cat("\n--- age audit (Control + N+WSH, months) ---\n")
cat("T2 blood draw            :", qi(aud$ageday_bt2 / 30.4167), "\n")
cat("T3 blood draw            :", qi(aud$ageday_bt3 / 30.4167), "\n")
cat("14-mo TAC stool          :", qi(aud$agem_tac), "\n")
cat("30-mo parasite specimen  :", qi(aud$agem_parasite), "\n")
cat("TAC stool - T2 draw      :", qi(aud$agem_tac - aud$ageday_bt2 / 30.4167), "\n")
cat("parasite  - T2 draw      :", qi(aud$agem_parasite - aud$ageday_bt2 / 30.4167), "\n")
cat("parasite  - T3 draw      :", qi(aud$agem_parasite - aud$ageday_bt3 / 30.4167), "\n")


#---------------------------------------
# Modifier coding
#---------------------------------------
## NOTE: qPCR-based STH from the 30-month specimen (qpcr.positive.Al /
## .Tt) remains excluded. That substudy was run only in the Control and
## WSH arms, while the immune outcomes exist only in Control and
## Nutrition+WSH - so it has no overlap with the intervention arm and
## the EMM contrast is not estimable. Kato-Katz was run in all arms.

path30_vars <- c("ch_pos_giardia", "ch_pos_entamoeba", "ch_pos_crypto",
                 "ch_pos_trichuris", "ch_pos_ascaris")
tac14_vars  <- c("tac14_pos_giardia", "tac14_pos_crypto", "tac14_pos_entamoeba",
                 "tac14_pos_ascaris", "tac14_pos_trichuris", "tac14_pos_hookworm")
illness_vars <- c("diar7d_t2", "diar7d_t3", "fever7d_t2", "fever7d_t3")

## Code every modifier as a 3-level factor with an explicit "Missing"
## level; washb_function() then drops "Missing" before fitting.
code_modifier <- function(x) {
  x <- as.character(x)
  x[is.na(x)] <- "Missing"
  factor(x, levels = c("0", "1", "Missing"))
}
for (v in c(path30_vars, tac14_vars, illness_vars)) {
  d[[v]] <- code_modifier(d[[v]])
}

d$monsoon_bt3 <- as.factor(d$monsoon_bt3)
d$monsoon_bt3 <- addNA(d$monsoon_bt3)
levels(d$monsoon_bt3)[length(levels(d$monsoon_bt3))] <- "Missing"


#---------------------------------------
# EMM fitting
#---------------------------------------

washb_function <- function(df, Y, v) {

  ## Restrict to children with KNOWN (non-missing) modifier status so
  ## the modifier enters the model as a binary 0/1 factor. Keeping the
  ## "Missing" level made V a 3-level factor, which corrupted two things:
  ##   1. lincom returned a third (Missing) stratum, so lincom[1:2,] no
  ##      longer reliably mapped to the 0 and 1 subgroups.
  ##   2. The tr:V interaction term shifted off the hard-coded fit row.
  ## Subsetting also makes the GLM sample match the N/means in the table.
  dsub <- df %>%
    filter(!is.na(!!sym(Y)), !!sym(v) != "Missing") %>%
    mutate(!!sym(v) := droplevels(factor(!!sym(v))))

  W <- subset(dsub, select = c(Wvars, v))

  mean_res <- dsub %>% group_by(tr, !!sym(v)) %>%
    summarise(n = n(), mean = mean(!!sym(Y)), sd = sd(!!sym(Y)), .groups = "drop")
  colnames(mean_res)[2] <- "subgroup"
  mean_res$subgroup <- as.character(mean_res$subgroup)

  mean_res <- mean_res %>%
    pivot_wider(id_cols = subgroup, names_from = tr,
                values_from = c(n, mean, sd), names_glue = "{.value}_{tr}")

  colnames(mean_res) <- gsub("Nutrition \\+ WSH", "NWSH", colnames(mean_res))

  temp <- washb_glm(Y = dsub[, Y], tr = dsub$tr, pair = NULL, W = W, V = v,
                    id = dsub$block, contrast = c("Control", "Nutrition + WSH"),
                    family = "gaussian", print = FALSE, verbose = FALSE)
  temp_metric <- as.data.frame(temp$lincom[1:2, ])
  colnames(temp_metric) <- c("subgroup", "RD", "SE", "ci.lb", "ci.ub", "z", "P-value")
  temp_metric$subgroup <- as.character(temp_metric$subgroup)

  ## Interaction p = the single tr:V interaction coefficient. With a
  ## binary modifier there is exactly one ":" term in the fit table;
  ## extract it by name rather than by a fixed row position.
  fit <- as.data.frame(temp$fit)
  int_row <- grep(":", rownames(fit))
  int_p <- if (length(int_row) >= 1) fit[int_row[1], 6] else NA_real_
  temp_metric$int_Pval <- c(int_p, NA)
  temp_metric$Y <- Y
  temp_metric$V <- v
  temp_metric <- left_join(mean_res, temp_metric, by = "subgroup")
  return(temp_metric)
}

## Ratio outcomes at each timepoint. Delta uses the d23_ prefix
## (same construction as src/06).
ratio_stems <- c(
  "ratio_th1_th2", "ratio_il12_il4", "ratio_ifn_il4", "ratio_il12_il5", "ratio_ifn_il5",  "ratio_il12_il13", "ratio_ifn_il13",
  "ratio_pro_il10", "ratio_il1_il10","ratio_il6_il10", "ratio_tnf_il10",
  "ratio_il2_il10",  "ratio_th1_il10",  "ratio_th2_il10", "ratio_il12_il10", "ratio_ifn_il10", "ratio_il4_il10",
  "ratio_il5_il10", "ratio_il13_il10","ratio_th17_il10", "ratio_il17_il10", "ratio_il21_il10",
  "ratio_gmc_il10",
  "ratio_th1_th17", "ratio_il12_il17", "ratio_ifn_il17", "ratio_il12_il21", "ratio_ifn_il21")

t2_Y    <- paste0("t2_",  ratio_stems)
t3_Y    <- paste0("t3_",  ratio_stems)
delta_Y <- paste0("d23_", ratio_stems)

## Modifier sets, assigned on TEMPORALITY (see header).
## T3 carries BOTH pathogen families on purpose:
##   - the 14-month TAC modifiers are measured a median of 13.9 months
##     BEFORE the 28-month draw (100% of children), giving a prospective
##     "prior infection" contrast;
##   - the ~30-month specimen is measured 2.0 months AFTER it, giving a
##     cross-sectional "concurrent infection" contrast.
## Reporting both lets a reader see whether the two agree. Table labels
## and footnotes state which is which - see modifier_heading()/
## modifier_note() in table scripts/_emm_helpers.R.
Vvars_t2    <- c("diar7d_t2", "fever7d_t2", tac14_vars)
Vvars_t3    <- c("diar7d_t2", "diar7d_t3", "fever7d_t2", "fever7d_t3",
                 tac14_vars, path30_vars)
Vvars_delta <- c("diar7d_t2", "diar7d_t3", "fever7d_t2", "fever7d_t3", tac14_vars)

#---------------------------------------
# Sparse-cell screen
#---------------------------------------
## Several modifiers have almost no positives in the two contrasted
## arms. At 30 months, Entamoeba histolytica has 0 positives in Control
## and 1 in N+WSH, and Cryptosporidium has 5 and 5. Fitting these still
## returns RDs, CIs and p-values - some nominally significant - but they
## are artefacts of a single child. Screen them out up front and print
## exactly what was dropped and why.

screen_modifiers <- function(Y_set, Vvars_set, label) {
  Y_present <- Y_set[Y_set %in% names(d)]
  has_out <- rowSums(!is.na(d[, Y_present, drop = FALSE])) > 0
  dd <- d[has_out & d$tr %in% c("Control", "Nutrition + WSH"), , drop = FALSE]
  dd$tr <- droplevels(factor(dd$tr))

  keep <- character(0)
  cat("\n--- modifier screen:", label, "(min cell =", MIN_CELL, ") ---\n")
  for (v in Vvars_set) {
    tb <- table(dd$tr, dd[[v]])
    tb <- tb[, colnames(tb) != "Missing", drop = FALSE]
    if (ncol(tb) < 2) {
      cat(sprintf("  DROP %-22s only one observed level\n", v))
      next
    }
    if (min(tb) < MIN_CELL) {
      cat(sprintf("  DROP %-22s min cell = %d  [C: %d neg / %d pos ; N+WSH: %d neg / %d pos]\n",
                  v, min(tb), tb[1,1], tb[1,2], tb[2,1], tb[2,2]))
      next
    }
    cat(sprintf("  KEEP %-22s min cell = %d  [C: %d neg / %d pos ; N+WSH: %d neg / %d pos]\n",
                v, min(tb), tb[1,1], tb[1,2], tb[2,1], tb[2,2]))
    keep <- c(keep, v)
  }
  keep
}

Vvars_t2    <- screen_modifiers(t2_Y,    Vvars_t2,    "T2 (14 months)")
Vvars_t3    <- screen_modifiers(t3_Y,    Vvars_t3,    "T3 (28 months)")
Vvars_delta <- screen_modifiers(delta_Y, Vvars_delta, "Delta (14 -> 28 months)")


run_emm_set <- function(Y_set, Vvars_set, label) {
  out <- NULL
  for (i in Y_set) {
    if (!i %in% names(d)) {
      message("Skipping ", i, " (not in dataset)")
      next
    }
    for (j in Vvars_set) {
      cat(label, ":", i, "x", j, "\n")
      temp_res <- tryCatch(
        washb_function(d, i, j),
        error = function(e) {
          warning("washb_function failed for ", i, " x ", j, ": ", conditionMessage(e))
          NULL
        }
      )
      out <- bind_rows(out, temp_res)
    }
  }
  out
}

res_t2    <- run_emm_set(t2_Y,    Vvars_t2,    "T2")
res_t3    <- run_emm_set(t3_Y,    Vvars_t3,    "T3")
res_delta <- run_emm_set(delta_Y, Vvars_delta, "Delta")

write.csv(res_t2,    here("results/bangladesh-immune-posthoc-subgroup-results-t2.csv"),    row.names = FALSE)
write.csv(res_t3,    here("results/bangladesh-immune-posthoc-subgroup-results-t3.csv"),    row.names = FALSE)
write.csv(res_delta, here("results/bangladesh-immune-posthoc-subgroup-results-delta.csv"), row.names = FALSE)
