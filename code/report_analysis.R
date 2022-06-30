## ------------------------------------------------------------------------------------------
##
## Script name: SOT-sepsis main analysis
##
## Purpose of script: Show an example of code to clean lmtp analysis results and compute
##      confidence intervals for effect modifier estimates.
## 
## Author: Katherine Hoffman
##
## Date Created: 2022-05-22
##
## Author: Katherine Hoffman, 2022
## Email: kathoffman.stats@gmail.com
##
## ------------------------------------------------------------------------------------------

library(lmtp)
library(tidyverse)
library(gt)

dat_lmtp <- read_rds(here::here("data/toy_data.rds")) 
dat_dx_nosep <- dat_lmtp %>% filter(transplant_type_proc != 1) # procedure transplants removed for sensitivity
dat_trans_immuno <- dat_lmtp %>% filter(immunosuppressed_no_sot_not_immunosuppressed != 1)

# read in cluster results

load(here::here("output/results.rdata"))

report_contrast_results <- function(dat_copy, ts, tns, nts, ntns, transplant_label, no_transplant_label){
  transplant_fit <- lmtp_contrast(ts, ref = tns)
  no_transplant_fit <- lmtp_contrast(nts, ref = ntns)
  
  cis <- data.frame(pop = c(transplant_label, no_transplant_label),
                    shift_low = c(ts$low, nts$low),
                    shift_high = c(ts$high, nts$high),
                    ref_low = c(tns$low, ntns$low),
                    ref_high = c(tns$high, ntns$high)
  )
  
  # influence curves for transplant and non transplants ---------------------
  
  ict <- ts$eif -  tns$eif
  icnot <- nts$eif -  ntns$eif
  
  # replace each obs w respective influence curve ---------------------------
  
  dat_copy[dat_copy$transplant_any == 1, 'ict'] <- ict # all transplant patients assigned their influence curve
  dat_copy[is.na(dat_copy$ict), 'ict'] <- 0 # patients without transplant have ic of 0
  dat_copy[, 'ict'] <- dat_copy$transplant_any * dat_copy[, 'ict'] / 
    mean(dat_copy$transplant_any)
  # transplant patients influence curve becomes their ic
  #   divided by the proportion of patients with a transplant in whole population
  
  # same thing for no transplant patients
  dat_copy[dat_copy$transplant_any == 0, 'icnot'] <- icnot
  dat_copy[is.na(dat_copy$icnot), 'icnot'] <- 0
  dat_copy[, 'icnot'] <- (1 - dat_copy$transplant_any) * dat_copy[, 'icnot'] /
    mean(1 - dat_copy$transplant_any)
  
  set <- sd(dat_copy$ict) / sqrt(nrow(dat_copy))
  senot <- sd(dat_copy$icnot) / sqrt(nrow(dat_copy))
  
  seeffmod <- sd(with(dat_copy, ict - icnot)) / sqrt(nrow(dat_copy))
  
  transplant_contrast <- lmtp_contrast(ts, ref = tns)
  no_transplant_contrast <- lmtp_contrast(nts, ref = ntns)
  
  
  effectt <- transplant_contrast$vals$theta
  effectnot <- no_transplant_contrast$vals$theta
  
  effectmod <- effectnot - effectt
  
  effects_lmtp_form <- data.frame("theta" = effectmod, "std.error" = seeffmod) %>%
    mutate(conf.low = theta + 1.96*std.error, conf.high=theta - 1.96*std.error,
           p.value = 2*pnorm(abs(theta/std.error), lower.tail = FALSE))
  
  return(list("transplant_fit" = transplant_fit, "no_transplant_fit" = no_transplant_fit,
              "cis" = cis, "effects_lmtp_form" = effects_lmtp_form,
              "transplant_label" = transplant_label, "no_transplant_label" = no_transplant_label))
  
}

results_tbl <- function(lmtp_form){
  
  
  tbl <-
    bind_rows(lmtp_form$transplant_fit$vals, lmtp_form$no_transplant_fit$vals) %>%
    mutate(pop = c(lmtp_form$transplant_label, lmtp_form$no_transplant_label)) %>%
    left_join(lmtp_form$cis) %>%
    full_join(lmtp_form$effects_lmtp_form)
  
  tbl %>%
    select(pop, shift, ref, theta, conf.low, conf.high, shift_low, shift_high, ref_low, ref_high, p.value) %>%
    mutate(p.value = case_when(p.value < .001 ~ "<0.001",
                               T ~ as.character(round(p.value, 2))),
           across(shift:ref_high, ~ round(.x * 100, 1))#,
           # across(everything(), ~replace_na(.x, "-"))
    ) %>%
    gt() %>%
    cols_label(pop = "Population",
               theta = "Difference",
               shift = "Sepsis",
               ref = "No Sepsis",
               p.value = "p-value") %>%
    tab_spanner("28 day Survival Probability",
                columns = c(shift, ref, theta, conf.low, conf.high)) %>%
    cols_merge(
      columns = c("theta", "conf.low", "conf.high"),
      pattern = "{1} ({2}, {3})"
    ) %>%
    cols_merge(
      columns = c("shift", "shift_low", "shift_high"),
      pattern = "{1} ({2}, {3})"
    ) %>%
    cols_merge(
      columns = c("ref", "ref_low", "ref_high"),
      pattern = "{1} ({2}, {3})"
    ) %>%
    tab_options(table.font.size = 'small', data_row.padding = gt::px(1))
}

# Main adjusted analysis
report_contrast_results(dat_lmtp, transplant_sepsis, transplant_no_sepsis,
                        no_transplant_sepsis, no_transplant_no_sepsis,
                        "SOT", "No SOT")  %>%
  results_tbl() %>%
  tab_header("Main Analysis - Adjusted")

report_contrast_results(dat_dx_nosep, dx_transplant_sepsis, dx_transplant_no_sepsis,
                        no_transplant_sepsis, no_transplant_no_sepsis,
                        "Dx SOT", "No SOT") %>%
  results_tbl() %>%
  tab_header("Sensitivity Analysis 1 (No Procedure) - Adjusted")

report_contrast_results(dat_trans_immuno, transplant_sepsis, transplant_no_sepsis,
                        immuno_sepsis, immuno_no_sepsis,
                        "SOT", "Immunocompromised") %>%
  results_tbl() %>%
  tab_header("Sensitivity Analysis 1 (Immunocompromised Controls) - Adjusted")
