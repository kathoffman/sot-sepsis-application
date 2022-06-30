## ------------------------------------------------------------------------------------------
##
## Script name: SOT-sepsis main analysis
##
## Purpose of script: Show an example of code to run an analysis for a 
##     time-to-event outcome (28 day mortality) with dozens of baseline confounders.
##     The exposure of interest in sepsis at the time of ICU admission and the effect
##     modifier is Solid Organ Transplant (SOT) status. Patients are lost to follow-up
##     via discharge from the hospital.
## 
## Author: Katherine Hoffman
##
## Date Created: 2022-05-22
##
## Author: Katherine Hoffman, 2022
## Email: kathoffman.stats@gmail.com
##
## ------------------------------------------------------------------------------------------
##
## Notes: The data set of n=5000 patients contains realistic values for patients but
## is a toy data set and cannot be used to obtain the result in the paper due to
## data-sharing restrictions. The contrasts are sepsis (yes vs. no) in two subgroups (SOT vs. no SOT)
##
## ------------------------------------------------------------------------------------------

## system set up -------------------------------------------------------------------------

set.seed(7)

# install.packages(c("tidyverse","earth","BART","glmnet"))
# devtools::install_github("nt-williams/lmtp@sl3") -- sl3 compatible branch (faster)
# devtools::install_github("tlverse/sl3")

library(sl3) 
library(lmtp) 
library(tidyverse)

## load data and define column names for lmtp  -------------------------------------------------------------------------

dat_lmtp <- read_rds(here::here("data/toy_data.rds"))

A <- "sep_any" # exposure indicator
Y <- dat_lmtp %>% select(starts_with("Y_")) %>% names() # outcome indicators
C <- dat_lmtp %>% select(starts_with("C_")) %>% names() # censoring indicators

# separate confounders by transplant/sepsis status
W_t <- dat_lmtp %>% select(age:micro_any_true, -transplant_type_no_sot) %>% names()
W_tp <- dat_lmtp %>% select(age:micro_any_true, -transplant_type_no_sot, -starts_with("dx_transplant")) %>% names() # proc transplant
W_td <- dat_lmtp %>% select(age:micro_any_true, -transplant_type_no_sot, -starts_with("proc_transplant")) %>% names() # dx transplant
W_nt <- dat_lmtp %>% select(age:micro_any_true, -contains("transplant")) %>% names()

## set up superlearner libraries  -------------------------------------------------------------------------

mars_grid_params <- list( # set up Earth grid
  degree = c(2,3),
  penalty = c(1,2,3)
)

mars_grid <- expand.grid(mars_grid_params, KEEP.OUT.ATTRS = FALSE)
mars_learners <- apply(mars_grid, MARGIN = 1, function(tuning_params) {
  do.call(Lrnr_earth$new, as.list(tuning_params))
})

lrn_lasso <- Lrnr_glmnet$new(alpha = 1)
lrn_ridge <- Lrnr_glmnet$new(alpha = 0)
lrn_enet <- Lrnr_glmnet$new(alpha = 0.5)
lrn_xgboost <- Lrnr_xgboost$new()
lrn_mean <- Lrnr_mean$new()

# For this demo repo, only LASSO and mean are run to speed up computational time
learners <- unlist(list(
  # mars_learners,
  lrn_lasso,
  # lrn_ridge,
  # lrn_enet,
  # lrn_xgboost,
  lrn_mean
),
recursive = TRUE
)

lrnrs <- make_learner(Stack, learners)

## set up data subsets for main and sensitivity analyses  -------------------------------------------------------------------------

dat_trans <- dat_lmtp %>% filter(transplant_any == 1) 
dat_no_trans <- dat_lmtp %>% filter(transplant_any == 0)

dat_proc <- dat_lmtp %>% filter(transplant_type_proc == 1)
dat_dx <- dat_lmtp %>% filter(transplant_type_no_sot == 0 & transplant_type_proc == 0)

dat_immuno <- dat_lmtp %>% filter(immunosuppressed_no_sot_immunosuppressed == 1)
dat_no_immuno <- dat_lmtp %>% filter(immunosuppressed_no_sot_not_immunosuppressed == 1)

transplant_sepsis <-
    lmtp_tmle(dat_trans,
              trt = A,
              outcome = Y,
              baseline = W_t,
              cens = C,
              outcome_type = "survival",
              folds = 5,
              .SL_folds = 5,
              shift = static_binary_on,
              learners_outcome = lrnrs,
              learners_trt = lrnrs)

transplant_no_sepsis <-
    lmtp_tmle(dat_trans,
              trt = A,
              outcome = Y,
              baseline = W_t,
              cens = C,
              outcome_type = "survival",
              folds = 5,
              .SL_folds = 5,
              shift = static_binary_off,
              learners_outcome = lrnrs,
              learners_trt = lrnrs)

no_transplant_sepsis <- 
    lmtp_tmle(dat_no_trans,
              trt = A,
              outcome = Y,
              baseline = W_nt, 
              cens = C,
              outcome_type = "survival",
              folds = 5,
              .SL_folds = 5,
              shift = static_binary_on,
              learners_outcome = lrnrs,
              learners_trt = lrnrs)


no_transplant_no_sepsis <-
    lmtp_tmle(dat_no_trans,
              trt = A,
              outcome = Y,
              baseline = W_nt,
              cens = C,
              outcome_type = "survival",
              folds = 5,
              .SL_folds = 5,
              shift = static_binary_off,
              learners_outcome = lrnrs,
              learners_trt = lrnrs)


dx_transplant_sepsis <- lmtp_tmle(dat_dx,
              trt = A,
              outcome = Y,
              baseline = W_td,
              cens = C,
              outcome_type = "survival",
              folds = 5,
              .SL_folds = 5,
              shift = static_binary_on,
              learners_outcome = lrnrs,
              learners_trt = lrnrs)

dx_transplant_no_sepsis <-
    lmtp_tmle(dat_dx,
              trt = A,
              outcome = Y,
              baseline = W_td,
              cens = C,
              outcome_type = "survival",
              folds = 5,
              .SL_folds = 5,
              shift = static_binary_off,
              learners_outcome = lrnrs,
              learners_trt = lrnrs)


immuno_sepsis <-
    lmtp_tmle(dat_immuno,
              trt = A,
              outcome = Y,
              baseline = W_nt,
              cens = C,
              outcome_type = "survival",
              folds = 5,
              .SL_folds = 5,
              shift = static_binary_on,
              learners_outcome = lrnrs,
              learners_trt = lrnrs)



immuno_no_sepsis <-
    lmtp_tmle(dat_immuno,
              trt = A,
              outcome = Y,
              baseline = W_nt,
              cens = C,
              outcome_type = "survival",
              folds = 5,
              .SL_folds = 5,
              shift = static_binary_off,
              learners_outcome = lrnrs,
              learners_trt = lrnrs)

save(transplant_sepsis,
     transplant_no_sepsis,
     no_transplant_sepsis,
     no_transplant_no_sepsis,
     dx_transplant_sepsis,
     dx_transplant_no_sepsis,
     immuno_sepsis,
     immuno_no_sepsis,
     file = here::here("output/results.rdata")
)