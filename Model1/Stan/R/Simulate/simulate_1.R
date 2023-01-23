rm(list = ls())
cat("\014")

library(mrgsolve)
library(tidybayes)
library(cmdstanr)
library(tidyverse)

set_cmdstan_path("~/Torsten/cmdstan")

model_simulate <- 
  cmdstan_model("Model1/Stan/Torsten/Simulate/depot_1cmt_ppa.stan") 

create_dosing_data <- function(trial_number, n_subjects_per_dose = 30, 
                               dose_amounts = c(10, 30, 60, 120), 
                               addl = 0, ii = 0, cmt = 1, tinf = 0){
  
  dosing_data <- expand.ev(ID = 1:n_subjects_per_dose, addl = addl, ii = ii, 
                           cmt = cmt, amt = dose_amounts, tinf = tinf,
                           evid = 1) %>%
    realize_addl() %>% 
    as_tibble() %>% 
    rename_all(toupper) %>%
    mutate(SS = 0,
           TRIAL = trial_number) %>% 
    select(TRIAL, ID, TIME, everything())
  
  return(dosing_data)
  
}

TVCL <- 4
TVVC <- 70
TVKA <- 1

omega_cl <- 0.3
omega_vc <- 0.3
omega_ka <- 0.3

R <- diag(rep(1, times = 3))

sigma_p <- 0.3
sigma_a <- 0
cor_p_a <- 0

n_trials <- 500

# dosing_data <- create_dosing_data(1) 
dosing_data <- map_dfr(1:n_trials, create_dosing_data)
  
nonmem_data_simulate <- dosing_data %>% 
  group_by(TRIAL, ID) %>% 
  slice(rep(1, times = 4)) %>% 
  mutate(AMT = 0,
         CMT = 2,
         EVID = 0,
         TIME = runif(n(), 0, 24)) %>% 
  arrange(TIME, .by_group = TRUE) %>% 
  ungroup() %>%
  bind_rows(dosing_data) %>% 
  arrange(TRIAL, ID, TIME) %>% 
  select(-TINF)
  
n_subjects <- nonmem_data_simulate %>%  # number of individuals to simulate
  distinct(TRIAL, ID) %>% 
  count() %>% 
  deframe()

n_total <- nrow(nonmem_data_simulate) # total number of time points at which to predict

subj_start <- nonmem_data_simulate %>% 
  mutate(row_num = 1:n()) %>% 
  group_by(TRIAL, ID) %>% 
  slice_head(n = 1) %>%
  ungroup() %>% 
  select(row_num) %>% 
  deframe()

subj_end <- c(subj_start[-1] - 1, n_total) 

stan_data <- list(n_subjects = n_subjects,
                  n_total = n_total,
                  amt = nonmem_data_simulate$AMT,
                  cmt = nonmem_data_simulate$CMT,
                  evid = nonmem_data_simulate$EVID,
                  rate = nonmem_data_simulate$RATE,
                  ii = nonmem_data_simulate$II,
                  addl = nonmem_data_simulate$ADDL,
                  ss = nonmem_data_simulate$SS,
                  time = nonmem_data_simulate$TIME,
                  subj_start = subj_start,
                  subj_end = subj_end,
                  TVCL = TVCL,
                  TVVC = TVVC,
                  TVKA = TVKA,
                  omega_cl = omega_cl,
                  omega_vc = omega_vc,
                  omega_ka = omega_ka,
                  R = R,
                  sigma_p = sigma_p,
                  sigma_a = sigma_a,
                  cor_p_a = cor_p_a)

simulated_data <- model_simulate$sample(data = stan_data,
                                        fixed_param = TRUE,
                                        seed = 112358,
                                        iter_warmup = 0,
                                        iter_sampling = 1,
                                        chains = 1,
                                        parallel_chains = 1)

params_ind <- simulated_data$draws(c("CL", "VC", "KA")) %>% 
  spread_draws(CL[i], VC[i], KA[i]) %>% 
  inner_join(dosing_data %>% 
               mutate(i = 1:n()),
             by = "i") %>% 
  ungroup() %>% 
  select(TRIAL, ID, CL, VC, KA)

data <- simulated_data$draws(c("dv"), format = "draws_df") %>% 
  spread_draws(dv[i]) %>% 
  ungroup() %>% 
  inner_join(nonmem_data_simulate %>% 
               mutate(i = 1:n()),
             by = "i") %>%
  select(TRIAL, ID, AMT, RATE, II, ADDL, CMT, EVID, SS, TIME, DV = "dv")
  
  
data %>% 
  group_split(TRIAL) %>% 
  walk(~ write_csv(.x, file.path("Model1", "Data", 
                                 str_c(str_pad(.x$TRIAL[1], 3, side = "left",
                                               pad = "0"), 
                                       ".csv"))))
