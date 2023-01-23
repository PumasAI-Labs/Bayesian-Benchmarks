rm(list = ls())
cat("\014")

library(mrgsolve)
library(tidybayes)
library(cmdstanr)
library(tidyverse)

set_cmdstan_path("~/Torsten/cmdstan")

model_simulate <- 
  cmdstan_model("Model1/Stan/Torsten/Simulate/depot_1cmt_ppa.stan") 

TVCL <- 4
TVVC <- 70
TVKA <- 1

omega_cl <- 0.3
omega_vc <- 0.3
omega_ka <- 0.3

R <- diag(rep(1, times = 3))

sigma_p <- 0.2
sigma_a <- 0
cor_p_a <- 0

n_trials <- 500
n_subjects_per_dose <- 30

dosing_data <- expand.ev(ID = 1:n_subjects_per_dose, addl = 0, ii = 0, 
                         cmt = 1, amt = c(10, 30, 60, 120), tinf = 0,
                         evid = 1) %>%
  as_tibble() %>% 
  rename_all(toupper) %>%
  mutate(SS = 0) %>% 
  select(ID, TIME, everything())

nonmem_data_simulate <- dosing_data %>% 
  group_by(ID) %>% 
  slice(rep(1, times = 4)) %>% 
  mutate(AMT = 0,
         CMT = 2,
         EVID = 0) %>% 
  arrange(TIME, .by_group = TRUE) %>% 
  ungroup() %>%
  bind_rows(dosing_data) %>% 
  arrange(ID, desc(AMT)) %>% 
  select(-TINF)

n_subjects <- nonmem_data_simulate %>%  # number of individuals to simulate
  distinct(ID) %>% 
  count() %>% 
  deframe()

n_total <- nrow(nonmem_data_simulate) # total number of time points at which to predict

subj_start <- nonmem_data_simulate %>% 
  mutate(row_num = 1:n()) %>% 
  group_by(ID) %>% 
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

simulate_data_and_write_csv <- function(trial_number, nonmem_data_simulate){
  
  nonmem_data_simulate <- nonmem_data_simulate %>% 
    group_by(ID) %>% 
    mutate(TIME = c(0, runif(4, 0, 24))) %>% 
    arrange(TIME, .by_group = TRUE) %>% 
    ungroup()
  
  stan_data$time <- nonmem_data_simulate$TIME
  
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
    mutate(TRIAL = trial_number) %>% 
    select(TRIAL, ID, CL, VC, KA)
  
  data <- simulated_data$draws(c("dv")) %>% 
    spread_draws(dv[i]) %>% 
    ungroup() %>% 
    inner_join(nonmem_data_simulate %>% 
                 mutate(i = 1:n()),
               by = "i") %>%
    mutate(TRIAL = trial_number) %>% 
    select(TRIAL, ID, AMT, RATE, II, ADDL, CMT, EVID, SS, TIME, DV = "dv")
  
  data %>% 
    write_csv(file.path("Model1", "Data", 
                        str_c(str_pad(trial_number, 3, side = "left", 
                                      pad = "0"), 
                              ".csv")))
  
  return(params_ind)
  
}


params_ind <- map_dfr(1:n_trials, simulate_data_and_write_csv,
                      nonmem_data_simulate = nonmem_data_simulate)

params_ind %>% 
  write_csv("Model1/Data/params_ind.csv")

