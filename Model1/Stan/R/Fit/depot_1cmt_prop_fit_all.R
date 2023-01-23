rm(list = ls())
cat("\014")

library(cmdstanr)
library(tidyverse)

set_cmdstan_path("~/Torsten/cmdstan")

read_data_and_fit <- function(trial_number, model){
  
  nonmem_data <- read_csv(file.path("Model1", "Data",
                                    str_c(str_pad(trial_number, width = 3,
                                                  side = "left", pad = "0"), 
                                          ".csv")),
                          na = ".") %>% 
    rename_all(tolower) %>% 
    rename(ID = "id",
           DV = "dv") %>% 
    mutate(DV = if_else(is.na(DV), 5555555, DV), # This value can be anything except NA. It'll be indexed away 
           mdv = if_else(evid == 0, 0, 1),
           lloq = 0,
           bloq = 0) 
  
  n_subjects <- nonmem_data %>%  # number of individuals
    distinct(ID) %>%
    count() %>%
    deframe()
  
  n_total <- nrow(nonmem_data)   # total number of records
  
  i_obs <- nonmem_data %>%
    mutate(row_num = 1:n()) %>%
    filter(evid == 0) %>%
    select(row_num) %>%
    deframe()
  
  n_obs <- length(i_obs)
  
  subj_start <- nonmem_data %>%
    mutate(row_num = 1:n()) %>%
    group_by(ID) %>%
    slice_head(n = 1) %>%
    ungroup() %>%
    select(row_num) %>%
    deframe()
  
  subj_end <- c(subj_start[-1] - 1, n_total)
  
  stan_data <- list(n_subjects = n_subjects,
                    n_total = n_total,
                    n_obs = n_obs,
                    i_obs = i_obs,
                    ID = nonmem_data$ID,
                    amt = nonmem_data$amt,
                    cmt = nonmem_data$cmt,
                    evid = nonmem_data$evid,
                    rate = nonmem_data$rate,
                    ii = nonmem_data$ii,
                    addl = nonmem_data$addl,
                    ss = nonmem_data$ss,
                    time = nonmem_data$time,
                    dv = nonmem_data$DV,
                    subj_start = subj_start,
                    subj_end = subj_end,
                    lloq = nonmem_data$lloq,
                    bloq = nonmem_data$bloq,
                    location_tvcl = 4,
                    location_tvvc = 70,
                    location_tvka = 1,
                    scale_tvcl = 1,
                    scale_tvvc = 1,
                    scale_tvka = 1,
                    scale_omega_cl = 0.4,
                    scale_omega_vc = 0.4,
                    scale_omega_ka = 0.4,
                    lkj_df_omega = 2,
                    scale_sigma_p = 0.5)
  
  fit <- model$sample(data = stan_data,
                      seed = 11235,
                      chains = 4,
                      parallel_chains = 4,
                      threads_per_chain = 2,
                      iter_warmup = 500,
                      iter_sampling = 1000,
                      adapt_delta = 0.8,
                      refresh = 0,
                      max_treedepth = 10,
                      init = function() list(TVCL = rlnorm(1, log(4), 0.3),
                                             TVVC = rlnorm(1, log(70), 0.3),
                                             TVKA = rlnorm(1, log(1), 0.3),
                                             omega = rlnorm(3, log(0.3), 0.3),
                                             sigma_p = rlnorm(1, log(0.2), 0.3)))
  

  fit$save_object(file.path("Model1", "Stan", "Torsten", "Fits",
                            str_c(str_pad(trial_number, width = 3,
                                          side = "left", pad = "0"), 
                                  ".rds")))
}

model <- cmdstan_model("Model1/Stan/Torsten/Fit/depot_1cmt_prop.stan",
                       cpp_options = list(stan_threads = TRUE))

walk(1:500, read_data_and_fit, model = model)


