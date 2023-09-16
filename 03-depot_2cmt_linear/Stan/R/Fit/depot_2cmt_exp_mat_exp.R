rm(list = ls())
cat("\014")

library(cmdstanr)
library(tidyverse)

set_cmdstan_path("cmdstan")

read_data_and_prepare_stan_data <- function(n_doses = c("single", "multiple")){
  
  n_doses <- match.arg(n_doses)
  
  nonmem_data <- read_csv(str_c("03-depot_2cmt_linear/data/", 
                                n_doses, "_dose.csv"),
                          na = ".") %>% 
    rename_all(tolower) %>% 
    rename(ID = "id",
           DV = "dv") %>% 
    mutate(DV = if_else(is.na(DV), 5555555, DV),    # This value can be anything except NA. It'll be indexed away 
           bloq = if_else(is.na(bloq), -999, bloq)) # This value can be anything except NA. It'll be indexed away 
  
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
                    location_tvq = 4,
                    location_tvvp = 40,
                    location_tvka = 1,
                    scale_tvcl = 1,
                    scale_tvvc = 1,
                    scale_tvq = 1,
                    scale_tvvp = 1,
                    scale_tvka = 1,
                    scale_omega_cl = 0.4,
                    scale_omega_vc = 0.4,
                    scale_omega_q = 0.4,
                    scale_omega_vp = 0.4,
                    scale_omega_ka = 0.4,
                    lkj_df_omega = 2,
                    scale_sigma = 0.5)
  
  return(stan_data)
  
}

sample_and_save_all <- function(n_doses = c("single", "multiple"), run_number){
  
  n_doses <- match.arg(n_doses)
  
  print(str_c(n_doses, " dose, run #", run_number))
  
  init_files <- str_c("03-depot_2cmt_linear/data/inits/inits_", run_number, "_", 
                      1:4, ".json")
  
  stan_data <- read_data_and_prepare_stan_data(n_doses = n_doses)
  
  model <- cmdstan_model(
    "03-depot_2cmt_linear/Stan/Torsten/Fit/depot_2cmt_exp_mat_exp.stan",
    cpp_options = list(stan_threads = TRUE))
  
  fit <- model$sample(data = stan_data,
                      seed = 112356,
                      chains = 4,
                      parallel_chains = 4,
                      threads_per_chain = parallel::detectCores()/4,
                      iter_warmup = 500,
                      iter_sampling = 1000,
                      adapt_delta = 0.8,
                      refresh = 500,
                      max_treedepth = 10,
                      init = init_files)
  
  fit$save_object(str_c("03-depot_2cmt_linear/Stan/Torsten/Fits/", 
                        n_doses, "_dose_mat_exp_", run_number, ".rds"))
  
  parameters_to_summarize <- c(str_subset(fit$metadata()$stan_variables, "TV"),
                               str_subset(fit$metadata()$stan_variables, "omega"),
                               "sigma")
  
  fit$draws(parameters_to_summarize, format = "draws_df") %>% 
    as_tibble() %>% 
    write_csv(str_c("03-depot_2cmt_linear/Stan/Torsten/Fits/", 
                    n_doses, "_dose_draws_df_mat_exp_", run_number, ".csv"))
  
}

expand_grid(n_doses = c("single", "multiple"), run_number = 1:5) %>% 
  pwalk(.f = sample_and_save_all)