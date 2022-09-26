rm(list = ls())
cat("\014")

library(cmdstanr)
library(tidyverse)

set_cmdstan_path("~/Torsten/cmdstan")

data <- read_rds("R/Data/iv_2cmt_macro_ir1_cl_vc_q_vp_kin_kout_ic50_ppa_lognormal_001_bloq.rds")

nonmem_data <- data$data_nonmem %>% 
  rename_all(tolower) %>% 
  rename(ID = "id",
         DV = "dv") %>% 
  mutate(amt = if_else(is.na(amt), 0, amt), 
         rate = if_else(evid == 0, 0, rate),
         ii = 0,
         addl = 0,
         ss = 0,
         DV = if_else(bloq %in% c(1, NA), 5555555, DV), # This value can be anything > 0. It'll be indexed away 
         bloq = if_else(is.na(bloq), -999, bloq), 
         cmt = cmt + 1) %>%                             # cmt = 2 for PK, cmt = 4 for PD
  filter(evid %in% c(0, 1)) 

tmp <- ggplot(nonmem_data %>%
                mutate(ID = factor(ID)) %>%
                group_by(ID) %>%
                mutate(Dose = factor(max(amt, na.rm = TRUE))) %>%
                ungroup() %>%
                filter(mdv == 0)) +
    geom_line(mapping = aes(x = time, y = DV, group = ID, color = Dose)) +
    geom_point(mapping = aes(x = time, y = DV, group = ID, color = Dose)) +
    scale_color_discrete(name = "Dose (mg)") +
    # scale_y_log10(name = "Drug Conc. (ug/mL)",
    scale_y_log10(name = "DV",
                  limits = c(NA, NA)) +
    scale_x_continuous(name = "Time (w)",
                       breaks = seq(0, 168, by = 28),
                       labels = seq(0, 168/7, by = 28/7)) +
    theme_bw(18) +
    theme(axis.text = element_text(size = 14, face = "bold"),
          axis.title = element_text(size = 18, face = "bold"),
          axis.line = element_line(size = 2),
          legend.position = "bottom") +
    geom_hline(aes(yintercept = lloq), size = 1.2, linetype = 2) +
    geom_vline(data = nonmem_data %>%
                 mutate(ID = factor(ID)) %>%
                 filter(evid == 1),
               mapping = aes(xintercept = time, group = ID),
               linetype = 2, color = "red", alpha = 0.5) +
    ggforce::facet_wrap_paginate(~ID + cmt, labeller = label_both,
                                 nrow = 3, ncol = 4,
                                 page = 1, scales = "free")


for(i in 1:ggforce::n_pages(tmp)){
  print(ggplot(nonmem_data %>%
                 mutate(ID = factor(ID)) %>%
                 group_by(ID) %>%
                 mutate(Dose = factor(max(amt, na.rm = TRUE))) %>%
                 ungroup() %>%
                 filter(mdv == 0)) +
          geom_line(mapping = aes(x = time, y = DV, group = ID, color = Dose)) +
          geom_point(mapping = aes(x = time, y = DV, group = ID, color = Dose)) +
          scale_color_discrete(name = "Dose (mg)") +
          # scale_y_log10(name = "Drug Conc. (ug/mL)",
          scale_y_log10(name = "DV",
                        limits = c(NA, NA)) +
          scale_x_continuous(name = "Time (w)",
                             breaks = seq(0, 168, by = 28),
                             labels = seq(0, 168/7, by = 28/7)) +
          theme_bw(18) +
          theme(axis.text = element_text(size = 14, face = "bold"),
                axis.title = element_text(size = 18, face = "bold"),
                axis.line = element_line(size = 2),
                legend.position = "bottom") +
          geom_hline(aes(yintercept = lloq), size = 1.2, linetype = 2) +
          geom_vline(data = nonmem_data %>%
                       mutate(ID = factor(ID)) %>%
                       filter(evid == 1),
                     mapping = aes(xintercept = time, group = ID),
                     linetype = 2, color = "red", alpha = 0.5) +
          ggforce::facet_wrap_paginate(~ID + cmt, labeller = label_both,
                                       nrow = 3, ncol = 4,
                                       page = i, scales = "free"))
  
}



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
                  scale_tvcl = 2,
                  scale_tvvc = 10,
                  scale_tvq = 2,
                  scale_tvvp = 20,
                  scale_omega_cl = 0.4,
                  scale_omega_vc = 0.4,
                  scale_omega_q = 0.4,
                  scale_omega_vp = 0.4,
                  lkj_df_omega = 2,
                  scale_sigma_p = 0.5,
                  scale_sigma_a = 1,
                  lkj_df_sigma = 2,
                  scale_tvkin = 100,
                  scale_tvkout = 2,
                  scale_tvic50 = 5,
                  scale_omega_kin = 0.4,
                  scale_omega_kout = 0.4,
                  scale_omega_ic50 = 0.4,
                  lkj_df_omega_pd = 0.2,
                  scale_sigma_pd = 0.4)

model <- cmdstan_model(
  "Torsten/Fit/iv_2cmt_macro_ir1_cl_vc_q_vp_kin_kout_ic50_ppa_lognormal_bloq_lower.stan",
  cpp_options = list(stan_threads = TRUE))

fit <- model$sample(data = stan_data,
                    chains = 4,
                    parallel_chains = 4,
                    threads_per_chain = 12,
                    iter_warmup = 100,
                    iter_sampling = 500,
                    adapt_delta = 0.98,
                    refresh = 5,
                    max_treedepth = 15,
                    # output_dir = "Torsten/Fits",
                    # output_basename = "ir1",
                    # save_warmup = TRUE,
                    init = function() list(TVCL = rlnorm(1, log(0.3), 0.3),
                                           TVVC = rlnorm(1, log(5), 0.3),
                                           TVQ = rlnorm(1, log(1), 0.3),
                                           TVVP = rlnorm(1, log(10), 0.3),
                                           omega = rlnorm(4, log(0.3), 0.3),
                                           sigma = rlnorm(2, log(0.5), 0.3),
                                           TVKIN = rlnorm(1, log(110), 0.3),
                                           TVKOUT = rlnorm(1, log(1), 0.3),
                                           TVIC50 = rlnorm(1, log(6), 0.3),
                                           omega_pd = rlnorm(3, log(0.3), 0.3),
                                           sigma_pd = rlnorm(1, log(0.1), 0.3)))

fit$save_object(
  "Torsten/Fits/iv_2cmt_macro_ir1_cl_vc_q_vp_kin_kout_ic50_ppa_lognormal_001_bloq_lower_general_500.rds")


model_coupled <- cmdstan_model(
  "Torsten/Fit/iv_2cmt_macro_ir1_cl_vc_q_vp_kin_kout_ic50_ppa_lognormal_bloq_lower_coupled.stan",
  cpp_options = list(stan_threads = TRUE))

fit_coupled <- model_coupled$sample(data = stan_data,
                    chains = 4,
                    parallel_chains = 4,
                    threads_per_chain = 12,
                    iter_warmup = 100,
                    iter_sampling = 500,
                    adapt_delta = 0.98,
                    refresh = 5,
                    max_treedepth = 15,
                    # output_dir = "Torsten/Fits",
                    # output_basename = "ir1",
                    # save_warmup = TRUE,
                    init = function() list(TVCL = rlnorm(1, log(0.3), 0.3),
                                           TVVC = rlnorm(1, log(5), 0.3),
                                           TVQ = rlnorm(1, log(1), 0.3),
                                           TVVP = rlnorm(1, log(10), 0.3),
                                           omega = rlnorm(4, log(0.3), 0.3),
                                           sigma = rlnorm(2, log(0.5), 0.3),
                                           TVKIN = rlnorm(1, log(110), 0.3),
                                           TVKOUT = rlnorm(1, log(1), 0.3),
                                           TVIC50 = rlnorm(1, log(6), 0.3),
                                           omega_pd = rlnorm(3, log(0.3), 0.3),
                                           sigma_pd = rlnorm(1, log(0.1), 0.3)))

fit_coupled$save_object(
  "Torsten/Fits/iv_2cmt_macro_ir1_cl_vc_q_vp_kin_kout_ic50_ppa_lognormal_001_bloq_lower_coupled_500.rds")
