rm(list = ls())
cat("\014")

library(patchwork)
library(cmdstanr)
library(tidyverse)

set_cmdstan_path("cmdstan")

nonmem_data <- read_csv("05-friberg/data/multiple_dose.csv",
                        na = ".") %>% 
  rename_all(tolower) %>% 
  rename(ID = "id",
         DV = "dv") %>% 
  mutate(DV = if_else(is.na(DV), 5555555, DV),    # This value can be anything except NA. It'll be indexed away 
         bloq = if_else(is.na(bloq), -999, bloq)) # This value can be anything except NA. It'll be indexed away 


(p_pk <- ggplot(nonmem_data %>%
                  mutate(ID = factor(ID)) %>%
                  group_by(ID) %>%
                  mutate(Dose = factor(max(amt, na.rm = TRUE))) %>%
                  ungroup() %>%
                  filter(mdv == 0, cmt == 2)) +
    geom_line(mapping = aes(x = time, y = DV, group = ID, color = Dose)) +
    geom_point(mapping = aes(x = time, y = DV, group = ID, color = Dose)) +
    scale_color_discrete(name = "Dose (ug)") +
    scale_y_continuous(name = "Drug Conc. (ng/mL)",
                       limits = c(NA, NA),
                       trans = "log10") +
    scale_x_continuous(name = "Time (d)",
                       breaks = seq(0, 216, by = 24),
                       labels = seq(0, 216/24, by = 24/24),
                       limits = c(0, NA)) +
    theme_bw(18) +
    theme(axis.text = element_text(size = 14, face = "bold"),
          axis.title = element_text(size = 18, face = "bold"),
          axis.line = element_line(size = 2),
          legend.position = "bottom"))

(p_pd <- ggplot(nonmem_data %>%
                  mutate(ID = factor(ID)) %>%
                  group_by(ID) %>%
                  mutate(Dose = factor(max(amt, na.rm = TRUE))) %>%
                  ungroup() %>%
                  filter(mdv == 0, cmt == 3)) +
    geom_line(mapping = aes(x = time, y = DV, group = ID, color = Dose)) +
    geom_point(mapping = aes(x = time, y = DV, group = ID, color = Dose)) +
    scale_color_discrete(name = "Dose (ug)") +
    scale_y_continuous(name = latex2exp::TeX("$Neutrophils\\;(10^9/L)"),
                       limits = c(NA, NA),
                       trans = "identity") +
    scale_x_continuous(name = "Time (d)",
                       breaks = seq(0, 672, by = 24),
                       labels = seq(0, 672/24, by = 24/24),
                       limits = c(0, NA)) +
    theme_bw(18) +
    theme(axis.text = element_text(size = 14, face = "bold"),
          axis.title = element_text(size = 18, face = "bold"),
          axis.line = element_line(size = 2),
          legend.position = "bottom"))

p_pk +
  p_pd +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

nonmem_data %>%
  filter(evid == 0) %>% 
  group_by(ID, cmt) %>% 
  summarize(lloq = unique(lloq),
            n_obs = n(),
            n_bloq = sum(bloq == 1)) %>% 
  filter(n_bloq > 0)

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
                  location_tvmtt = 125,
                  location_tvcirc0 = 5,
                  location_tvgamma = 0.17,
                  location_tvalpha = 3e-4,
                  scale_tvcl = 1,
                  scale_tvvc = 1,
                  scale_tvq = 1,
                  scale_tvvp = 1,
                  scale_tvka = 1,
                  scale_tvmtt = 1,
                  scale_tvcirc0 = 1,
                  scale_tvgamma = 1,
                  scale_tvalpha = 1,
                  scale_omega_cl = 0.4,
                  scale_omega_vc = 0.4,
                  scale_omega_q = 0.4,
                  scale_omega_vp = 0.4,
                  scale_omega_ka = 0.4,
                  scale_omega_mtt = 0.4,
                  scale_omega_circ0 = 0.4,
                  scale_omega_gamma = 0.4,
                  scale_omega_alpha = 0.4,
                  lkj_df_omega = 2,
                  scale_sigma_p = 0.5,
                  scale_sigma_p_pd = 0.5)

sample_and_save_all <- function(n_doses = "multiple", run_number){
  
  n_doses <- match.arg(n_doses)  
  
  print(str_c(n_doses, " dose, run #", run_number))
  
  init_files <- str_c("05-friberg/data/inits/inits_", run_number, "_", 
                      1:4, ".json")
  
  model <- cmdstan_model(
    "05-friberg/Stan/Torsten/Fit/depot_2cmt_friberg_prop_coupled.stan",
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
  
  fit$save_object(str_c("05-friberg/Stan/Torsten/Fits/", 
                        n_doses, "_dose_coupled_", run_number, ".rds"))
  
  parameters_to_summarize <- c(str_subset(fit$metadata()$stan_variables, "TV"),
                               str_subset(fit$metadata()$stan_variables, "omega"),
                               str_subset(fit$metadata()$stan_variables, "sigma"))
  
  fit$draws(parameters_to_summarize, format = "draws_df") %>% 
    as_tibble() %>% 
    write_csv(str_c("05-friberg/Stan/Torsten/Fits/", 
                    n_doses, "_dose_draws_df_coupled_", run_number, ".csv"))
  
}

expand_grid(n_doses = "multiple", run_number = 1:5) %>% 
  pwalk(.f = sample_and_save_all)
