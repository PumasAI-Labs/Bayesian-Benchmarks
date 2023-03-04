rm(list = ls())
cat("\014")

library(trelliscopejs)
library(mrgsolve)
library(tidybayes)
library(cmdstanr)
library(tidyverse)

set_cmdstan_path("cmdstan")

model_simulate <- 
  cmdstan_model(
    "Model2/iv_1cmt_linear/Stan/Torsten/Simulate/iv_1cmt_ppa.stan") 

TVCL <- 4
TVVC <- 70
# TVKA <- 1

omega_cl <- 0.3
omega_vc <- 0.3
# omega_ka <- 0.3

R <- diag(rep(1, times = 2))

sigma_p <- 0.2
sigma_a <- 0
cor_p_a <- 0

n_subjects_per_dose <- 30

dosing_data <- expand.ev(ID = 1:n_subjects_per_dose, addl = 0, ii = 0, 
                         cmt = 2, amt = c(10, 30, 60, 120), tinf = 1,
                         evid = 1) %>%
  as_tibble() %>% 
  rename_all(toupper) %>%
  mutate(SS = 0) %>% 
  select(ID, TIME, everything())


times_to_simulate <- c(0.25, 0.5, 0.75, 1, 2, 4, 8, 12, 
                       16, 20, 24, 36, 48, 60, 72)

nonmem_data_simulate <- dosing_data %>% 
  group_by(ID) %>% 
  slice(rep(1, times = length(times_to_simulate))) %>% 
  mutate(AMT = 0,
         CMT = 2,
         EVID = 0,
         TIME = times_to_simulate,
         RATE = 0,
         TINF = 0) %>% 
  ungroup() %>%
  bind_rows(dosing_data) %>% 
  arrange(ID, TIME, AMT) %>% 
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
                  time = nonmem_data_simulate$TIME,
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
                  omega_cl = omega_cl,
                  omega_vc = omega_vc,
                  R = R,
                  sigma_p = sigma_p,
                  sigma_a = sigma_a,
                  cor_p_a = cor_p_a)

simulated_data <- model_simulate$sample(data = stan_data,
                                        fixed_param = TRUE,
                                        seed = 1234,
                                        iter_warmup = 0,
                                        iter_sampling = 1,
                                        chains = 1,
                                        parallel_chains = 1)

params_ind <- simulated_data$draws(c("CL", "VC")) %>% 
  spread_draws(CL[i], VC[i]) %>% 
  inner_join(dosing_data %>% 
               mutate(i = 1:n()),
             by = "i") %>% 
  ungroup() %>%
  select(ID, CL, VC)

data <- simulated_data$draws(c("dv")) %>% 
  spread_draws(dv[i]) %>% 
  ungroup() %>% 
  inner_join(nonmem_data_simulate %>% 
               mutate(i = 1:n()),
             by = "i") %>%
  select(ID, AMT, RATE, II, ADDL, CMT, EVID, SS, TIME, DV = "dv") %>% 
  mutate(LLOQ = 0.0001, 
         BLOQ = case_when(EVID == 1 ~ NA_real_,
                          DV <= LLOQ ~ 1,
                          DV > LLOQ ~ 0,
                          TRUE ~ NA_real_),
         DV = if_else(EVID == 1 | BLOQ == 1, NA_real_, DV),
         MDV = if_else(is.na(DV), 1, 0),
         CMT = 1) %>% 
  relocate(DV, .after = last_col()) %>% 
  relocate(TIME, .before = DV)

ggplot(data %>% 
         group_by(ID) %>% 
         mutate(Dose = factor(max(AMT)))) +
  geom_point(aes(x = TIME, y = DV, group = ID, color = Dose)) +
  geom_line(aes(x = TIME, y = DV, group = ID, color = Dose)) +
  theme_bw(18) +
  scale_y_continuous(name = latex2exp::TeX("Drug Conc. $(\\mu g/mL)"),
                     trans = "log10") + 
  scale_x_continuous(name = "Time (h)",
                     breaks = seq(0, 72, by = 12),
                     labels = seq(0, 72, by = 12)) +
  geom_point(data = data %>% 
               filter(BLOQ == 1),
             mapping = aes(x = TIME, y = LLOQ), color = "limegreen") #+ 
# facet_trelliscope(~ID, nrow = 4, ncol = 5)

data %>% 
  write_csv(file.path("Model2", "iv_1cmt_linear", "Data", "single_dose.csv"),
            na = ".")

params_ind %>% 
  write_csv("Model2/iv_1cmt_linear/Data/single_dose_params_ind.csv")

################ Now Simulate Multiple Doses ####################

dosing_data <- expand.ev(ID = 1:n_subjects_per_dose, addl = 6, ii = 24, 
                         cmt = 2, amt = c(10, 30, 60, 120), tinf = 1,
                         evid = 1) %>%
  as_tibble() %>% 
  rename_all(toupper) %>%
  mutate(SS = 0) %>% 
  select(ID, TIME, everything())

times_to_simulate <- c(24, 48, 72, 72.25, 72.5, 72.75, 73, 73.5, 74, 74.5, 75,
                       76, 78, 80, 84, 88, 92, 96, 120, 144, 144.25, 144.5, 
                       144.75, 145, 145.5, 146, 146.5, 147, 148, 150, 152, 156, 
                       160, 164, 168, 180, 192, 204, 216)

nonmem_data_simulate <- dosing_data %>% 
  group_by(ID) %>% 
  slice(rep(1, times = length(times_to_simulate))) %>% 
  mutate(AMT = 0,
         ADDL = 0,
         II = 0,
         CMT = 2,
         EVID = 0,
         TIME = times_to_simulate,
         RATE = 0,
         TINF = 0) %>% 
  ungroup() %>%
  bind_rows(dosing_data) %>% 
  arrange(ID, TIME, AMT) %>% 
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
                  time = nonmem_data_simulate$TIME,
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
                  omega_cl = omega_cl,
                  omega_vc = omega_vc,
                  R = R,
                  sigma_p = sigma_p,
                  sigma_a = sigma_a,
                  cor_p_a = cor_p_a)

simulated_data <- model_simulate$sample(data = stan_data,
                                        fixed_param = TRUE,
                                        seed = 1234,
                                        iter_warmup = 0,
                                        iter_sampling = 1,
                                        chains = 1,
                                        parallel_chains = 1)

params_ind <- simulated_data$draws(c("CL", "VC")) %>% 
  spread_draws(CL[i], VC[i]) %>% 
  inner_join(dosing_data %>% 
               mutate(i = 1:n()),
             by = "i") %>% 
  ungroup() %>%
  select(ID, CL, VC)

data <- simulated_data$draws(c("dv")) %>% 
  spread_draws(dv[i]) %>% 
  ungroup() %>% 
  inner_join(nonmem_data_simulate %>% 
               mutate(i = 1:n()),
             by = "i") %>%
  select(ID, AMT, RATE, II, ADDL, CMT, EVID, SS, TIME, DV = "dv") %>% 
  mutate(LLOQ = 0.0001, 
         BLOQ = case_when(EVID == 1 ~ NA_real_,
                          DV <= LLOQ ~ 1,
                          DV > LLOQ ~ 0,
                          TRUE ~ NA_real_),
         DV = if_else(EVID == 1 | BLOQ == 1, NA_real_, DV),
         MDV = if_else(is.na(DV), 1, 0),
         CMT = 1) %>% 
  relocate(DV, .after = last_col()) %>% 
  relocate(TIME, .before = DV)

ggplot(data %>% 
         filter(!is.na(DV))) +
  geom_point(mapping = aes(x = TIME, y = DV, group = ID)) +
  geom_line(mapping = aes(x = TIME, y = DV, group = ID)) +
  theme_bw(18) +
  scale_y_continuous(name = latex2exp::TeX("Drug Conc. $(\\mu g/mL)"),
                     trans = "log10") + 
  scale_x_continuous(name = "Time (h)",
                     breaks = seq(0, 216, by = 48),
                     labels = seq(0, 216, by = 48),
                     limits = c(0, 216)) + 
  geom_vline(data = data %>% 
               filter(EVID == 1) %>% 
               mrgsolve::realize_addl(), 
             mapping = aes(xintercept = TIME, group = ID), linetype = "dashed",
             color = "magenta") +
  geom_point(data = data %>% 
               filter(BLOQ == 1),
             mapping = aes(x = TIME, y = LLOQ), color = "limegreen") #+
# facet_trelliscope(~ID, nrow = 3, ncol = 4, scales = "free_y")

data %>% 
  write_csv(file.path("Model2", "iv_1cmt_linear", "Data", "multiple_dose.csv"),
            na = ".")

params_ind %>% 
  write_csv("Model2/iv_1cmt_linear/Data/multiple_dose_params_ind.csv")


