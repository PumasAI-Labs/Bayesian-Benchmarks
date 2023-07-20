rm(list = ls())
cat("\014")

library(mrgsolve)
library(tidybayes)
library(cmdstanr)
library(tidyverse)

set_cmdstan_path("cmdstan")

model_simulate <- 
  cmdstan_model("05-friberg/Stan/Torsten/Simulate/depot_2cmt_exp_friberg.stan") 

TVCL <- 4        # L/h
TVVC <- 70       # L
TVQ <- 4         # L/h
TVVP <- 50       # L
TVKA <- 1        # h^-1
TVMTT <- 125     # h
TVCIRC0 <- 5     # x/L
TVGAMMA <- 0.17  # 
TVALPHA <- 3e-4  # 

omega_cl <- 0.3
omega_vc <- 0.3
omega_q <- 0.3
omega_vp <- 0.3
omega_ka <- 0.3
omega_mtt <- 0.3
omega_circ_0 <- 0.3
omega_gamma <- 0.3
omega_alpha <- 0.3

# omega_cl <- 0
# omega_vc <- 0
# omega_q <- 0
# omega_vp <- 0
# omega_ka <- 0
# omega_mtt <- 0
# omega_circ_0 <- 0
# omega_gamma <- 0
# omega_alpha <- 0

R <- diag(rep(1, times = 9))

sigma <- 0.2

sigma_pd <- 0.2

n_subjects_per_dose <- 3

dosing_data <- expand.ev(ID = 1:n_subjects_per_dose, addl = 6, ii = 24,
                         cmt = 1, amt = c(10, 20, 40, 80)*1000, tinf = 0,
                         evid = 1) %>%
  as_tibble() %>%
  rename_all(toupper) %>%
  mutate(SS = 0) %>%
  select(ID, TIME, everything())

times_to_simulate_pk <- c(24, 48, 72, 72.25, 72.5, 72.75, 73, 73.5, 74, 74.5, 
                          75, 76, 78, 80, 84, 88, 92, 96, 120, 144, 144.25, 
                          144.5, 144.75, 145, 145.5, 146, 146.5, 147, 148, 150, 
                          152, 156, 160, 164, 168, 180, 192, 204, 216)

times_to_simulate_pd <- seq(0, 672, by = 24)

nonmem_data_pk <- dosing_data %>% 
  group_by(ID) %>% 
  slice(rep(1, times = length(times_to_simulate_pk))) %>% 
  mutate(AMT = 0,
         CMT = 2,
         EVID = 0,
         TIME = times_to_simulate_pk,
         RATE = 0,
         TINF = 0,
         ADDL = 0,
         II = 0) %>% 
  ungroup() 

nonmem_data_pd <- dosing_data %>% 
  group_by(ID) %>% 
  slice(rep(1, times = length(times_to_simulate_pd))) %>% 
  mutate(AMT = 0,
         CMT = 3,
         EVID = 0,
         TIME = times_to_simulate_pd,
         RATE = 0,
         TINF = 0,
         ADDL = 0,
         II = 0) %>% 
  ungroup()

nonmem_data_simulate <- nonmem_data_pk %>%
  bind_rows(nonmem_data_pd, dosing_data) %>% 
  arrange(ID, TIME, AMT, CMT) %>% 
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
                  TVQ = TVQ,
                  TVVP = TVVP,
                  TVKA = TVKA,
                  TVMTT = TVMTT,
                  TVCIRC0 = TVCIRC0,
                  TVALPHA = TVALPHA,
                  TVGAMMA = TVGAMMA,
                  omega_cl = omega_cl,
                  omega_vc = omega_vc,
                  omega_q = omega_q,
                  omega_vp = omega_vp,
                  omega_ka = omega_ka,
                  omega_mtt = omega_mtt,
                  omega_circ_0 = omega_circ_0,
                  omega_gamma = omega_gamma,
                  omega_alpha = omega_alpha,
                  R = R,
                  sigma = sigma,
                  sigma_pd = sigma_pd)

simulated_data <- model_simulate$sample(data = stan_data,
                                        fixed_param = TRUE,
                                        seed = 123,
                                        iter_warmup = 0,
                                        iter_sampling = 1,
                                        chains = 1,
                                        parallel_chains = 1)

params_ind <- simulated_data$draws(c("CL", "VC", "Q", "VP", "KA",
                                     "MTT", "CIRC0", "ALPHA", "GAMMA")) %>% 
  spread_draws(CL[i], VC[i], Q[i], VP[i], KA[i],
               MTT[i], CIRC0[i], ALPHA[i], GAMMA[i]) %>% 
  inner_join(dosing_data %>% 
               mutate(i = 1:n()),
             by = "i") %>% 
  ungroup() %>%
  select(ID, CL, VC, Q, VP, KA, MTT, CIRC0, ALPHA, GAMMA)

data <- simulated_data$draws(c("cp", "dv")) %>% 
  spread_draws(cp[i], dv[i]) %>% 
  ungroup() %>% 
  inner_join(nonmem_data_simulate %>% 
               mutate(i = 1:n()),
             by = "i") %>%
  select(ID, AMT, RATE, II, ADDL, CMT, EVID, SS, TIME, CP = "cp", DV = "dv") %>% 
  mutate(LLOQ = if_else(CMT == 2, 0, 0), 
         BLOQ = case_when(EVID == 1 ~ NA_real_,
                          DV <= LLOQ ~ 1,
                          DV > LLOQ ~ 0,
                          TRUE ~ NA_real_),
         DV = if_else(EVID == 1 | BLOQ == 1, NA_real_, DV),
         DVLN = log(DV),
         MDV = if_else(is.na(DV), 1, 0)) %>%
  relocate(CP, .after = last_col()) %>% 
  relocate(DV, .after = last_col()) %>% 
  relocate(DVLN, .after = last_col()) %>% 
  relocate(TIME, .before = CP)

ggplot(data %>% 
         group_by(ID) %>% 
         mutate(Dose = factor(max(AMT))) %>% 
         ungroup() %>% 
         filter(EVID == 0) %>% 
         inner_join(params_ind, by = "ID")) +
  geom_point(aes(x = TIME, y = DV, group = ID, color = Dose)) +
  geom_line(aes(x = TIME, y = DV, group = ID, color = Dose)) +
  # geom_line(aes(x = TIME, y = CP, group = ID, color = Dose)) +
  theme_bw(18) +
  # scale_y_continuous(name = latex2exp::TeX("Drug Conc. $(\\mu g/mL)"),
  #                    trans = "log10") + 
  scale_y_continuous(name = latex2exp::TeX("Dependent Variable"),
                     trans = "log10") + 
  scale_x_continuous(name = "Time (h)",
                     breaks = seq(0, max(data$TIME), by = 24),
                     labels = seq(0, max(data$TIME), by = 24),
                     limits = c(0, NA)) +
  # geom_point(data = data %>% 
  #              filter(BLOQ == 1),
  #            mapping = aes(x = TIME, y = LLOQ), color = "limegreen") +
  facet_wrap(~CMT, labeller = label_both, scales = "free_y", 
             ncol = 1) #+ 
# facet_trelliscope(~ID, nrow = 3, ncol = 4)

data %>%
  write_csv(file.path("05-friberg", "data", "multiple_dose.csv"),
            na = ".")

params_ind %>%
  write_csv("05-friberg/data/multiple_dose_params_ind.csv")

