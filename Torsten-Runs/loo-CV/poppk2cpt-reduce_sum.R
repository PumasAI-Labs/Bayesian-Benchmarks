library(dplyr)
library(readr)
library(tibble)
library(cmdstanr)

set_cmdstan_path("Torsten/cmdstan")

data <- rstan::read_rdump("Torsten/example-models/poppk2cpt/poppk2cpt.data.R")

nonmem_data <- tibble(
  c = NA_real_,
  time = data$time,
  amt = data$amt,
  evid = data$evid,
  rate = data$rate,
  addl = data$addl,
  ii = data$ii,
  cmt = data$cmt,
  bloq = 0,
  lloq = 0,
  ss = data$ss
) %>%
  mutate(
    i = 1:n(),
    ID = cut(1:n(),
      breaks = c(data$start, last(data$end)),
      labels = FALSE, right = FALSE, include.lowest = TRUE
    )
  ) %>%
  left_join(tibble(
    DV = data$cObs,
    i = data$iObs
  ),
  by = "i"
  ) %>%
  mutate(DV = if_else(is.na(DV),
    5555555,
    DV
  )) %>% # This value can be anything > 0. It'll be indexed away
  select(c, ID, everything(), -i)

n_subjects <- nonmem_data %>% # number of individuals
  distinct(ID) %>%
  count() %>%
  deframe()

n_total <- nrow(nonmem_data) # total number of records

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

## Fit the same model that Metrum fits (exponential error, lognormal priors).
## The few differences are:
## 1) Half-normal prior with no upper bound on omegas rather than uniform on
##        (0.01, 2)
## 2) No upper bound on the population parameters
##

stan_data <- with(
  nonmem_data,
  list(
    n_subjects = n_subjects,
    n_total = n_total,
    n_obs = n_obs,
    i_obs = i_obs,
    ID = ID,
    amt = amt,
    cmt = cmt,
    evid = evid,
    rate = rate,
    ii = ii,
    addl = addl,
    ss = ss,
    time = time,
    dv = DV,
    subj_start = subj_start,
    subj_end = subj_end,
    lloq = lloq,
    bloq = bloq,
    location_tvcl = 10,
    location_tvvc = 35,
    location_tvq = 15,
    location_tvvp = 105,
    location_tvka = 2.5,
    scale_tvcl = 0.25,
    scale_tvvc = 0.25,
    scale_tvq = 0.5,
    scale_tvvp = 0.5,
    scale_tvka = 1,
    scale_omega_cl = 0.4,
    scale_omega_vc = 0.4,
    scale_omega_q = 0.4,
    scale_omega_vp = 0.4,
    scale_omega_ka = 0.4,
    lkj_df_omega = 2,
    scale_sigma = 5
  )
)

model_normal <- cmdstan_model(
  "Torsten/example-models/poppk2cpt/depot_2cmt_match_metrum_half_normal_omega.stan", # nolint
  cpp_options = list(stan_threads = TRUE)
)

fit_normal <- model_normal$sample(
  data = stan_data,
  chains = 4,
  parallel_chains = 4,
  threads_per_chain = 8,
  iter_warmup = 1000,
  iter_sampling = 1000,
  adapt_delta = 0.8,
  refresh = 100,
  max_treedepth = 10,
  init = file.path(
    "Torsten", "example-models", "poppk2cpt",
    "poppk2cpt.init_cd.R"
  )
)

# ELPD
fit_normal_loo <- fit_normal$loo()
loo_df <- fit_normal_loo$pointwise %>%
  as_data_frame() %>%
  mutate(
    row_number = row_number()
  ) %>%
  left_join(
    nonmem_data %>%
      filter(evid == 0) %>%
      transmute(
        row_number = row_number(),
        ID, time, DV, bloq
      ),
    by = "row_number"
  ) %>%
  relocate(
    row_number,
    .before = elpd_loo
  )

loo_df %>% write_csv("Torsten-Runs/loo-CV/loo_poppk2cpt-reduce_sum.csv")


# Computed from 4000 by 530 log-likelihood matrix

#          Estimate   SE
# elpd_loo  -3196.0 22.8
# p_loo        15.6  1.2
# looic      6391.9 45.7
# ------
# Monte Carlo SE of elpd_loo is 0.1.

# All Pareto k estimates are good (k < 0.5).
# See help('pareto-k-diagnostic') for details.

model_uniform <- cmdstan_model(
  "Torsten/example-models/poppk2cpt/depot_2cmt_match_metrum_uniform_omega.stan",
  cpp_options = list(stan_threads = TRUE)
)

fit_uniform <- model_uniform$sample(
  data = stan_data,
  chains = 4,
  parallel_chains = 4,
  threads_per_chain = 8,
  iter_warmup = 1000,
  iter_sampling = 1000,
  adapt_delta = 0.8,
  refresh = 100,
  max_treedepth = 10,
  init = file.path(
    "Torsten", "example-models", "poppk2cpt",
    "poppk2cpt.init_cd.R"
  )
)

# ELPD
fit_uniform$loo()
