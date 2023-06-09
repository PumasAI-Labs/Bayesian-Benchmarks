rm(list = ls())
cat("\014")

library(cmdstanr)
library(tidyverse)

set_cmdstan_path("cmdstan")

nonmem_data <- read_csv("02-depot_1cmt_linear/data/single_dose.csv",
# nonmem_data <- read_csv("02-depot_1cmt_linear/data/multiple_dose.csv",
                        na = ".") %>% 
  rename_all(tolower) %>% 
  rename(ID = "id",
         DV = "dv") %>% 
  mutate(DV = if_else(is.na(DV), 5555555, DV),    # This value can be anything except NA. It'll be indexed away 
         bloq = if_else(is.na(bloq), -999, bloq)) # This value can be anything except NA. It'll be indexed away 


(p1 <- ggplot(nonmem_data %>%
                mutate(ID = factor(ID)) %>%
                group_by(ID) %>%
                mutate(Dose = factor(max(amt, na.rm = TRUE))) %>%
                ungroup() %>%
                filter(mdv == 0)) +
    geom_line(mapping = aes(x = time, y = DV, group = ID, color = Dose)) +
    geom_point(mapping = aes(x = time, y = DV, group = ID, color = Dose)) +
    scale_color_discrete(name = "Dose (mg)") +
    scale_y_continuous(name = "Drug Conc. (ug/mL)",
                       limits = c(NA, NA),
                       trans = "identity") +
    scale_x_continuous(name = "Time (d)",
                       breaks = seq(0, 216, by = 24),
                       labels = seq(0, 216/24, by = 24/24),
                       limits = c(0, NA)) +
    theme_bw(18) +
    theme(axis.text = element_text(size = 14, face = "bold"),
          axis.title = element_text(size = 18, face = "bold"),
          axis.line = element_line(size = 2),
          legend.position = "bottom"))

# p1 +
#   geom_vline(data = nonmem_data %>%
#                mutate(ID = factor(ID)) %>%
#                filter(evid == 1),
#              mapping = aes(xintercept = time, group = ID),
#              linetype = 2, color = "red", alpha = 0.5) +
#   facet_wrap(~ID, labeller = label_both, scales = "free_y")

nonmem_data %>%
  filter(evid== 0) %>% 
  group_by(ID) %>% 
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
                  location_tvka = 1,
                  scale_tvcl = 1,
                  scale_tvvc = 1,
                  scale_tvka = 1,
                  scale_omega_cl = 0.4,
                  scale_omega_vc = 0.4,
                  scale_omega_ka = 0.4,
                  lkj_df_omega = 2,
                  scale_sigma_p = 0.5)

model <- cmdstan_model(
  "02-depot_1cmt_linear/Stan/Torsten/Fit/depot_1cmt_prop_mat_exp.stan",
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
                    init = str_c("02-depot_1cmt_linear/data/inits/inits_1_", 
                                 1:4, ".json"))

fit$save_object("02-depot_1cmt_linear/Stan/Torsten/Fits/single_dose_mat_exp.rds")
# fit$save_object("02-depot_1cmt_linear/Stan/Torsten/Fits/multiple_dose_mat_exp.rds")


parameters_to_summarize <- c(str_subset(fit$metadata()$stan_variables, "TV"),
                             str_subset(fit$metadata()$stan_variables, "omega"),
                             "sigma_p")

fit$draws(parameters_to_summarize, format = "draws_df") %>% 
  as_tibble() %>% 
  write_csv("02-depot_1cmt_linear/Stan/Torsten/Fits/single_dose_draws_df_mat_exp.csv")
# write_csv("02-depot_1cmt_linear/Stan/Torsten/Fits/multiple_dose_draws_df_mat_exp.csv")