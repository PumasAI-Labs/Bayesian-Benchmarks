rm(list = ls())
cat("\014")

library(cmdstanr)
library(tidybayes)
library(posterior)
library(tidyverse)

set_cmdstan_path("Torsten/cmdstan")

# fit <- read_rds("Model2/depot_2cmt_linear/Stan/Torsten/Fits/single_dose.rds")
fit <- read_rds("Model2/depot_2cmt_linear/Stan/Torsten/Fits/multiple_dose.rds")

# For this example, let's simulate 10 mg, 30 mg, 60 mg, 120 mg, 240 mg, 480 mg 
# dosing_data <- mrgsolve::expand.ev(addl = 0, ii = 0, cmt = 1,
dosing_data <- mrgsolve::expand.ev(addl = 6, ii = 24, cmt = 1,
                                   amt = c(10, 30, 60, 120, 240, 480), tinf = 0,
                                   evid = 1, mdv = 1) %>%
  as_tibble() %>%
  mutate(ss = 0) %>%
  select(ID, time, everything())

t1 <- dosing_data %>% 
  select(time) %>% 
  distinct() %>% 
  deframe()

# times_new <- tibble(time = sort(unique(c(t1, 0.25, seq(0, 72, by = 0.5)))))
times_new <- tibble(time = sort(unique(c(t1, 0.25, seq(0, 216, by = 0.5)))))

new_data <- bind_rows(replicate(max(dosing_data$ID), times_new, 
                                simplify = FALSE)) %>% 
  mutate(ID = rep(1:max(dosing_data$ID), each = nrow(times_new)),
         amt = 0, 
         evid = 0, 
         rate = 0, 
         addl = 0, 
         ii = 0, 
         cmt = 2, 
         mdv = 1, 
         ss = 0) %>%
  filter(time != 0) %>% 
  select(ID, time, everything()) %>% 
  bind_rows(dosing_data) %>% 
  arrange(ID, time)


# number of individuals in the original dataset
n_subjects <- fit$metadata()$stan_variable_sizes$Z[2]

n_subjects_new <- new_data %>%  # number of new individuals
  distinct(ID) %>% 
  count() %>% 
  deframe()

n_time_new <- nrow(new_data) # total number of time points at which to predict

subj_start <- new_data %>% 
  mutate(row_num = 1:n()) %>% 
  group_by(ID) %>% 
  slice_head(n = 1) %>%
  ungroup() %>% 
  select(row_num) %>% 
  deframe()

subj_end <- c(subj_start[-1] - 1, n_time_new)  

stan_data <- list(n_subjects = n_subjects,
                  n_subjects_new = n_subjects_new,
                  n_time_new = n_time_new,
                  time = new_data$time,
                  amt = new_data$amt,
                  cmt = new_data$cmt,
                  evid = new_data$evid,
                  rate = new_data$rate,
                  ii = new_data$ii,
                  addl = new_data$addl,
                  ss = new_data$ss,
                  subj_start = subj_start,
                  subj_end = subj_end)

model <- cmdstan_model(
  "Model2/depot_2cmt_linear/Stan/Torsten/Predict/depot_2cmt_prop_predict_new_subjects.stan")

preds <- model$generate_quantities(fit,
                                   data = stan_data,
                                   parallel_chains = 4,
                                   seed = 1234) 

# preds$save_object("Torsten/Preds/iv_2cmt_ppa_m4_predict_new_subjects.rds")

preds_df <- preds$draws(format = "draws_df")

post_preds_summary <- preds_df %>%
  spread_draws(cp[i], dv[i]) %>%
  median_qi(cp, dv)

regimens <- str_c(c(10, 30, 60, 120, 240, 480), " mg")

preds_ind <- post_preds_summary %>%
  mutate(ID = new_data$ID[i],
         time = new_data$time[i]) %>%
  select(ID, time, everything(), -i) %>% 
  mutate(regimen = factor(regimens[ID],
                          levels = regimens))

tmp <- ggplot(preds_ind, aes(x = time, group = ID)) +
  geom_line(aes(y = cp), linetype = 1, size = 1.15) +
  geom_line(aes(y = dv), linetype = 2, size = 1.05) +
  ggforce::facet_wrap_paginate(~ID, labeller = label_both,
                               nrow = 3, ncol = 2,
                               page = 1)

for(i in 1:ggforce::n_pages(tmp)){
  print(ggplot(preds_ind, aes(x = time, group = ID)) +
          geom_ribbon(aes(ymin = dv.lower, ymax = dv.upper),
                      fill = "blue", alpha = 0.25, show.legend = FALSE) +
          geom_ribbon(aes(ymin = cp.lower, ymax = cp.upper),
                      fill = "blue", alpha = 0.5, show.legend = FALSE) +
          geom_line(aes(y = cp), linetype = 1, size = 1.15) +
          geom_line(aes(y = dv), linetype = 2, size = 1.05) +
          scale_y_continuous(name = "Drug Conc. (ug/mL)",
                             trans = "log10",
                             limits = c(NA, NA)) +
          scale_x_continuous(name = "Time (h)",
                             breaks = seq(0, 216, by = 24),
                             labels = seq(0, 216, by = 24),
                             limits = c(0, 216)) +
          # scale_x_continuous(name = "Time (d)") +
          theme_bw() +
          theme(axis.text = element_text(size = 14, face = "bold"),
                axis.title = element_text(size = 18, face = "bold"),
                legend.position = "bottom") +
          ggforce::facet_wrap_paginate(~ regimen,
                                       nrow = 3, ncol = 2,
                                       page = i))
  
}

# data <- read_csv("Model2/depot_2cmt_linear/Data/single_dose.csv",
data <- read_csv("Model2/depot_2cmt_linear/Data/multiple_dose.csv",
                 na = ".") %>% 
  rename_all(tolower) %>% 
  rename(ID = "id",
         DV = "dv") %>% 
  group_by(ID) %>% 
  mutate(Dose = max(amt, na.rm = TRUE),
         mdv = evid) %>% 
  ungroup() %>% 
  mutate(regimen = str_c(Dose, " mg"))


preds_ind %>% 
  filter(regimen %in% str_c(c(10, 30, 60, 120), " mg")) %>% 
  ggplot(aes(x = time, group = ID)) +
  geom_ribbon(aes(ymin = dv.lower, ymax = dv.upper),
              fill = "blue", alpha = 0.25, show.legend = FALSE) +
  geom_ribbon(aes(ymin = cp.lower, ymax = cp.upper),
              fill = "blue", alpha = 0.5, show.legend = FALSE) +
  geom_line(aes(y = cp), linetype = 1, size = 1.15) +
  geom_line(aes(y = dv), linetype = 2, size = 1.05) +
  geom_point(data = data %>% 
               filter(regimen %in% str_c(c(10, 30, 60, 120), " mg"), 
                      mdv == 0),
             mapping = aes(x = time, y = DV), color = "red", 
             inherit.aes = FALSE) +
  scale_y_log10(name = "Drug Conc. (ug/mL)",
                limits = c(NA, NA)) +
  scale_x_continuous(name = "Time (h)",
                     breaks = seq(0, 216, by = 24),
                     labels = seq(0, 216, by = 24),
                     limits = c(0, 216)) +
  theme_bw() +
  theme(axis.text = element_text(size = 14, face = "bold"),
        axis.title = element_text(size = 18, face = "bold"),
        legend.position = "bottom") +
  facet_wrap(~ regimen, scales = "free_y")

