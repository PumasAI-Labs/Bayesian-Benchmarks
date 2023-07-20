rm(list = ls())
cat("\014")

library(trelliscopejs)
library(cmdstanr)
library(tidybayes)
library(posterior)
library(tidyverse)

set_cmdstan_path("cmdstan")

# fit <- read_rds("05-friberg/Stan/Torsten/Fits/multiple_dose_1.rds")
fit <- read_rds("05-friberg/Stan/Torsten/Fits/multiple_dose_coupled_1.rds")

nonmem_data <- read_csv("05-friberg/data/multiple_dose.csv",
                        na = ".") %>% 
  rename_all(tolower) %>% 
  rename(ID = "id",
         DV = "dv") %>% 
  mutate(DV = if_else(is.na(DV), 5555555, DV),    # This value can be anything except NA. It'll be indexed away 
         bloq = if_else(is.na(bloq), -999, bloq)) # This value can be anything except NA. It'll be indexed away 

new_data_to_simulate_pk <- nonmem_data %>% 
  filter(cmt == 2 | evid == 1) %>%
  group_by(ID) %>% 
  slice(c(1, n())) %>% 
  expand(time = seq(time[1], time[2], by = 0.5)) %>%
  # expand(time = seq(time[1], time[2], by = 2)) %>%
  ungroup() %>% 
  mutate(amt = 0,
         evid = 2,
         rate = 0,
         addl = 0,
         ii = 0,
         cmt = 2,
         mdv = 1, 
         ss = 0,
         lloq = 0.0001,
         bloq = 0,
         DV = NA_real_,
         c = NA_character_) %>% 
  select(c, ID, time, everything())

new_data_to_simulate_pd <- nonmem_data %>% 
  filter(cmt == 3 | evid == 1) %>%
  group_by(ID) %>% 
  slice(c(1, n())) %>% 
  expand(time = seq(time[1], time[2], by = 1)) %>%
  # expand(time = seq(time[1], time[2], by = 2)) %>%
  ungroup() %>% 
  mutate(amt = 0,
         evid = 2,
         rate = 0,
         addl = 0,
         ii = 0,
         cmt = 3,
         mdv = 1, 
         ss = 0,
         lloq = 0.0001,
         bloq = 0,
         DV = NA_real_,
         c = NA_character_) %>% 
  select(c, ID, time, everything())

observed_data_to_simulate <- nonmem_data %>% 
  mutate(evid = if_else(evid == 1, 1, 2))

new_data <- new_data_to_simulate_pk %>% 
  bind_rows(new_data_to_simulate_pd,
            observed_data_to_simulate) %>% 
  arrange(ID, time, cmt, evid) %>% 
  distinct(ID, time, cmt, .keep_all = TRUE) %>% 
  select(-DV, -lloq, -bloq)


n_subjects <- new_data %>%  # number of inddepotiduals
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
  "05-friberg/Stan/Torsten/Predict/depot_2cmt_friberg_exp_predict_observed_subjects.stan")

preds <- model$generate_quantities(fit,
                                   data = stan_data,
                                   parallel_chains = 4,
                                   seed = 1234) 

preds_df <- preds$draws(format = "draws_df")

rm(list = setdiff(ls(), c("preds_df", "new_data", "nonmem_data")))
gc()

post_preds_summary <- preds_df %>%
  spread_draws(pred[i], ipred[i], dv[i]) %>%
  mean_qi(pred, ipred, dv) 

preds_ind <- post_preds_summary %>% 
  mutate(ID = new_data$ID[i],
         time = new_data$time[i],
         cmt = new_data$cmt[i]) %>% 
  left_join(nonmem_data %>% 
              filter(mdv == 0) %>%
              select(ID, time, cmt, DV), 
            by = c("ID", "time", "cmt")) %>% 
  select(ID, time, cmt, everything(), -i)

ggplot(preds_ind %>% 
         filter(cmt == 2), 
       aes(x = time, group = ID)) +
  geom_ribbon(aes(ymin = dv.lower, ymax = dv.upper),
              fill = "blue", alpha = 0.25, show.legend = FALSE) +
  geom_ribbon(aes(ymin = ipred.lower, ymax = ipred.upper),
              fill = "blue", alpha = 0.5, show.legend = FALSE) +
  geom_line(aes(y = ipred), linetype = 1, linewidth = 1.15) +
  geom_line(aes(y = pred), linetype = 2, linewidth = 1.05) +
  geom_point(aes(y = DV), size = 2, show.legend = FALSE, 
             color = "red") +
  scale_y_continuous(name = latex2exp::TeX("Drug Conc. $(\\mu g/mL)$"),
                     trans = "identity",
                     limits = c(NA, NA)) +
  scale_x_continuous(name = "Time (h)",
                     breaks = seq(0, 216, by = 24),
                     labels = seq(0, 216, by = 24),
                     limits = c(0, NA)) +
  theme_bw() +
  theme(axis.text = element_text(size = 14, face = "bold"),
        axis.title = element_text(size = 18, face = "bold"),
        legend.position = "bottom") +
  facet_trelliscope(~ ID, nrow = 3, ncol = 4, scales = "free")

ggplot(preds_ind %>% 
         filter(cmt == 3), 
       aes(x = time, group = ID)) +
  geom_ribbon(aes(ymin = dv.lower, ymax = dv.upper),
              fill = "blue", alpha = 0.25, show.legend = FALSE) +
  geom_ribbon(aes(ymin = ipred.lower, ymax = ipred.upper),
              fill = "blue", alpha = 0.5, show.legend = FALSE) +
  geom_line(aes(y = ipred), linetype = 1, linewidth = 1.15) +
  geom_line(aes(y = pred), linetype = 2, linewidth = 1.05) +
  geom_point(aes(y = DV), size = 2, show.legend = FALSE, 
             color = "red") +
  scale_y_continuous(name = latex2exp::TeX("$Neutrophils\\;(10^9/L)"),
                     trans = "identity",
                     limits = c(NA, NA)) +
  scale_x_continuous(name = "Time (h)",
                     breaks = seq(0, 672, by = 24),
                     labels = seq(0, 672, by = 24),
                     limits = c(0, NA)) +
  theme_bw() +
  theme(axis.text = element_text(size = 14, face = "bold"),
        axis.title = element_text(size = 18, face = "bold"),
        legend.position = "bottom") +
  facet_trelliscope(~ ID, nrow = 3, ncol = 4, scales = "free")


tmp <- ggplot(preds_ind, aes(x = time, group = ID)) +
  geom_line(aes(y = ipred), linetype = 1, linewidth = 1.15) +
  ggforce::facet_wrap_paginate(~ID, labeller = label_both,
                               nrow = 1, ncol = 1,
                               page = 1, scales = "free")

for(i in 1:ggforce::n_pages(tmp)){
  print(ggplot(preds_ind %>% 
                 filter(cmt == 2), 
               aes(x = time, group = ID)) +
          geom_ribbon(aes(ymin = dv.lower, ymax = dv.upper),
                      fill = "blue", alpha = 0.25, show.legend = FALSE) +
          geom_ribbon(aes(ymin = ipred.lower, ymax = ipred.upper),
                      fill = "blue", alpha = 0.5, show.legend = FALSE) +
          geom_line(aes(y = ipred), linetype = 1, linewidth = 1.15) +
          geom_line(aes(y = pred), linetype = 2, linewidth = 1.05) +
          geom_point(aes(y = DV), size = 2, show.legend = FALSE, 
                     color = "red") +
          geom_point(data = nonmem_data %>% 
                       filter(bloq == 1),
                     mapping = aes(x = time, y = lloq, group = ID), 
                     color = "limegreen") +
          scale_y_continuous(name = latex2exp::TeX("Drug Conc. $(\\mu g/mL)$"),
                             trans = "log10",
                             limits = c(NA, NA)) +
          scale_x_continuous(name = "Time (h)",
                             breaks = seq(0, 216, by = 24),
                             labels = seq(0, 216, by = 24),
                             limits = c(0, NA)) +
          theme_bw() +
          theme(axis.text = element_text(size = 14, face = "bold"),
                axis.title = element_text(size = 18, face = "bold"),
                legend.position = "bottom") +
          ggforce::facet_wrap_paginate(~ ID, labeller = label_both,
                                       nrow = 1, ncol = 1,
                                       page = i, scales = "free"))
  
}

for(i in 1:ggforce::n_pages(tmp)){
  print(ggplot(preds_ind %>% 
                 filter(cmt == 3), 
               aes(x = time, group = ID)) +
          geom_ribbon(aes(ymin = dv.lower, ymax = dv.upper),
                      fill = "blue", alpha = 0.25, show.legend = FALSE) +
          geom_ribbon(aes(ymin = ipred.lower, ymax = ipred.upper),
                      fill = "blue", alpha = 0.5, show.legend = FALSE) +
          geom_line(aes(y = ipred), linetype = 1, linewidth = 1.15) +
          geom_line(aes(y = pred), linetype = 2, linewidth = 1.05) +
          geom_point(aes(y = DV), size = 2, show.legend = FALSE, 
                     color = "red") +
          geom_point(data = nonmem_data %>% 
                       filter(bloq == 1),
                     mapping = aes(x = time, y = lloq, group = ID), 
                     color = "limegreen") +
          scale_y_continuous(name = latex2exp::TeX("$Neutrophils\\;(10^9/L)"),
                             trans = "identity",
                             limits = c(NA, NA)) +
          scale_x_continuous(name = "Time (d)",
                             breaks = seq(0, 672, by = 24),
                             labels = seq(0, 672/24, by = 24/24),
                             limits = c(0, NA)) +
          theme_bw() +
          theme(axis.text = element_text(size = 14, face = "bold"),
                axis.title = element_text(size = 18, face = "bold"),
                legend.position = "bottom") +
          ggforce::facet_wrap_paginate(~ ID, labeller = label_both,
                                       nrow = 1, ncol = 1,
                                       page = i, scales = "free"))
  
}
