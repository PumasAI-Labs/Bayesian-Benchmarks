rm(list = ls())
cat("\014")

library(cmdstanr)
library(tidybayes)
library(posterior)
library(tidyverse)

set_cmdstan_path("cmdstan")


fit <- read_rds("05-friberg/Stan/Torsten/Fits/multiple_dose.rds")

# For this example, let's simulate 10 mg, 20 mg, 40 mg, 80 mg, 160 mg, 320 mg 
dosing_data <- mrgsolve::expand.ev(addl = 6, ii = 24, cmt = 1,
                                   amt = c(10, 20, 40, 80, 160, 320)*1000, 
                                   tinf = 0, evid = 1, mdv = 1) %>%
  as_tibble() %>% 
  mutate(ss = 0) %>% 
  select(ID, time, everything())

t1 <- dosing_data %>% 
  select(time) %>% 
  distinct() %>% 
  deframe()

times_new_pk <- tibble(time = sort(unique(c(t1, 0.25, seq(0, 216, by = 0.5)))))
times_new_pd <- tibble(time = sort(unique(c(t1, 1, seq(1, 672, by = 1)))))

new_data_pk <- bind_rows(replicate(max(dosing_data$ID), times_new_pk, 
                                simplify = FALSE)) %>% 
  mutate(ID = rep(1:max(dosing_data$ID), each = nrow(times_new_pk)),
         amt = 0, 
         evid = 0, 
         rate = 0, 
         addl = 0, 
         ii = 0, 
         cmt = 2, 
         mdv = 1, 
         ss = 0) %>%
  filter(time != 0) %>% 
  select(ID, time, everything()) 

new_data_pd <- bind_rows(replicate(max(dosing_data$ID), times_new_pd, 
                                simplify = FALSE)) %>% 
  mutate(ID = rep(1:max(dosing_data$ID), each = nrow(times_new_pd)),
         amt = 0, 
         evid = 0, 
         rate = 0, 
         addl = 0, 
         ii = 0, 
         cmt = 3, 
         mdv = 1, 
         ss = 0) %>%
  filter(time != 0) %>% 
  select(ID, time, everything())


new_data <- new_data_pk %>% 
  bind_rows(new_data_pd, dosing_data) %>% 
  arrange(ID, time, cmt)


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
  "05-friberg/Stan/Torsten/Predict/depot_2cmt_friberg_prop_predict_new_subjects.stan")

preds <- model$generate_quantities(fit,
                                   data = stan_data,
                                   parallel_chains = 4,
                                   seed = 1234) 

preds_df <- preds$draws(format = "draws_df")



# post_preds_summary <- preds_df %>%
#   spread_draws(pred[i], ipred[i], dv[i]) %>%
#   mean_qi(pred, ipred, dv) 
# 
# preds_ind <- post_preds_summary %>% 
#   mutate(ID = new_data$ID[i],
#          time = new_data$time[i],
#          cmt = new_data$cmt[i]) %>% 
#   left_join(nonmem_data %>% 
#               filter(mdv == 0) %>%
#               select(ID, time, cmt, DV), 
#             by = c("ID", "time", "cmt")) %>% 
#   select(ID, time, cmt, everything(), -i)



post_preds_summary <- preds_df %>%
  spread_draws(cp[i], dv[i]) %>%
  median_qi(cp, dv)

regimens <- str_c(c(10, 20, 40, 80, 160, 320), " mg")

preds_ind <- post_preds_summary %>%
  mutate(ID = new_data$ID[i],
         time = new_data$time[i],
         cmt = new_data$cmt[i]) %>%
  select(ID, time, cmt, everything(), -i) %>% 
  mutate(regimen = factor(regimens[ID],
                          levels = regimens))

tmp <- ggplot(preds_ind, aes(x = time, group = ID)) +
  geom_line(aes(y = cp), linetype = 1, size = 1.15) +
  geom_line(aes(y = dv), linetype = 2, size = 1.05) +
  ggforce::facet_wrap_paginate(~ID, labeller = label_both,
                               nrow = 3, ncol = 2,
                               page = 1)

for(i in 1:ggforce::n_pages(tmp)){
  print(ggplot(preds_ind %>% 
                 filter(cmt == 2), 
               aes(x = time, group = ID)) +
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
                             limits = c(0, NA)) +
          # scale_x_continuous(name = "Time (d)") +
          theme_bw() +
          theme(axis.text = element_text(size = 14, face = "bold"),
                axis.title = element_text(size = 18, face = "bold"),
                legend.position = "bottom") +
          ggforce::facet_wrap_paginate(~ regimen,
                                       nrow = 1, ncol = 6,
                                       page = i))
  
}

for(i in 1:ggforce::n_pages(tmp)){
  print(ggplot(preds_ind %>% 
                 filter(cmt == 3), 
               aes(x = time, group = ID)) +
          geom_ribbon(aes(ymin = dv.lower, ymax = dv.upper),
                      fill = "blue", alpha = 0.25, show.legend = FALSE) +
          geom_ribbon(aes(ymin = cp.lower, ymax = cp.upper),
                      fill = "blue", alpha = 0.5, show.legend = FALSE) +
          geom_line(aes(y = cp), linetype = 1, size = 1.15) +
          geom_line(aes(y = dv), linetype = 2, size = 1.05) +
          scale_y_continuous(name = latex2exp::TeX("$Neutrophils\\;(10^9/L)"),
                             trans = "log10",
                             limits = c(NA, NA)) +
          scale_x_continuous(name = "Time (d)",
                             breaks = seq(0, 672, by = 24),
                             labels = seq(0, 672/24, by = 24/24),
                             limits = c(0, NA)) +
          # scale_x_continuous(name = "Time (d)") +
          theme_bw() +
          theme(axis.text = element_text(size = 14, face = "bold"),
                axis.title = element_text(size = 18, face = "bold"),
                legend.position = "bottom") +
          ggforce::facet_wrap_paginate(~ regimen,
                                       nrow = 1, ncol = 6,
                                       page = i))
  
}


data <- read_csv("05-friberg/data/multiple_dose.csv",
                 na = ".") %>% 
  rename_all(tolower) %>% 
  rename(ID = "id",
         DV = "dv") %>% 
  group_by(ID) %>% 
  mutate(Dose = max(amt, na.rm = TRUE)/1000,
         mdv = evid) %>% 
  ungroup() %>% 
  mutate(regimen = str_c(Dose, " mg"))


preds_ind %>% 
  filter(regimen %in% str_c(c(10, 20, 40, 80), " mg"), 
         cmt == 2) %>% 
  ggplot(aes(x = time, group = ID)) +
  geom_ribbon(aes(ymin = dv.lower, ymax = dv.upper),
              fill = "blue", alpha = 0.25, show.legend = FALSE) +
  geom_ribbon(aes(ymin = cp.lower, ymax = cp.upper),
              fill = "blue", alpha = 0.5, show.legend = FALSE) +
  geom_line(aes(y = cp), linetype = 1, size = 1.15) +
  geom_line(aes(y = dv), linetype = 2, size = 1.05) +
  geom_point(data = data %>% 
               filter(regimen %in% str_c(c(10, 20, 40, 80), " mg"), 
                      cmt == 2,
                      mdv == 0) %>% 
               mutate(regimen = fct_inorder(regimen)),
             mapping = aes(x = time, y = DV), color = "red", 
             inherit.aes = FALSE) +
  scale_y_continuous(name = "Drug Conc. (ug/mL)",
                     limits = c(NA, NA),
                     trans = "log10") +
  scale_x_continuous(name = "Time (h)",
                     breaks = seq(0, 216, by = 24),
                     labels = seq(0, 216, by = 24),
                     limits = c(0, NA)) +
  theme_bw() +
  theme(axis.text = element_text(size = 14, face = "bold"),
        axis.title = element_text(size = 18, face = "bold"),
        legend.position = "bottom") +
  facet_wrap(~ regimen, nrow = 2, ncol = 2)

preds_ind %>% 
  filter(regimen %in% str_c(c(10, 20, 40, 80), " mg"), 
         cmt == 3) %>% 
  ggplot(aes(x = time, group = ID)) +
  geom_ribbon(aes(ymin = dv.lower, ymax = dv.upper),
              fill = "blue", alpha = 0.25, show.legend = FALSE) +
  geom_ribbon(aes(ymin = cp.lower, ymax = cp.upper),
              fill = "blue", alpha = 0.5, show.legend = FALSE) +
  geom_line(aes(y = cp), linetype = 1, size = 1.15) +
  geom_line(aes(y = dv), linetype = 2, size = 1.05) +
  geom_point(data = data %>% 
               filter(regimen %in% str_c(c(10, 20, 40, 80), " mg"), 
                      cmt == 3,
                      mdv == 0) %>% 
               mutate(regimen = fct_inorder(regimen)),
             mapping = aes(x = time, y = DV), color = "red", 
             inherit.aes = FALSE) +
  scale_y_continuous(name = "Drug Conc. (ug/mL)",
                     limits = c(NA, NA),
                     trans = "log10") +
  scale_x_continuous(name = "Time (h)",
                     breaks = seq(0, 216, by = 24),
                     labels = seq(0, 216, by = 24),
                     limits = c(0, NA)) +
  theme_bw() +
  theme(axis.text = element_text(size = 14, face = "bold"),
        axis.title = element_text(size = 18, face = "bold"),
        legend.position = "bottom") +
  facet_wrap(~ regimen, nrow = 2, ncol = 2, scales = "free_y")
