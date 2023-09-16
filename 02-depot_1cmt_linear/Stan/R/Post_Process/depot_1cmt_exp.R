rm(list = ls())
cat("\014")

library(gganimate)
library(loo)
library(cmdstanr)
library(tidybayes)
library(bayesplot)
library(posterior)
library(tidyverse)

## Read in and visualize observed data
nonmem_data <- read_csv("02-depot_1cmt_linear/data/single_dose.csv",
# nonmem_data <- read_csv("02-depot_1cmt_linear/data/multiple_dose.csv",
                        na = ".") %>% 
  rename_all(tolower) %>% 
  rename(ID = "id",
         DV = "dv") %>% 
  mutate(DV = if_else(is.na(DV), 5555555, DV),    # This value can be anything except NA. It'll be indexed away 
         bloq = if_else(is.na(bloq), -999, bloq)) # This value can be anything except NA. It'll be indexed away 

# (p1 <- ggplot(nonmem_data %>%
#                 mutate(ID = factor(ID)) %>%
#                 group_by(ID) %>%
#                 mutate(Dose = factor(max(amt, na.rm = TRUE))) %>%
#                 ungroup() %>%
#                 filter(mdv == 0)) +
#     geom_line(mapping = aes(x = time, y = DV, group = ID, color = Dose)) +
#     geom_point(mapping = aes(x = time, y = DV, group = ID, color = Dose)) +
#     scale_color_discrete(name = "Dose (mg)") +
#     scale_y_log10(name = latex2exp::TeX("$Drug Conc. \\; (\\mu g/mL)$"),
#                   limits = c(NA, NA)) +
#     scale_x_continuous(name = "Time (h)",
#                        breaks = seq(0, 24, by = 4),
#                        labels = seq(0, 24, by = 4)) +
#     theme_bw(18) +
#     theme(axis.text = element_text(size = 14, face = "bold"),
#           axis.title = element_text(size = 18, face = "bold"),
#           axis.line = element_line(size = 2),
#           legend.position = "bottom") + 
#     geom_hline(aes(yintercept = lloq, group = ID), size = 1.2, linetype = 2) +
#     geom_vline(data = nonmem_data %>%
#                  mutate(ID = factor(ID)) %>%
#                  filter(evid == 1),
#                mapping = aes(xintercept = time, group = ID), 
#                linetype = 2, color = "magenta", alpha = 0.5) +   
#     geom_point(data = nonmem_data %>%
#                  mutate(ID = factor(ID),
#                         DV = lloq) %>%
#                  filter(bloq == 1, time > 0),
#                aes(x = time, y = DV, group = ID), inherit.aes = TRUE,
#                shape = 18, color = "limegreen", size = 3) + 
#     facet_wrap(~ID, labeller = label_both, ncol = 3, scales = "free_y")) 
# 
# 
# nonmem_data %>%
#   filter(evid == 0) %>% 
#   group_by(ID) %>% 
#   summarize(lloq = unique(lloq),
#             n_obs = n(),
#             n_bloq = sum(bloq)) %>% 
#   print(n = 100)


## Read in fit
fit <- read_rds("02-depot_1cmt_linear/Stan/Torsten/Fits/single_dose_1.rds")
# fit <- read_rds("02-depot_1cmt_linear/Stan/Torsten/Fits/multiple_dose_1.rds")


## Summary of parameter estimates
parameters_to_summarize <- c(str_subset(fit$metadata()$stan_variables, "TV"),
                             str_subset(fit$metadata()$stan_variables, "omega"),
                             "sigma")

summary <- summarize_draws(fit$draws(parameters_to_summarize), 
                           mean, median, sd, mcse_mean,
                           ~quantile2(.x, probs = c(0.025, 0.975)), rhat,
                           ess_bulk, ess_tail)

summary %>%
  mutate(rse = sd/mean*100) %>% 
  select(variable, mean, sd, q2.5, median, q97.5, rse, rhat, 
         starts_with("ess")) %>% 
  print(n = 50)


## Check the sampler (this is very non-comprehensive)
sampler_diagnostics <- fit$sampler_diagnostics()
sum(as_draws_df(sampler_diagnostics)$divergent__)

fit$time()
