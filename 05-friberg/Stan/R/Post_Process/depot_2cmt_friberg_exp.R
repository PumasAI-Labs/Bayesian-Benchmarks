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
nonmem_data <- read_csv("05-friberg/data/multiple_dose.csv",
                        na = ".") %>%
  rename_all(tolower) %>%
  rename(ID = "id",
         DV = "dv") %>%
  mutate(DV = if_else(is.na(DV), 5555555, DV), # This value can be anything except NA. It'll be indexed away
         bloq = if_else(is.na(bloq), -999, bloq)) # This value can be anything except NA. It'll be indexed away

# (p_pk <- ggplot(nonmem_data %>%
#                   mutate(ID = factor(ID)) %>%
#                   group_by(ID) %>%
#                   mutate(Dose = factor(max(amt, na.rm = TRUE))) %>%
#                   ungroup() %>%
#                   filter(mdv == 0, cmt == 2)) +
#     geom_line(mapping = aes(x = time, y = DV, group = ID, color = Dose)) +
#     geom_point(mapping = aes(x = time, y = DV, group = ID, color = Dose)) +
#     scale_color_discrete(name = "Dose (ug)") +
#     scale_y_continuous(name = "Drug Conc. (ng/mL)",
#                        limits = c(NA, NA),
#                        trans = "log10") +
#     scale_x_continuous(name = "Time (d)",
#                        breaks = seq(0, 216, by = 24),
#                        labels = seq(0, 216/24, by = 24/24),
#                        limits = c(0, NA)) +
#     theme_bw(18) +
#     theme(axis.text = element_text(size = 14, face = "bold"),
#           axis.title = element_text(size = 18, face = "bold"),
#           axis.line = element_line(size = 2),
#           legend.position = "bottom"))
# 
# (p_pd <- ggplot(nonmem_data %>%
#                   mutate(ID = factor(ID)) %>%
#                   group_by(ID) %>%
#                   mutate(Dose = factor(max(amt, na.rm = TRUE))) %>%
#                   ungroup() %>%
#                   filter(mdv == 0, cmt == 3)) +
#     geom_line(mapping = aes(x = time, y = DV, group = ID, color = Dose)) +
#     geom_point(mapping = aes(x = time, y = DV, group = ID, color = Dose)) +
#     scale_color_discrete(name = "Dose (ug)") +
#     scale_y_continuous(name = latex2exp::TeX("$Neutrophils\\;(10^9/L)"),
#                        limits = c(NA, NA),
#                        trans = "identity") +
#     scale_x_continuous(name = "Time (d)",
#                        breaks = seq(0, 672, by = 24),
#                        labels = seq(0, 672/24, by = 24/24),
#                        limits = c(0, NA)) +
#     theme_bw(18) +
#     theme(axis.text = element_text(size = 14, face = "bold"),
#           axis.title = element_text(size = 18, face = "bold"),
#           axis.line = element_line(size = 2),
#           legend.position = "bottom"))
# 
# p_pk +
#   p_pd +
#   plot_layout(guides = "collect") &
#   theme(legend.position = "bottom")
# 
# nonmem_data %>%
#   filter(evid == 0) %>% 
#   group_by(ID, cmt) %>% 
#   summarize(lloq = unique(lloq),
#             n_obs = n(),
#             n_bloq = sum(bloq == 1)) %>% 
#   filter(n_bloq > 0)

## Read in fit
fit <- read_rds("05-friberg/Stan/Torsten/Fits/multiple_dose_1.rds")


## Summary of parameter estimates
parameters_to_summarize <- c(str_subset(fit$metadata()$stan_variables, "TV"),
                             str_subset(fit$metadata()$stan_variables, "omega"),
                             str_subset(fit$metadata()$stan_variables, "sigma"))

summary <- summarize_draws(fit$draws(parameters_to_summarize),
                           mean, median, sd, mcse_mean,
                           ~ quantile2(.x, probs = c(0.025, 0.975)), rhat,
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
