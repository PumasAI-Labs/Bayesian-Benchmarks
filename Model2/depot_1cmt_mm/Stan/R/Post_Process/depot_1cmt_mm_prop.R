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
# nonmem_data <- read_csv("Model2/depot_1cmt_mm/Data/single_dose.csv",
nonmem_data <- read_csv("Model2/iv_1cmt_linear/Data/multiple_dose.csv",
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

nonmem_data %>%	
  filter(evid== 0) %>% 	
  group_by(ID) %>% 	
  summarize(lloq = unique(lloq),	
            n_obs = n(),	
            n_bloq = sum(bloq == 1)) %>% 	
  filter(n_bloq > 0)


## Read in fit
# fit <- read_rds("Model2/depot_1cmt_mm/Stan/Torsten/Fits/single_dose.rds")
fit <- read_rds("Model2/depot_1cmt_mm/Stan/Torsten/Fits/multiple_dose.rds")


## Summary of parameter estimates
parameters_to_summarize <- c(str_subset(fit$metadata()$stan_variables, "TV"),
                             str_c("omega_", c("vc", "vmax", "km", "ka")),
                             "sigma_p")

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

# Check your largest R-hats. These should be close to 1 (< 1.05 for sure, < 1.02
# ideally)
summary %>%
  arrange(-rhat)

# Density Plots and Traceplots
mcmc_combo(fit$draws(c("TVVC", "TVVMAX", "TVKM", "TVKA")),
           combo = c("dens_overlay", "trace"))
mcmc_combo(fit$draws(c("omega_vc", "omega_vmax", "omega_km", "omega_ka")),
           combo = c("dens_overlay", "trace"))
mcmc_combo(fit$draws(c("sigma_p")),
           combo = c("dens_overlay", "trace"))


## Check Leave-One-Out Cross-Validation
fit_loo <- fit$loo()
fit_loo
plot(fit_loo, label_points = TRUE)

draws_df <- fit$draws(format = "draws_df")

## Individual estimates (posterior mean)
est_stan <- draws_df %>%
  spread_draws(VC[ID], VMAX[ID], KM[ID],
               eta_vc[ID], eta_vmax[ID], eta_km[ID]) %>% 
  mean_qi() %>% 
  select(ID, VC, VMAX, KM, 
         eta_vc, eta_vmax, eta_km)

post_preds_summary <- draws_df %>%
  spread_draws(pred[i], ipred[i], dv_ppc[i]) %>%
  mean_qi(pred, ipred, dv_ppc) %>%
  mutate(DV = nonmem_data$DV[nonmem_data$evid == 0][i],
         bloq = nonmem_data$bloq[nonmem_data$evid == 0][i])


ppc_ind <- post_preds_summary %>% 
  mutate(ID = nonmem_data$ID[nonmem_data$evid == 0][i],
         time = nonmem_data$time[nonmem_data$evid == 0][i])

p_dv_vs_pred <- ggplot(post_preds_summary %>% 
                         filter(bloq == 0), aes(x = pred, y = DV)) +
  geom_point() +
  theme_bw() +
  geom_abline(slope = 1, intercept = 0, color = "blue", size = 1.5) +
  xlab("Population Predictions") +
  ylab("Observed Concentration") +
  theme(axis.text = element_text(size = 14, face = "bold"),
        axis.title = element_text(size = 18, face = "bold"),
        axis.line = element_line(size = 2)) +
  scale_y_log10() +
  scale_x_log10()

p_dv_vs_ipred <- ggplot(post_preds_summary %>% 
                          filter(bloq == 0), aes(x = ipred, y = DV)) +
  geom_point() +
  theme_bw() +
  geom_abline(slope = 1, intercept = 0, color = "blue", size = 1.5) +
  xlab("Individual Predictions") +
  ylab("Observed Concentration") +
  theme(axis.text = element_text(size = 14, face = "bold"),
        axis.title = element_text(size = 18, face = "bold"),
        axis.line = element_line(size = 2)) +
  scale_y_log10() +
  scale_x_log10()

ggpubr::ggarrange(p_dv_vs_pred, p_dv_vs_ipred)

tmp <- ggplot(ppc_ind) +
  ggforce::facet_wrap_paginate(~ID, labeller = label_both,
                               nrow = 3, ncol = 4,
                               page = 1, scales = "free_y")

for(i in 1:ggforce::n_pages(tmp)){
  print(ggplot() +
          geom_ribbon(data = ppc_ind, 
                      mapping = aes(x = time, ymin = ipred.lower, 
                                    ymax = ipred.upper, group = ID),
                      fill = "blue", alpha = 0.5, show.legend = FALSE) +
          geom_ribbon(data = ppc_ind, 
                      mapping = aes(x = time, ymin = dv_ppc.lower, 
                                    ymax = dv_ppc.upper, group = ID),
                      fill = "blue", alpha = 0.25, show.legend = FALSE) +
          geom_line(data = ppc_ind, mapping = aes(x = time, y = ipred, 
                                                  group = ID),
                    linetype = 1, size = 1.15) +
          geom_line(data = ppc_ind, mapping = aes(x = time, y = pred, 
                                                  group = ID),
                    linetype = 2, size = 1.05) +
          geom_point(data = ppc_ind %>% 
                       filter(bloq == 0), 
                     mapping = aes(x = time, y = DV, group = ID),
                     size = 2, color = "red", show.legend = FALSE) +
          scale_y_continuous(name = "Drug Conc. (ug/mL)",
                             limits = c(NA, NA),
                             trans = "identity") +
          scale_x_continuous(name = "Time (d)",
                             breaks = seq(0, 168, by = 24),
                             labels = seq(0, 168/24, by = 1)) +
          theme_bw() +
          theme(axis.text = element_text(size = 14, face = "bold"),
                axis.title = element_text(size = 18, face = "bold"),
                legend.position = "bottom") +
          ggforce::facet_wrap_paginate(~ ID, labeller = label_both,
                                       nrow = 3, ncol = 4,
                                       page = i, scales = "free_y"))
  
}

## Standardized Random Effects (posterior mean)
eta_std <- est_stan %>% 
  select(ID, starts_with("eta")) %>% 
  pivot_longer(cols = starts_with("eta"), 
               names_to = "parameter", 
               values_to = "eta",
               names_prefix = "eta_") %>% 
  group_by(parameter) %>% 
  mutate(eta_std = (eta - mean(eta))/sd(eta)) %>% 
  ungroup() %>% 
  mutate(parameter = toupper(parameter))

## Standardized Random Effects (posterior mean) with Standard Normal overlayed
eta_std %>% 
  ggplot(aes(x = eta_std, group = parameter)) + 
  geom_histogram(aes(y = ..density..)) + 
  geom_density(color = "blue", size = 1.5) +
  stat_function(fun = dnorm, args = list(mean = 0, sd = 1),
                color = "red", size = 1) +
  scale_x_continuous(name = "Standarized Indiv. Effect",
                     limits = c(-2.5, 2.5)) +
  theme_bw(18) +
  facet_wrap(~parameter, scales = "free_x")

eta_std %>% 
  ggplot(aes(x = parameter, y = eta_std)) + 
  geom_boxplot() +
  scale_x_discrete(name = "Parameter") +
  scale_y_continuous(name = "Standarized Indiv. Effect") +
  theme_bw(18) +
  geom_hline(yintercept = 0, linetype = "dashed")

eta_std %>% 
  ggplot(aes(x = parameter, y = eta_std)) + 
  geom_boxplot(notch = TRUE) +
  scale_x_discrete(name = "Parameter") +
  scale_y_continuous(name = "Standarized Indiv. Effect") +
  theme_bw(18) +
  geom_hline(yintercept = 0, linetype = "dashed")

residuals <- draws_df %>%
  spread_draws(res[i], wres[i], ires[i], iwres[i], ipred[i]) %>% 
  mutate(time = nonmem_data$time[nonmem_data$evid == 0][i],
         bloq = nonmem_data$bloq[nonmem_data$evid == 0][i])

residuals %>% 
  filter(bloq == 0) %>%
  sample_draws(4) %>% 
  ggplot(aes(sample = iwres)) + 
  geom_qq() +
  geom_abline() +
  theme_bw() +
  facet_wrap(~.draw, labeller = label_both)


some_residuals_qq <- residuals %>% 
  filter(bloq == 0) %>%
  sample_draws(100) %>% 
  ggplot(aes(sample = iwres)) + 
  geom_qq() +
  geom_abline() +
  theme_bw() +
  transition_manual(.draw)

animate(some_residuals_qq, nframes = 100, width = 384, height = 384, res = 96, 
        dev = "png", type = "cairo")

residuals %>% 
  filter(bloq == 0) %>%
  sample_draws(4) %>%
  ggplot(aes(x = iwres)) + 
  geom_histogram(aes(y = ..density..)) +
  geom_density(color = "blue", size = 1.5) +
  stat_function(fun = dnorm, args = list(mean = 0, sd = 1),
                color = "red", size = 1) +
  theme_bw() +
  facet_wrap(~ .draw, labeller = label_both)


residuals %>% 
  filter(bloq == 0) %>%
  sample_draws(9) %>%
  mutate(qn_lower = qnorm(0.025),
         qn_upper = qnorm(0.975)) %>% 
  ggplot(aes(x = time, y = iwres)) + 
  geom_point() +
  geom_hline(yintercept = 0, color = "blue", size = 1.5) +
  geom_hline(aes(yintercept = qn_lower), linetype = "dashed", color = "blue",
             size = 1.25) +
  geom_hline(aes(yintercept = qn_upper), linetype = "dashed", color = "blue",
             size = 1.25) +
  theme_bw() +
  xlab("Time (h)") +
  facet_wrap(~ .draw, labeller = label_both)

residuals %>% 
  filter(bloq == 0) %>%
  sample_draws(9) %>%
  mutate(qn_lower = qnorm(0.025),
         qn_upper = qnorm(0.975)) %>% 
  ggplot(aes(x = ipred, y = iwres)) + 
  geom_point() +
  geom_smooth(se = FALSE, color = "blue", size = 1.5) +
  geom_hline(yintercept = 0, color = "red", size = 1.5) +
  geom_hline(aes(yintercept = qn_lower), linetype = "dashed", color = "red",
             size = 1.25) +
  geom_hline(aes(yintercept = qn_upper), linetype = "dashed", color = "red",
             size = 1.25) +
  theme_bw() +
  facet_wrap(~ .draw, labeller = label_both)

some_residuals_scatter <- residuals %>% 
  filter(bloq == 0) %>%
  sample_draws(30) %>% 
  mutate(qn_lower = qnorm(0.025),
         qn_upper = qnorm(0.975)) %>%
  ggplot(aes(x = time, y = iwres)) + 
  geom_point() +
  geom_smooth(se = FALSE, color = "blue", size = 1.5) +
  geom_hline(yintercept = 0, color = "red", size = 1.5) +
  geom_hline(aes(yintercept = qn_lower), linetype = "dashed", color = "red",
             size = 1.25) +
  geom_hline(aes(yintercept = qn_upper), linetype = "dashed", color = "red",
             size = 1.25) +
  theme_bw() +
  transition_manual(.draw)

animate(some_residuals_scatter, nframes = 30, width = 384, height = 384, 
        res = 96, dev = "png", type = "cairo", fps = 5)


# This one is hard to see when each individual has very similar timepoints
residuals %>% 
  filter(bloq == 0) %>%
  ggplot(aes(x = time, y = iwres, group = i)) +
  stat_pointinterval(.width = c(.95), point_size = 3) +
  theme_bw() +
  geom_hline(yintercept = 0, color = "blue", size = 1.5) +
  scale_x_continuous(name = "Time (h)",
                     limits = c(0.01, 24))

