## alpha beta calculator

library(tidyverse)
library(GGally)


n_sims <- 10000

tvcl <- 4
tvvc <- 70
tvq <- 4
tvvp <- 40
tvka <- 1

parameters <- tibble(cl = rlnorm(n_sims, log(tvcl), 0.4),
                     vc = rlnorm(n_sims, log(tvvc), 0.4),
                     q = rlnorm(n_sims, log(tvq), 0.4),
                     vp = rlnorm(n_sims, log(tvvp), 0.4),
                     ka = rlnorm(n_sims, log(tvka), 0.2)) %>% 
  mutate(ke = cl/vc,
         k_cp = q/vc,
         k_pc = q/vp,
         k_sum = k_cp + k_pc + ke,
         alpha = 0.5*(k_sum + sqrt(k_sum^2 - 4*k_pc*ke)),
         beta = 0.5*(k_sum - sqrt(k_sum^2 - 4*k_pc*ke)),
         flip_flop = ka < alpha,
         flip_flop_2 = ka < beta)

parameters$flip_flop %>% 
  mean

parameters %>% 
  group_by(flip_flop) %>% 
  summarize(across(everything(), mean))

ggpairs(parameters, columns = 1:5, 
        mapping = aes(color = flip_flop, alpha = 0.05)) +
  scale_y_log10() +
  scale_x_log10()
  
parameters %>% 
  select(ka, alpha) %>% 
  pivot_longer(ka:alpha) %>% 
  ggplot() +
  geom_density(aes(x = value, group = name, color = name))


parameters %>% 
  select(ka, alpha) %>% 
  ggplot() +
  geom_point(aes(x = alpha, y = ka, color = ka < alpha))

cl <- 4
vc <- 70
q <- 4
vp <- 50
ka <- 1

ke <- cl/vc
k_cp <- q/vc
k_pc <- q/vp

k_sum <- k_cp + k_pc + ke

beta <- 0.5*(k_sum - sqrt(k_sum^2 - 4*k_pc*ke))
alpha <- k_pc*ke/beta

round(c(alpha, beta), 4)
  
alpha <- 0.5*(k_sum + sqrt(k_sum^2 - 4*k_pc*ke))
beta <- 0.5*(k_sum - sqrt(k_sum^2 - 4*k_pc*ke))

round(c(alpha, beta), 4)

