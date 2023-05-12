rm(list = ls())
cat("\014")

library(tidyjson)
library(tidyverse)

n_subjects <- read_csv("05-friberg/data/multiple_dose.csv",
                       na = ".") %>% 
  distinct(ID) %>%
  count() %>%
  deframe()

write_inits <- function(run, chain){
  list(TVCL = rlnorm(1, log(4), 0.3),
       TVVC = rlnorm(1, log(70), 0.3),
       TVQ = rlnorm(1, log(4), 0.3),
       TVVP = rlnorm(1, log(40), 0.3),
       TVKA = rlnorm(1, log(1), 0.3),
       TVMTT = rlnorm(1, log(125), 0.3),
       TVCIRC0 = rlnorm(1, log(5), 0.3),
       TVGAMMA = rlnorm(1, log(0.17), 0.3),
       TVALPHA = rlnorm(1, log(3e-4), 0.3),
       omega = rlnorm(9, log(0.3), 0.3),
       sigma_p = rlnorm(1, log(0.2), 0.3),
       sigma_p_pd = rlnorm(1, log(0.2), 0.3),
       L = diag(9),
       Z = matrix(0, nrow = n_subjects, ncol = 9)) %>% 
    map(round, 5) %>% 
    write_stan_json(str_c("05-friberg/data/inits/inits_", 
                          run, "_", chain, ".json"))
}

expand_grid(run = 1:10, chain = 1:4) %>% 
  {map2(.$run, .$chain, write_inits)}


