rm(list = ls())
cat("\014")

library(tidyjson)
library(tidyverse)

n_subjects <- read_csv("04-depot_1cmt_mm/data/single_dose.csv",
# nonmem_data <- read_csv("04-depot_1cmt_mm/data/multiple_dose.csv",
                       na = ".") %>% 
  distinct(ID) %>%
  count() %>%
  deframe()

write_inits <- function(run, chain){
  list(TVVC = rlnorm(1, log(70), 0.3),
       TVVMAX = rlnorm(1, log(1.000), 0.3),
       TVKM = rlnorm(1, log(0.250), 0.3),
       TVKA = rlnorm(1, log(1), 0.3),
       omega = rlnorm(4, log(0.3), 0.3),
       sigma_p = rlnorm(1, log(0.2), 0.3),
       L = diag(4),
       Z = matrix(rnorm(n_subjects*4), nrow = n_subjects, ncol = 4)) %>% 
    map(round, 3) %>% 
    write_stan_json(str_c("04-depot_1cmt_mm/data/inits/inits_", 
                          run, "_", chain, ".json"))
}

expand_grid(run = 1:10, chain = 1:4) %>% 
  {map2(.$run, .$chain, write_inits)}


