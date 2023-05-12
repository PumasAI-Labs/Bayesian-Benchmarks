rm(list = ls())
cat("\014")

library(tidyjson)
library(tidyverse)

n_subjects <- read_csv("03-depot_2cmt_linear/data/single_dose.csv",
# nonmem_data <- read_csv("03-depot_2cmt_linear/data/multiple_dose.csv",
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
       omega = rlnorm(5, log(0.3), 0.3),
       sigma_p = rlnorm(1, log(0.2), 0.3),
       L = diag(5),
       Z = matrix(rnorm(n_subjects*5), nrow = n_subjects, ncol = 5)) %>% 
    map(round, 3) %>% 
    write_stan_json(str_c("03-depot_2cmt_linear/data/inits/inits_", 
                          run, "_", chain, ".json"))
}

expand_grid(run = 1:10, chain = 1:4) %>% 
  {map2(.$run, .$chain, write_inits)}


