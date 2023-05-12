rm(list = ls())
cat("\014")

library(tidyjson)
library(tidyverse)

n_subjects <- read_csv("02-depot_1cmt_linear/data/single_dose.csv",
# nonmem_data <- read_csv("02-depot_1cmt_linear/data/multiple_dose.csv",
                       na = ".") %>% 
  distinct(ID) %>%
  count() %>%
  deframe()

write_inits <- function(run, chain){
  list(TVCL = rlnorm(1, log(4), 0.3),
       TVVC = rlnorm(1, log(70), 0.3),
       TVKA = rlnorm(1, log(1), 0.3),
       omega = rlnorm(3, log(0.3), 0.3),
       sigma_p = rlnorm(1, log(0.2), 0.3),
       L = diag(3),
       Z = matrix(rnorm(n_subjects*3), nrow = n_subjects, ncol = 3)) %>% 
    map(round, 3) %>% 
    write_stan_json(str_c("02-depot_1cmt_linear/data/inits/inits_", 
                          run, "_", chain, ".json"))
}

expand_grid(run = 1:10, chain = 1:4) %>% 
  {map2(.$run, .$chain, write_inits)}


