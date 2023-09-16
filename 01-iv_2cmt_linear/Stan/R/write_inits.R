rm(list = ls())
cat("\014")

library(tidyjson)
library(tidyverse)

n_subjects <- read_csv("01-iv_2cmt_linear/data/single_dose.csv",
# nonmem_data <- read_csv("01-iv_2cmt_linear/data/multiple_dose.csv",
                        na = ".") %>% 
  distinct(ID) %>%
  count() %>%
  deframe()


write_inits <- function(run, chain){
  list(TVCL = rlnorm(1, log(4), 0.3),
       TVVC = rlnorm(1, log(70), 0.3),
       TVQ = rlnorm(1, log(4), 0.3),
       TVVP = rlnorm(1, log(50), 0.3),
       omega = rlnorm(4, log(0.3), 0.3),
       sigma = rlnorm(1, log(0.2), 0.3),
       L = diag(4),
       Z = matrix(rnorm(n_subjects*4), nrow = n_subjects, ncol = 4)) %>% 
    map(round, 3) %>% 
    write_stan_json(str_c("01-iv_2cmt_linear/data/inits/inits_", 
                          run, "_", chain, ".json"))
}

expand_grid(run = 1:10, chain = 1:4) %>% 
  {map2(.$run, .$chain, write_inits)}


