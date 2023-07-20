source("utils/stan2arrow.R")
library(furrr)
library(stringr)

env_multi = str_c("05-friberg/Stan/Torsten/Fits/multiple_dose_", 1:5, ".rds")
env_multi_coupled = str_c("05-friberg/Stan/Torsten/Fits/multiple_dose_coupled_", 1:5, ".rds")

parameters = c(
  "TVCL",
  "TVVC",
  "TVQ",
  "TVVP",
  "TVKA",
  "TVMTT",
  "TVCIRC0",
  "TVGAMMA",
  "TVALPHA",
  "omega",
  "sigma",
  "sigma_pd"
)

future_walk(
  c(
    env_multi,
    env_multi_coupled
  ),
  ~ convert_env_to_arrow(.x, parameters=parameters)
)
