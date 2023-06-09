source("utils/stan2arrow.R")
library(furrr)
library(stringr)

env = str_c("04-depot_1cmt_mm/Stan/Torsten/Fits/single_dose_", 1:5, ".rds")
env_multi = str_c("04-depot_1cmt_mm/Stan/Torsten/Fits/multiple_dose_", 1:5, ".rds")

parameters = c(
  "TVVC",
  "TVVMAX",
  "TVKM",
  "TVKA",
  "omega",
  "sigma_p"
)

future_walk(
  c(
    env,
    env_multi,
  ),
  ~ convert_env_to_arrow(.x, parameters=parameters)
)
