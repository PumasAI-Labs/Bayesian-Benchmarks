source("utils/stan2arrow.R")
library(furrr)
library(stringr)

env = str_c("02-depot_1cmt_linear/Stan/Torsten/Fits/single_dose_", 1:5, ".rds")
env_multi = str_c("02-depot_1cmt_linear/Stan/Torsten/Fits/multiple_dose_", 1:5, ".rds")
env_mat_exp = str_c("02-depot_1cmt_linear/Stan/Torsten/Fits/single_dose_mat_exp_", 1:5, ".rds")
env_multi_mat_exp = str_c("02-depot_1cmt_linear/Stan/Torsten/Fits/multiple_dose_mat_exp_", 1:5, ".rds")

parameters = c(
  "TVCL",
  "TVVC",
  "TVKA",
  "omega",
  "sigma_p"
)

future_walk(
  c(
    env,
    env_multi,
    env_mat_exp,
    env_multi_mat_exp
  ),
  ~ convert_env_to_arrow(.x, parameters=parameters)
)
