source("utils/stan2arrow.R")
library(furrr)

env = "03-depot_2cmt_linear/Stan/Torsten/Fits/single_dose.rds"
env_multi = "03-depot_2cmt_linear/Stan/Torsten/Fits/multiple_dose.rds"
env_mat_exp = "03-depot_2cmt_linear/Stan/Torsten/Fits/single_dose_mat_exp.rds"
env_multi_mat_exp = "03-depot_2cmt_linear/Stan/Torsten/Fits/multiple_dose_mat_exp.rds"

parameters = c(
  "TVCL",
  "TVVC",
  "TVQ",
  "TVVP",
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
