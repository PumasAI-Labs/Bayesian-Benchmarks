source("utils/stan2arrow.R")
library(furrr)
library(stringr)

env = str_c("01-iv_2cmt_linear/Stan/Torsten/Fits/single_dose_", 1:5, ".rds")
env_multi = str_c("01-iv_2cmt_linear/Stan/Torsten/Fits/multiple_dose_", 1:5, ".rds")


parameters = c(
  "TVCL",
  "TVVC",
  "TVQ",
  "TVVP",
  "omega",
  "sigma"
)

future_walk(
  c(
    env,
    env_multi,
  ),
  ~ convert_env_to_arrow(.x, parameters=parameters)
)
