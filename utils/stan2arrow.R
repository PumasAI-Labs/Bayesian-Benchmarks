library(arrow)
library(readr)
library(cmdstanr)
library(dplyr)
library(magrittr)

convert_csv_to_arrow <- function(csv_files, arrow_dir, time_stamp) {
  # Create output directory if it doesn't exist
  if (!dir.exists(arrow_dir)) {
    dir.create(arrow_dir)
  }
  
  # Loop through CSV files
  for (csv_file in csv_files) {
    # Read CSV file
    df <- read_csv(csv_file, comment = "#")
    df$time <- time_stamp
    
    # Define Arrow file name and path
    arrow_file <- file.path(arrow_dir, paste0(basename(csv_file), ".arrow"))
    
    # Write Arrow file
    arrow::write_feather(df, arrow_file)
  }
}

convert_env_to_arrow <- function(env, parameters) {
  model_name <- basename(env)
  dir_name <- dirname(env)
  env <- read_rds(env)
  chains <- env$draws(c("lp__", parameters), format = "draws_df") %>% 
    as_tibble()
  times <- env$time()$chains
  chains %<>% left_join(
    times,
    by = c(".chain" = "chain_id")
  )
  
  arrow_file <- file.path(dir_name, gsub(".rds", ".arrow", model_name))
    
  # Write Arrow file
  arrow::write_feather(chains, arrow_file)
}
