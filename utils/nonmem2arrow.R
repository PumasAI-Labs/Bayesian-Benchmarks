library(arrow)
library(readr)
library(dplyr)
library(magrittr)
library(janitor)
library(purrr)
library(stringr)

extract_time <- function(lst_file) {
  lines <- readLines(lst_file)
  line <- grep("Elapsed estimation  time in seconds: ", lines)
  time <- as.numeric(str_extract(lines[line], "[[:digit:]]+\\.*[[:digit:]]*"))
  return(time)
}

# TODO: fix the new run thing
convert_nonmem_to_arrow <- function(nonmem_dir, run) {
  ext_files <- dir(nonmem_dir, pattern = paste0("run", run, "-\\d{1}", "\\.ext"), full.names = TRUE)
  lst_files <- dir(nonmem_dir, pattern = paste0("run", run, "-\\d{1}", "\\.lst"), full.names = TRUE)

  chains <- map_dfr(
    ext_files,
    ~ read_table(.x, skip = 1) %>%  clean_names() %>% filter(iteration > 0),
    .id = "chain"
  ) %>% 
    mutate(
      chain = as.numeric(chain)
    ) %>% 
    select(-mcmcobj)
  
  times <- map_dbl(
    lst_files,
    extract_time
  )
  
  times <- tibble(
    chain = seq_along(times),
    time = times
  )
  
  chains %<>% left_join(times, by = "chain")

  arrow_file <- file.path(nonmem_dir, paste0("chains-run", run, ".arrow"))
    
  # Write Arrow file
  arrow::write_feather(chains, arrow_file)
}
