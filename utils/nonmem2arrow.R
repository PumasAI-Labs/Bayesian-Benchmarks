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

convert_nonmem_to_arrow <- function(nonmem_dir) {
  ext_files <- dir(nonmem_dir, ".ext", full.names = TRUE)
  lst_files <- dir(nonmem_dir, ".lst", full.names = TRUE)

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

  arrow_file <- file.path(nonmem_dir, "chains.arrow")
    
  # Write Arrow file
  arrow::write_feather(chains, arrow_file)
}
