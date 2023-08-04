library(furrr)
library(future)
plan(multisession)

run_chains <- function(folder, nonmem_model, nchains=4, threads_per_chain=8,
                       n_runs) {
  has_md_folder <- check_folder(nonmem_model)
  
  # nonmem_dir = if sd nonmem_model/NONMEM/nonmem_model/chains/
  #              if md nonmem_model-md/NONMEM/nonmem_model/chains/
  # Construct the filepath based on the condition
  if(has_md_folder) {
    nonmem_dir <- paste0(folder, "/NONMEM/", nonmem_model, "/chains/", nonmem_model, "-md")
  } else {
    nonmem_dir <- paste0(folder, "/NONMEM/", nonmem_model, "/chains/", nonmem_model)
  }

  for(run in 1:n_runs){
    future_map(
      1:nchains,
      ~ system(
        paste0("execute ",
               paste0(nonmem_dir, "-run", run, "-", .x, ".mod"),
               " -parafile=/opt/NONMEM/nm751/run/mpilinux8.pnm",
               " -nodes=",
               threads_per_chain),
        wait = FALSE
      )
    )
  }
}

check_folder <- function(filepath) {
  # Check if "-md" is part of any directory in the path
  has_md_folder <- grepl("/-md/", filepath)
  
  return(has_md_folder)
}


