rm(list = ls())

source("utils/run_chains.R")
source("utils/write_nonmem_files_with_inits.R")

library(tidyjson)
library(tidyverse)

# Functions

modelName <- "depot-1cmt-mm-md"
root_folder <- "04-depot_1cmt_mm"

base_name <- paste0(root_folder, "/NONMEM/", modelName, "/chains/", modelName)
dir_name <- paste0(root_folder, "/NONMEM/", modelName)
dir.create(paste0(dir_name, "/chains"))

grid_of_runs_and_chains <- expand_grid(run = 1:5, chain = 1:4)

pwalk(grid_of_runs_and_chains, create_model_file)

# Submit chains in parallel
run_chains(root_folder, modelName, nchains = 1, threads_per_chain = 1, run = 1)
run_chains(root_folder, modelName, nchains = 4, threads_per_chain = 8, run = 2)
run_chains(root_folder, modelName, nchains = 4, threads_per_chain = 8, run = 3)
run_chains(root_folder, modelName, nchains = 4, threads_per_chain = 8, run = 4)
run_chains(root_folder, modelName, nchains = 4, threads_per_chain = 8, run = 5)
run_chains(root_folder, modelName, nchains = 4, threads_per_chain = 8, run = 6)
run_chains(root_folder, modelName, nchains = 4, threads_per_chain = 8, run = 7)
run_chains(root_folder, modelName, nchains = 4, threads_per_chain = 8, run = 8)
run_chains(root_folder, modelName, nchains = 4, threads_per_chain = 8, run = 9)
run_chains(root_folder, modelName, nchains = 4, threads_per_chain = 8, run = 10)

clean_nonmem_dir(root_folder)
