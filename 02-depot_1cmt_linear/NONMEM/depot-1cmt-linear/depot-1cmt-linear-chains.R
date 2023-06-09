rm(list=ls())
source("utils/run_chains.R")

library(tidyverse)

# Functions
lowTriMat <- function(x) {
  diag(x)[upper.tri(diag(x), diag = T)]
}

modelName <- "depot-1cmt-linear"
root_folder <- "02-depot_1cmt_linear"

nChains <- 4
nTheta <- 3
nOmega <- 3


# Generate initial estimates
out <- list()
set.seed(11235)
for(i in 1:nChains) {
  out[[i]] <- list(signif(c(rnorm(1, log(4), 0.3),
                            rnorm(1, log(70), 0.3),
                            rnorm(1, log(1), 0.3),
                            rlnorm(3, log(0.3), 0.3),
                            rlnorm(1, log(0.2), 0.3)), 4))
}

omegaMat <- lowTriMat(out[i][[1]][[1]][(nTheta+1):(nTheta+nOmega)])
omegaMat[omegaMat==0] <- 0.01

base_name <- paste0(root_folder, "/NONMEM/", modelName, "/chains/", modelName)
dir_name <- paste0(root_folder, "/NONMEM/", modelName)
dir.create(paste0(dir_name, "/chains"))

# Create model file for each chain
for(i in 1:nChains) {
  modelText <- readLines(paste0(dir_name, "/", modelName, "-template.mod"))
  modelTextLine <- grep("THETAUPDATE", modelText)
  modelText[modelTextLine] <- paste(out[i][[1]][[1]][1:nTheta], collapse=" ")
  
  modelTextLine <- grep("OMEGAUPDATE", modelText)
  modelText[modelTextLine] <- paste(omegaMat, collapse=" ")
  
  modelTextLine <- grep("SIGMAUPDATE", modelText)
  modelText[modelTextLine] <- paste(out[i][[1]][[1]][(nTheta+nOmega+1):length(out[i][[1]][[1]])], collapse=" ")
  
  modelTextLine <- grep("ITERFILE", modelText)
  modelText[modelTextLine] <- gsub("ITERFILE",  paste0(modelName, "-iter-",i,".tab"), modelText[modelTextLine])
  
  modelTextLine <- grep("SDFILE", modelText)
  modelText[modelTextLine] <- gsub("SDFILE",  paste0(modelName, "-sdtab-",i,".tab"), modelText[modelTextLine])
  
  modelTextLine <- grep("PAFILE", modelText)
  modelText[modelTextLine] <- gsub("PAFILE",  paste0(modelName, "-patab-",i,".tab"), modelText[modelTextLine])
  
  write(modelText, paste0(base_name, "-",i,".mod"))
}

# Submit chains  in parallel
run_chains(root_folder, modelName, nchains=4, threads_per_chain=2)
