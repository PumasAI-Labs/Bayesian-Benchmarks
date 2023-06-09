rm(list=ls())

library(tidyverse)

# Functions
lowTriMat <- function(x) {
  diag(x)[upper.tri(diag(x), diag = T)]
}

projDir <- "/data/BayesianBenchmarks/clean"
modelName <- "depot-2cmt-linear-md"

setwd(projDir)
dir.create(paste0(modelName, "/chains"))

nChains <- 4
nTheta <- 5
nOmega <- 5


# Generate initial estimates
out <- list()
set.seed(11235)
for(i in 1:nChains) {
  out[[i]] <- list(signif(c(rnorm(1, log(4), 0.3),
                            rnorm(1, log(70), 0.3),
                            rnorm(1, log(4), 0.3),
                            rnorm(1, log(40), 0.3),
                            rnorm(1, log(1), 0.3),
                            rlnorm(5, log(0.3), 0.3),
                            rlnorm(1, log(0.2), 0.3)), 4))
}

omegaMat <- lowTriMat(out[i][[1]][[1]][(nTheta+1):(nTheta+nOmega)])
omegaMat[omegaMat==0] <- 0.01

# Create model file for each chain
for(i in 1:nChains) {
  modelText <- readLines(paste0(modelName, "/", modelName, "-template.mod"))
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
  
  write(modelText, paste0(modelName, "/chains/", modelName, "-",i,".mod"))
}

# Submit chains  in parallel
setwd(paste0(projDir, "/", modelName, "/chains"))
for(i in 1:nChains){
  system(paste0("qsub /data/bash/execute_8_nm75.sh ", modelName, "-", i, ".mod"))
}
