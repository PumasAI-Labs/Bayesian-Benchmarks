rm(list=ls())

library(tidyverse)
library(bayesplot)
library(kableExtra)
library(grid)
library(gridExtra)
library(rstan)
color_scheme_set("brightblue")

projDir <- "/data/BayesianBenchmarks/clean"
setwd(projDir)

modelName <- "depot-2cmt-linear"
thetas <- c("CL","V2", "Q", "V3", "KA")
randoms <- c("OMEGA_CL", "OMEGA_V2", "OMEGA_Q", "OMEGA_V3", "OMEGA_KA", "SIGMA")

modelDir <- paste0(modelName, "/chains")
nChains <- 4
nIter <- 1000

## Get population parameters
out <- NULL
for(i in 1:nChains) {
params <- read.delim(paste0(modelDir, "/", modelName,"-", i, ".ext"), sep="", skip=1) %>% 
  select(ITERATION, CL=THETA1, V2=THETA2, Q=THETA3, V3=THETA4, KA=THETA5,
         OMEGA_CL=OMEGA.1.1., OMEGA_V2=OMEGA.2.2., OMEGA_Q=OMEGA.3.3., OMEGA_V3=OMEGA.4.4., OMEGA_KA=OMEGA.5.5., 
         SIGMA=SIGMA.1.1.) %>% 
  filter(ITERATION>0) %>% 
  mutate(Chain=i) %>% 
  mutate(across(CL:V2, exp)) %>% 
  mutate(across(OMEGA_CL:SIGMA, sqrt))

  out <- rbind(out, params)
}

plotData <- out %>% 
  select(-ITERATION) 

# Model Diagnostics
pdf(paste0(modelName,"/", modelName, "-Diagnostics.pdf"), paper="a4r", height=0, width=0)

grid.arrange(textGrob(paste0("Diagnostics Report\n", modelName), gp=gpar(fontsize=20)), nrow=1)

mcmc_combo(plotData, par=thetas, combo = c("dens_overlay", "trace", "hist"))
mcmc_combo(plotData, par=randoms, combo = c("dens_overlay", "trace", "hist"))

mcmc_pairs(plotData, par=thetas, diag_fun = "dens")
mcmc_pairs(plotData, par=randoms, diag_fun = "dens")

mcmc_acf(plotData, par=thetas)
mcmc_acf(plotData, par=randoms)

# Posterior Parameter Tables
parArray <- array(double(nIter * nChains * (ncol(plotData) - 1)), 
                  dim = c(nIter, nChains, ncol(plotData) - 1), 
                  dimnames =  list(NULL, NULL, setdiff(names(plotData), c("Chain"))))

for(i in 1:nChains){
  parArray[,i,] <- plotData %>%
    filter(Chain == i) %>%
    select(-Chain) %>%
    as.matrix
}

parTable <- monitor(parArray, warmup = 500, print = FALSE) %>% 
  as.matrix() %>%
  formatC(3) %>%
  as.data.frame


t1 <- parTable %>% 
  rename(SEmean = se_mean, SD = sd, pct2.5 = "2.5%", pct25 = "25%", median = "50%",
         pct75 = "75%", pct97.5 = "97.5%", Neff = "n_eff") %>%
  mutate(parameter = rownames(.), 
         "95% CI" = paste("(", pct2.5, ", ", pct97.5, ")", sep = "")) %>%
  select(parameter, mean, SD, median, "95% CI", Neff, Rhat, Bulk_ESS, Tail_ESS) 

grid.arrange(arrangeGrob(tableGrob(t1, rows=NULL), 
                         top = textGrob("Summary of Posterior Parameter Estimates\n",
                                        gp=gpar(fontsize=20)), nrow=2))

## Time for estimation
outTime <- data.frame(Model=modelName, Chain=1:nChains, Time=NA)
for(i in 1:nChains){
  lstText <- readLines(paste0(modelDir, "/", modelName, "-", i, ".lst"))
  estLine <- grep("Elapsed estimation  time in seconds: ", lstText)
  estTime <- as.numeric(str_extract(lstText[estLine], "[[:digit:]]+\\.*[[:digit:]]*"))
  outTime$Time[outTime$Chain==i] <- estTime
}

t2 <- outTime %>% 
  mutate(Avg = signif(mean(Time), 4)) %>% 
  select(-Model)

grid.arrange(arrangeGrob(tableGrob(t2, rows=NULL), 
                         top = textGrob("Summary of Estimation Time\n",
                                        gp=gpar(fontsize=20)), nrow=2))

dev.off()

