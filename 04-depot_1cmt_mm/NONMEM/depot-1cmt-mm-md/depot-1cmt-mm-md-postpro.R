rm(list=ls())

library(tidyverse)
library(bayesplot)
library(kableExtra)
library(grid)
library(gridExtra)
library(rstan)

setwd("/data/BayesianBenchmarks")
color_scheme_set("brightblue")

modelName <- "depot-1cmt-mm"
modelDir <- paste0("/data/BayesianBenchmarks/", modelName, "/chains")
nChains <- 4
nIter <- 1000

## Get population parameters
out <- NULL
for(i in 1:nChains) {
params <- read.delim(paste0(modelDir, "/", modelName,"-", i, ".ext"), sep="", skip=1) %>% 
  select(ITERATION, V=THETA1, VMAX=THETA2, KM=THETA3, KA=THETA4,
         OMEGA_V=OMEGA.1.1., OMEGA_VMAX=OMEGA.2.2., OMEGA_KM=OMEGA.3.3., OMEGA_KA=OMEGA.4.4., 
         SIGMA=SIGMA.1.1.) %>% 
  filter(ITERATION>0) %>% 
  mutate(Chain=i)

out <- rbind(out, params)
}

plotData <- out %>% 
  select(-ITERATION) 

# Model Diagnostics
pdf(paste0("/data/BayesianBenchmarks/", modelName,"/", modelName, "-Diagnostics.pdf"), paper="a4r", height=0, width=0)

mcmc_trace(plotData, par=c("V","VMAX", "KM", "KA"))
mcmc_trace(plotData, par=c("OMEGA_V", "OMEGA_VMAX", "OMEGA_KM", "OMEGA_KA", "SIGMA"))

mcmc_hist_by_chain(plotData, par=c("V","VMAX", "KM", "KA"))
mcmc_hist_by_chain(plotData, par=c("OMEGA_V", "OMEGA_VMAX", "OMEGA_KM", "OMEGA_KA", "SIGMA"))

mcmc_dens_overlay(plotData, par=c("V","VMAX", "KM", "KA"))
mcmc_dens_overlay(plotData, par=c("OMEGA_V", "OMEGA_VMAX", "OMEGA_KM", "OMEGA_KA", "SIGMA"))

mcmc_areas(plotData, par=c("V"))
mcmc_areas(plotData, par=c("VMAX"))
mcmc_areas(plotData, par=c("KM"))
mcmc_areas(plotData, par=c("KA"))
mcmc_areas(plotData, par=c("OMEGA_V", "OMEGA_VMAX", "OMEGA_KM", "OMEGA_KA"))
mcmc_areas(plotData, par=c("SIGMA"))

mcmc_pairs(plotData, par=c("V","VMAX", "KM", "KA"))
mcmc_pairs(plotData, par=c("OMEGA_V", "OMEGA_VMAX", "OMEGA_KM", "OMEGA_KA", "SIGMA"))

mcmc_acf(plotData, par=c("V","VMAX", "KM", "KA"), lags=10)
mcmc_acf(plotData, par=c("OMEGA_V", "OMEGA_VMAX", "OMEGA_KM", "OMEGA_KA", "SIGMA"), lags=10)

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
  select(parameter, mean, SD, median, "95% CI", Neff, Rhat) # %>%
  # kbl(caption = "Summary of Posterior Parameter Estimates", row.names = F) %>%
  # kable_styling(bootstrap_options = "striped", full_width = F) 


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

