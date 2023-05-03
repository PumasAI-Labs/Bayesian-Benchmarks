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

modelName <- "friberg"
thetas <- c("CL", "V1", "Q", "V2", "KA", "MTT", "CIRC0", "GAMMA", "ALPHA")
randoms <- c("OMEGA_CL", "OMEGA_V1", "OMEGA_Q", "OMEGA_V2", "OMEGA_KA", 
             "OMEGA_MTT", "OMEGA_CIRC0", "OMEGA_GAMMA", "OMEGA_ALPHA", "SIGMA")

modelDir <- paste0(modelName, "/chains")
nChains <- 4
nIter <- 400

## Get population parameters
out <- NULL
for(i in 1:nChains) {
  params <- read.delim(paste0(modelDir, "/", modelName,"-", i, ".ext"), sep="", skip=4012) %>% 
    select(ITERATION, CL=THETA1, V1=THETA2, Q=THETA3, V2=THETA4, KA=THETA5,  
           MTT=THETA6, CIRC0=THETA7, GAMMA=THETA8, ALPHA=THETA9,
           OMEGA_CL=OMEGA.1.1., OMEGA_V1=OMEGA.2.2., OMEGA_Q=OMEGA.3.3., OMEGA_V2=OMEGA.4.4., OMEGA_KA=OMEGA.5.5., 
           OMEGA_MTT=OMEGA.6.6., OMEGA_CIRC0=OMEGA.7.7., OMEGA_GAMMA=OMEGA.8.8., OMEGA_ALPHA=OMEGA.9.9., SIGMA=SIGMA.1.1.) %>% 
    filter(ITERATION>0) %>% 
    mutate(Chain=i)
  
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

parTable <- monitor(parArray, warmup = 200, print = FALSE) %>% 
  as.matrix() %>%
  formatC(3) %>%
  as.data.frame


t1 <- parTable %>% 
  rename(SEmean = se_mean, SD = sd, pct2.5 = "2.5%", pct25 = "25%", median = "50%",
         pct75 = "75%", pct97.5 = "97.5%", Neff = "n_eff") %>%
  mutate(parameter = rownames(.), 
         "95% CI" = paste("(", pct2.5, ", ", pct97.5, ")", sep = "")) %>%
  select(parameter, mean, SD, median, "95% CI", Neff, Rhat) 

grid.arrange(arrangeGrob(tableGrob(t1, rows=NULL), 
                         top = textGrob("Summary of Posterior Parameter Estimates\n",
                                        gp=gpar(fontsize=20)), nrow=2, heights = c(3, 1)))

## Time for estimation
outTime <- data.frame(Model=modelName, Chain=1:nChains, Time=NA)
outTimeAll <- NULL
Est <- c("BAYES", "NUTS")
for(j in 1:length(Est)) {
  for(i in 1:nChains){
    lstText <- readLines(paste0(modelDir, "/", modelName, "-", i, ".lst"))
    estLine <- grep("Elapsed estimation  time in seconds: ", lstText)
    estTime <- as.numeric(str_extract(lstText[estLine[j]], "[[:digit:]]+\\.*[[:digit:]]*"))
    outTime$Time[outTime$Chain==i] <- estTime
    outTime$Est <- Est[j]
  }
  outTimeAll <- rbind(outTimeAll, outTime)
}

t2 <- outTimeAll %>% 
  group_by(Chain) %>% 
  summarize(TotalTime=sum(Time)) %>%
  ungroup() %>% 
  mutate(Avg=signif(mean(TotalTime), 4))

grid.arrange(arrangeGrob(tableGrob(t2, rows=NULL), 
                         top = textGrob("Summary of Estimation Time\n",
                                        gp=gpar(fontsize=20)), nrow=2))
dev.off()

