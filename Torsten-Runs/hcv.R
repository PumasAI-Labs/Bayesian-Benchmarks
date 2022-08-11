library(cmdstanr)
library(dplyr)
library(readr)

# Setting CmdStan path to Torsten intall
set_cmdstan_path("Torsten/cmdstan")

# Compile Stan model
m <- cmdstan_model("Torsten/example-models/hcv/hcv.stan")

# data
df <- read_csv("data/hcv.csv")
# all are vector of length nt
nt <- length(df$time) # 45
time <- df$time
cmt <- df$cmt %>% replace(is.na(.), 0) + 1
amt <- df$amt %>% replace(is.na(.), 0)
rate <- df$rate %>% replace(is.na(.), 0)
evid <- df$evid
addl <- rep(0, nt)
ss <- df$ss
ii <- df$ii %>% replace(is.na(.), 0)
# EVIDs 33 0's and 12 1's
nObs <- df %>% filter(evid == 0) %>% nrow
iObs <- df %>% mutate(iObs = row_number()) %>% filter(evid == 0) %>% pull(iObs)
# 3 subjects
nSubjects <- df %>% pull(id) %>% unique %>% length
step <- 15
start <- c(1, 1+step, 1+step+step)
end <- c(15, 15+step, 15+step+step)
yPK <- df$yPK %>% na.omit
yPD <- df$yPD %>% na.omit

stan_data <- list(
  nt = nt,
  nObs = nObs,
  nSubjects = nSubjects,
  iObs = iObs,
  start = start,
  end = end,
  cmt = cmt,
  evid = evid,
  addl = addl,
  ss = ss,
  amt = amt,
  ii = ii,
  time = time,
  rate = rate,
  yPK = yPK,
  yPD = yPD
)

stan_inits <- function() {
  list(
    logthetaKa = -0.2231435513142097,
    logthetaKe = -1.8971199848858813,
    logthetaVd = 4.605170185988092,
    logthetan = 0.6931471805599453,
    logthetaδ = -1.6094379124341003,
    logthetac = 1.9459101490553132,
    logthetaEC50 = -2.120263536200091,
    omegaKa = 0.776261792561888,
    omegaKe = 0.776261792561888,
    omegaVd = 0.776261792561888,
    omegan = 0.776261792561888,
    omegaδ = 0.776261792561888,
    omegac = 0.776261792561888,
    omegaEC50 = 0.776261792561888,
    sigmaPK = 0.6895956605642456,
    sigmaPD = 0.6895956605642456)
}

# Stan fit
fit_poppk2cpt <- m$sample(
    data = stan_data,
    init = stan_inits,
    parallel_chains = 4,
    adapt_delta = 0.8,
    iter_sampling = 1000,
    iter_warmup = 1000
)

results <- fit_poppk2cpt$summary()
# All 4 chains finished successfully.
# Mean chain execution time: 5938.1 seconds.
# Total execution time: 8135.3 seconds.

# Warning: 1444 of 4000 (36.0%) transitions ended with a divergence.
# See https://mc-stan.org/misc/warnings for details.

# Warning: 2556 of 4000 (64.0%) transitions hit the maximum treedepth limit of 10.
# See https://mc-stan.org/misc/warnings for details.
results %>% write_csv("Torsten-Runs/hcv_results.csv")
