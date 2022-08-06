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


# Stan fit
fit_poppk2cpt <- m$sample(
    data=stan_data,
    # init=stan_inits,
    parallel_chains=4,
    adapt_delta=0.8,
    iter_sampling=1000,
    iter_warmup=1000
)

results <- fit_poppk2cpt$summary()
