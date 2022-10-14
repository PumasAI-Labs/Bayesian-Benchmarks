library(cmdstanr)

# Setting CmdStan path to Torsten intall
set_cmdstan_path("Torsten/cmdstan")

# Compile Stan model
m <- cmdstan_model("Torsten/example-models/pk2cpt_ode/pk2cpt_ode.stan")

# data
addl <- c(14, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
amt <- c(80000, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
cmt <- c(1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2)
cObs <- c(359.239725273613, 662.647763577673, 1106.23947626728,
1185.25862586227, 1802.42630338229, 2296.47896937277, 2008.04412122616,
2000.94001020581, 1115.28848387488, 902.769414322048, 445.989900210258,
285.638744806589, 333.509854464029, 664.522440159846, 960.044660816181,
1157.24318831826, 1593.58511079002, 2170.22768767046, 2175.56839105518,
2168.6215871781)
evid <- c(1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
ii <- c(12, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
iObs <- c(2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21)
nObs <- 20
nt <- 21
rate <- c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
ss <- c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
time <- c(0, 0.083, 0.167, 0.25, 0.5, 0.75, 1, 1.5, 2, 3, 4, 6, 8, 12,
12.083, 12.167, 12.25, 12.5, 12.75, 13, 13.5)

stan_data <- list(addl=addl, amt=amt, cmt=cmt, cObs=cObs, evid=evid, ii=ii,
                  iObs=iObs, nObs=nObs, nt=nt, rate=rate, ss=ss, time=time)

# inits
CL <- 7.4367958406427
ka <- 1.0811298754049
Q <- 28.0799996152587
sigma <- 0.589695154260051
V1 <- 78.4460632446725
V2 <- 68.1255965629187

stan_inits <- function() list(CL=CL, ka=ka, Q=Q, sigma=sigma, V1=V1, V2=V2)

# Stan fit
fit_pk2cpt_ode <- m$sample(
    data=stan_data,
    init=stan_inits,
    parallel_chains=4,
    adapt_delta=0.8,
    iter_sampling=1000,
    iter_warmup=1000
)
