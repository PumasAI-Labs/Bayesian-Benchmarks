## create initial estimates
init <- function(){
  list(TVCL = exp(rnorm(1, log(10), 0.2)),
       TVQ = exp(rnorm(1, log(20), 0.2)),
       TVVC = exp(rnorm(1, log(70), 0.2)),
       TVVP = exp(rnorm(1, log(70), 0.2)),
       TVKA = exp(rnorm(1, log(1), 0.2)),
       sigma = runif(1, 0.5, 2),
       L = diag(5),
       Z = matrix(rep(0, 5*10), nrow = 5),
       omega = runif(5, 0.01, 2))
}

inits <- init()
with(inits, 
     rstan::stan_rdump(ls(inits), 
                file = "Torsten/example-models/poppk2cpt/poppk2cpt.init_cd.R"))
