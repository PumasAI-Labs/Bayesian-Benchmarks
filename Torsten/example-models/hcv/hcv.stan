functions{
  vector hcvODE(real t, vector x, array[] real parms, array[] real x_r, array[] int x_i) {
    // parameters
    real logKa = parms[1];
    real logKe = parms[2];
    real logVd = parms[3];
    real logn = parms[4];
    real logdelta = parms[5];
    real logc = parms[6];
    real logEC50 = parms[7];

    // constants
    real p = 100.0;
    real d = 0.001;
    real e = 1e-7;
    real s = 20000.0;

    vector[5] y;

    // custom inits
    real init_T = fmax(machine_precision(),
                       exp(logc + logdelta) / (p * e));
    real init_I = fmax(machine_precision(),
                       (s * e * p - d * exp(logc + logdelta)) / (p * exp(logdelta) * e));
    real init_W = fmax(machine_precision(),
                       (s * e * p - d * exp(logc + logdelta)) / (exp(logc + logdelta) * e));

    y[1] = -exp(logKa) * x[1];                                                      // X'
    y[2] = exp(logKa) * x[1] - exp(logKe) * x[2];                                   // A'
    y[3] = s - (x[3] + init_T) * (e * (x[5] + init_W) + d);                         // T'
    y[4] = e * (x[5] + init_W) * (x[3] + init_T) - exp(logdelta) * (x[4] + init_I); // I'
    y[5] =
        p / ((x[1] / exp(logVd) / exp(logEC50) + 1e-100)^exp(logn) + 1) *           // W'
          (x[4] + init_I) - exp(logc) * (x[5] + init_W);

    return y;
  }
}
data{
  int<lower=1> nt;
  int<lower=1> nObs;
  int<lower=1> nSubjects;
  array[nObs] int<lower=1> iObs;
  array[nSubjects] int<lower=1> start;
  array[nSubjects] int<lower=1> end;
  array[nt] int<lower=1> cmt;
  array[nt] int evid;
  array[nt] int addl;
  array[nt] int ss;
  array[nt] real amt;
  array[nt] real time;
  array[nt] real rate;
  array[nt] real ii;
  vector<lower=0>[nObs] yPK;
  vector<lower=0>[nObs] yPD;
}
parameters{
  // fixed effects
  real logthetaKa;
  real logthetaKe;
  real logthetaVd;
  real logthetan;
  real logthetadelta;
  real logthetac;
  real logthetaEC50;
  // random effects variance parameters, must be posisitive
  real<lower=0> omegaKa;
  real<lower=0> omegaKe;
  real<lower=0> omegaVd;
  real<lower=0> omegan;
  real<lower=0> omegadelta;
  real<lower=0> omegac;
  real<lower=0> omegaEC50;
  // variance parameter in proportional error model
  real<lower=0> sigmaPK;
  real<lower=0> sigmaPD;
}

transformed parameters{
  matrix<lower=0>[5, nt] x;
  row_vector<lower=0>[nt] conc;
  row_vector<lower=0>[nObs] concObs;
  row_vector<lower=0>[nt] log10W;
  row_vector<lower=0>[nObs] log10WObs;

  for(j in 1:nSubjects)
  {
    // empty vector
    array[7] real parms;
    // parameters
    parms[1] = logthetaKa + omegaKa; // Ka
    parms[2] = logthetaKe + omegaKe; // Ke
    parms[3] = logthetaVd + omegaVd; // Vd
    parms[4] = logthetan + omegan; // n
    parms[5] = logthetadelta + omegadelta; // delta
    parms[6] = logthetac + omegac; // c
    parms[7] = logthetaEC50 + omegaEC50; // EC50

    x = pmx_solve_bdf(hcvODE, 5, time,
                      amt, rate, ii, evid, cmt, addl, ss, parms,
                      1e-3, // rel_tol
                      1e-6, // abs_tol
                      1e5   // max_step
                      );

    // constants
    real p = 100.0;
    real d = 0.001;
    real e = 1e-7;
    real s = 20000.0;

    // inital conditions to W (x[5])
    real init_W = fmax(machine_precision(),
                       (s * e * p - d * exp(parms[6] + parms[5])) / (exp(parms[6] + parms[7]) * e));
    conc[start[j]:end[j]] = x[2, start[j]:end[j]] ./ exp(parms[3]);     // divide by exp(logVd)
    log10W[start[j]:end[j]] = log10(x[5, start[j]:end[j]] + init_W);    // log10(W)
  }

  concObs  = conc[iObs];
  log10WObs  = log10W[iObs];
}

model{
  // Prior
  logthetaKa ~ normal(log(0.8), 1);
  logthetaKe ~ normal(log(0.15), 1);
  logthetaVd ~ normal(log(100), 1);
  logthetan ~ normal(log(2.0), 1);
  logthetadelta ~ normal(log(0.20), 1);
  logthetac ~ normal(log(7.0), 1);
  logthetaEC50 ~ normal(log(0.12), 1);
  // random effects variance parameters, must be posisitive
  omegaKa ~ normal(0.25, 1);
  omegaKe ~ normal(0.25, 1);
  omegaVd ~ normal(0.25, 1);
  omegan ~ normal(0.25, 1);
  omegadelta ~ normal(0.25, 1);
  omegac ~ normal(0.25, 1);
  omegaEC50 ~ normal(0.25, 1);
  // variance parameter in proportional error model
  sigmaPK ~ normal(0.04, 1);
  sigmaPD ~ normal(0.04, 1);
  yPK ~ normal(concObs, sqrt(sigmaPK));
  yPD ~ normal(log10WObs, sqrt(sigmaPK));
}
