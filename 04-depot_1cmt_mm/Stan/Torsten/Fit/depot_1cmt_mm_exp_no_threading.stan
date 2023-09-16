// First Order Absorption (oral/subcutaneous)
// One-compartment PK Model with MM elimination
// IIV on VC, VMAX, KM, KA
// exponential error - DV = CP*exp(eps)
// ODE solution using Torsten
// Deals with BLOQ values by the "CDF trick" (M3)

functions{
  
  vector depot_1cmt_mm_ode(real t, vector y, array[] real params, 
                           array[] real x_r, array[] int x_i){
  
    real vc = params[1];
    real vmax = params[2];
    real km = params[3];
    real ka = params[4];
    
    real conc = y[2]/vc;
    
    vector[2] dydt;

    dydt[1] = -ka*y[1];
    dydt[2] = ka*y[1] - vmax*conc/(km + conc);
    
    return dydt;
  }
  
}
data{
  
  int n_subjects;
  int n_total;
  int n_obs;
  array[n_obs] int i_obs;
  array[n_total] int ID;
  array[n_total] real amt;
  array[n_total] int cmt;
  array[n_total] int evid;
  array[n_total] real rate;
  array[n_total] real ii;
  array[n_total] int addl;
  array[n_total] int ss;
  array[n_total] real time;
  vector<lower = 0>[n_total] dv;
  array[n_subjects] int subj_start;
  array[n_subjects] int subj_end;
  vector[n_total] lloq;
  array[n_total] int bloq;
  
  // real<lower = 0> location_tvcl;    // Prior Location parameter for CL
  real<lower = 0> location_tvvc;    // Prior Location parameter for VC
  real<lower = 0> location_tvvmax;  // Prior Location parameter for VMAX
  real<lower = 0> location_tvkm;    // Prior Location parameter for KM
  real<lower = 0> location_tvka;    // Prior Location parameter for KA
  
  // real<lower = 0> scale_tvcl;       // Prior Scale parameter for CL
  real<lower = 0> scale_tvvc;       // Prior Scale parameter for VC
  real<lower = 0> scale_tvvmax;     // Prior Scale parameter for VMAX
  real<lower = 0> scale_tvkm;       // Prior Scale parameter for KM
  real<lower = 0> scale_tvka;       // Prior Scale parameter for KA
  
  // real<lower = 0> scale_omega_cl;   // Prior scale parameter for omega_cl
  real<lower = 0> scale_omega_vc;   // Prior scale parameter for omega_vc
  real<lower = 0> scale_omega_vmax; // Prior scale parameter for omega_vmax
  real<lower = 0> scale_omega_km;   // Prior scale parameter for omega_km
  real<lower = 0> scale_omega_ka;   // Prior scale parameter for omega_ka
  
  real<lower = 0> lkj_df_omega;     // Prior degrees of freedom for omega cor mat
  
  real<lower = 0> scale_sigma;      // Prior Scale parameter for exponential error
 
}
transformed data{ 
  
  int grainsize = 1;
  
  vector<lower = 0>[n_obs] dv_obs = dv[i_obs];
  array[n_obs] int dv_obs_id = ID[i_obs];
  
  vector[n_obs] lloq_obs = lloq[i_obs];
  array[n_obs] int bloq_obs = bloq[i_obs];
  
  int n_random = 4;                    // Number of random effects
  int n_cmt = 2;                       // Number of compartments in ODEs
  
  array[n_random] real scale_omega = {scale_omega_vc, scale_omega_vmax, 
                                      scale_omega_km, scale_omega_ka}; 
  
  array[n_cmt] real bioav = rep_array(1.0, n_cmt); // Hardcoding, but could be data or a parameter in another situation
  array[n_cmt] real tlag = rep_array(0.0, n_cmt);
  
}
parameters{ 
  
  // real<lower = 0> TVCL;       
  real<lower = 0> TVVC;
  real<lower = 0> TVVMAX;       
  real<lower = 0> TVKM;
  real<lower = 0> TVKA;
  
  vector<lower = 0>[n_random] omega;
  cholesky_factor_corr[n_random] L;
  
  real<lower = 0> sigma;
  
  matrix[n_random, n_subjects] Z;
  
}
transformed parameters{
  
  vector[n_subjects] VC;
  vector[n_subjects] VMAX;
  vector[n_subjects] KM;
  vector[n_subjects] KA;
  
  vector[n_obs] ipred;

  {

    row_vector[n_random] typical_values = 
                                      to_row_vector({TVVC, TVVMAX, TVKM, TVKA});
    
    matrix[n_subjects, n_random] eta = diag_pre_multiply(omega, L * Z)';

    matrix[n_subjects, n_random] theta =
                          (rep_matrix(typical_values, n_subjects) .* exp(eta));
                          
    vector[n_total] dv_ipred;
    matrix[n_total, n_cmt] x_ipred;                           
                          
    VC = col(theta, 1);
    VMAX = col(theta, 2);
    KM = col(theta, 3);
    KA = col(theta, 4);
    
    for(j in 1:n_subjects){
      
      x_ipred[subj_start[j]:subj_end[j],] =
        pmx_solve_rk45(depot_1cmt_mm_ode,
                       n_cmt,
                       time[subj_start[j]:subj_end[j]],
                       amt[subj_start[j]:subj_end[j]],
                       rate[subj_start[j]:subj_end[j]],
                       ii[subj_start[j]:subj_end[j]],
                       evid[subj_start[j]:subj_end[j]],
                       cmt[subj_start[j]:subj_end[j]],
                       addl[subj_start[j]:subj_end[j]],
                       ss[subj_start[j]:subj_end[j]],
                       {VC[j], VMAX[j], KM[j], KA[j]}, bioav, tlag)';
                      
      dv_ipred[subj_start[j]:subj_end[j]] = 
        x_ipred[subj_start[j]:subj_end[j], 2] ./ VC[j];
    
    }
    
    ipred = dv_ipred[i_obs];
  
  }
  
}
model{ 
  
  // Priors
  // TVCL ~ lognormal(log(location_tvcl), scale_tvcl);
  TVVC ~ lognormal(log(location_tvvc), scale_tvvc);
  TVVMAX ~ lognormal(log(location_tvvmax), scale_tvvmax);
  TVKM ~ lognormal(log(location_tvkm), scale_tvkm);
  TVKA ~ lognormal(log(location_tvka), scale_tvka);

  omega ~ normal(0, scale_omega);
  L ~ lkj_corr_cholesky(lkj_df_omega);
  
  sigma ~ normal(0, scale_sigma);
  
  to_vector(Z) ~ std_normal();
  
  // Likelihood
  for(i in 1:n_obs){
    if(bloq_obs[i] == 1){
      target += lognormal_lcdf(lloq_obs[i] | log(ipred[i]), sigma);
    }else{
      target += lognormal_lpdf(dv_obs[i] | log(ipred[i]), sigma);
    }
  }
}
