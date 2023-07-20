// First Order Absorption (oral/subcutaneous)
// One-compartment PK Model with parallel linear and MM elimination
// IIV on CL, VC, VMAX, KM, KA
// exponential error - DV = CP*exp(eps)
// ODE solution using Torsten

functions{
  
  vector one_cmt_lin_plus_mm_ode(real t, vector y, array[] real params, 
                                 array[] real x_r, array[] int x_i){
    
    real cl = params[1];
    real vc = params[2];
    real vmax = params[3];
    real km = params[4];
    real ka = params[5];
    
    real ke = cl/vc;
    real conc = y[2]/vc;
    
    vector[2] dydt;

    dydt[1] = -ka*y[1];
    dydt[2] = ka*y[1] - ke*y[2] - vmax*conc/(km + conc);
    
    return dydt;
  }
  
}
data{
  
  int n_subjects;
  int n_total;                  
  
  array[n_total] real amt;
  array[n_total] int cmt;
  array[n_total] int evid;
  array[n_total] real rate;
  array[n_total] real ii;
  array[n_total] int addl;
  array[n_total] int ss;
  array[n_total] real time;
  
  array[n_subjects] int subj_start;
  array[n_subjects] int subj_end;
  
  real<lower = 0> TVCL;
  real<lower = 0> TVVC;
  real<lower = 0> TVVMAX;
  real<lower = 0> TVKM;
  real<lower = 0> TVKA;
  
  real<lower = 0> omega_cl;
  real<lower = 0> omega_vc;
  real<lower = 0> omega_vmax;
  real<lower = 0> omega_km;
  real<lower = 0> omega_ka;
  
  corr_matrix[5] R;  // Correlation matrix before transforming to Omega.
                     // Can in theory change this to having inputs for
                     // cor_cl_vc, ... and then construct the 
                     // correlation matrix in transformed data, but it's easy
                     // enough to do in R
  
  real<lower = 0> sigma;
  
}
transformed data{
  
  int n_random = 5;
  int n_cmt = 2;

  vector[n_random] omega = [omega_cl, omega_vc, omega_vmax, omega_km, omega_ka]';
  
  matrix[n_random, n_random] L = cholesky_decompose(R);
  
}
model{
  
}
generated quantities{
  
  vector[n_total] dv; // concentration with residual error
  
  vector[n_subjects] CL;
  vector[n_subjects] VC;
  vector[n_subjects] VMAX;
  vector[n_subjects] KM;
  vector[n_subjects] KA;
  
  {
  
    vector[n_random] typical_values = 
                                    to_vector({TVCL, TVVC, TVVMAX, TVKM, TVKA});
    
    matrix[n_random, n_subjects] eta;   
    matrix[n_subjects, n_random] theta; 
    vector[n_total] cp; // concentration with no residual error
    matrix[n_total, 2] x_cp;
    
    for(i in 1:n_subjects){
      eta[, i] = multi_normal_cholesky_rng(rep_vector(0, n_random),
                                           diag_pre_multiply(omega, L));
    }
    theta = (rep_matrix(typical_values, n_subjects) .* exp(eta))';

    CL = col(theta, 1);
    VC = col(theta, 2);
    VMAX = col(theta, 3);
    KM = col(theta, 4);
    KA = col(theta, 5);
    
    for(j in 1:n_subjects){
      
      array[n_random] real params = {CL[j], VC[j], VMAX[j], KM[j], KA[j]};  // The 0 is for KA. Skip the absorption
                              
      x_cp[subj_start[j]:subj_end[j],] =
        pmx_solve_rk45(one_cmt_lin_plus_mm_ode,
                       n_cmt,
                       time[subj_start[j]:subj_end[j]],
                       amt[subj_start[j]:subj_end[j]],
                       rate[subj_start[j]:subj_end[j]],
                       ii[subj_start[j]:subj_end[j]],
                       evid[subj_start[j]:subj_end[j]],
                       cmt[subj_start[j]:subj_end[j]],
                       addl[subj_start[j]:subj_end[j]],
                       ss[subj_start[j]:subj_end[j]],
                       params)';
      
      cp[subj_start[j]:subj_end[j]] = x_cp[subj_start[j]:subj_end[j], 2] ./ VC[j];
    
    }

    for(i in 1:n_total){
      if(cp[i] == 0){
         dv[i] = 0;
      }else{
        dv[i] = lognormal_rng(log(cp[i]), sigma);
      }
    }
  }
}

