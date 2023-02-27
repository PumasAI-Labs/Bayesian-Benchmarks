// IV Infusion
// One-compartment PK Model with parallel linear and MM elimination
// IIV on CL, VC, VMAX, KM
// proportional plus additive error - DV = CP(1 + eps_p) + eps_a
// ODE solution using Torsten
// Observations are generated from a normal that is truncated below at 0
// Since we have a normal distribution on the error, but the DV must be > 0, it
//   generates values from a normal that is truncated below at 0

functions{
  
  real normal_lb_rng(real mu, real sigma, real lb){
    
    real p_lb = normal_cdf(lb | mu, sigma);
    real u = uniform_rng(p_lb, 1);
    real y = mu + sigma * inv_Phi(u);
    return y;

  }
  
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
  
  real<lower = 0> omega_cl;
  real<lower = 0> omega_vc;
  real<lower = 0> omega_vmax;
  real<lower = 0> omega_km;
  
  corr_matrix[4] R;  // Correlation matrix before transforming to Omega.
                     // Can in theory change this to having inputs for
                     // cor_cl_vc, ... and then construct the 
                     // correlation matrix in transformed data, but it's easy
                     // enough to do in R
  
  real<lower = 0> sigma_p;
  real<lower = 0> sigma_a;
  real<lower = -1, upper = 1> cor_p_a;
  
}
transformed data{
  
  int n_random = 4;
  int n_cmt = 2;

  vector[n_random] omega = [omega_cl, omega_vc, omega_vmax, omega_km]';
  
  matrix[n_random, n_random] L = cholesky_decompose(R);

  vector[2] sigma = [sigma_p, sigma_a]';
  matrix[2, 2] R_Sigma = rep_matrix(1, 2, 2);
  R_Sigma[1, 2] = cor_p_a;
  R_Sigma[2, 1] = cor_p_a;
  
  matrix[2, 2] Sigma = quad_form_diag(R_Sigma, sigma);
  
}
model{
  
}
generated quantities{
  
  vector[n_total] dv; // concentration with residual error
  
  vector[n_subjects] CL;
  vector[n_subjects] VC;
  vector[n_subjects] VMAX;
  vector[n_subjects] KM;
  
  {
  
    vector[n_random] typical_values = to_vector({TVCL, TVVC, TVVMAX, TVKM});
    
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
    
    for(j in 1:n_subjects){
      
      array[n_random + 1] real params = {CL[j], VC[j], VMAX[j], KM[j], 0};  // The 0 is for KA. Skip the absorption
                              
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
        real cp_tmp = cp[i];
        real sigma_tmp = sqrt(square(cp_tmp) * Sigma[1, 1] + Sigma[2, 2] + 
                              2*cp_tmp*Sigma[2, 1]);
        dv[i] = normal_lb_rng(cp_tmp, sigma_tmp, 0.0);
        
      }
    }
  }
}

