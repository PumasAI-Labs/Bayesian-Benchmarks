// First Order Absorption (oral/subcutaneous)
// One-compartment PK Model with MM elimination
// IIV on VC, VMAX, KM, KA
// proportional error - DV = CP(1 + eps_p)
// ODE solution using Torsten
// Predictions are generated from a normal that is truncated below at 0

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
  int n_time_new;
  array[n_time_new] real time;
  array[n_time_new] real amt;
  array[n_time_new] int cmt;
  array[n_time_new] int evid;
  array[n_time_new] real rate;
  array[n_time_new] real ii;
  array[n_time_new] int addl;
  array[n_time_new] int ss;
  array[n_subjects] int subj_start;
  array[n_subjects] int subj_end;
 
}
transformed data{ 
  
  int n_random = 4;                    // Number of random effects
  int n_cmt = 2;                       // Number of compartments in ODEs

}
parameters{ 
  
  // real<lower = 0> TVCL;       
  real<lower = 0> TVVC;
  real<lower = 0> TVVMAX;       
  real<lower = 0> TVKM;
  real<lower = 0> TVKA;
  
  vector<lower = 0>[n_random] omega;
  cholesky_factor_corr[n_random] L;
  
  real<lower = 0> sigma_p;
  
  matrix[n_random, n_subjects] Z;
  
}
generated quantities{
  
  vector[n_time_new] ipred; // ipred for the observed individuals at the new timepoints
  vector[n_time_new] pred;  // pred for the observed individuals at the new timepoints
  vector[n_time_new] dv;    // dv for the observed individuals at the new timepoints
 
  {
    row_vector[n_random] typical_values = 
                                      to_row_vector({TVVC, TVVMAX, TVKM, TVKA});

    matrix[n_random, n_random] R = multiply_lower_tri_self_transpose(L);
    matrix[n_random, n_random] Omega = quad_form_diag(R, omega);

    matrix[n_subjects, n_random] eta = diag_pre_multiply(omega, L * Z)';

    matrix[n_subjects, n_random] theta =
                          (rep_matrix(typical_values, n_subjects) .* exp(eta));
                          
    matrix[n_time_new, n_cmt] x_pred;
    matrix[n_time_new, n_cmt] x_ipred;
    
    array[n_random + 1] real theta_params_tv = {0, TVVC, TVVMAX, TVKM, TVKA};
    
    vector[n_subjects] VC = col(theta, 1);
    vector[n_subjects] VMAX = col(theta, 2);
    vector[n_subjects] KM = col(theta, 3);
    vector[n_subjects] KA = col(theta, 4);
    
    for(j in 1:n_subjects){

      array[n_random + 1] real theta_params = {0, VC[j], VMAX[j], KM[j], KA[j]};

      x_ipred[subj_start[j]:subj_end[j],] =
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
                      theta_params)';



     ipred[subj_start[j]:subj_end[j]] =
       x_ipred[subj_start[j]:subj_end[j], 2] ./ VC[j];

     x_pred[subj_start[j]:subj_end[j],] =
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
                      theta_params_tv)';

      pred[subj_start[j]:subj_end[j]] =
        x_pred[subj_start[j]:subj_end[j], 2] ./ TVVC;
    }
    
    for(i in 1:n_time_new){
      if(ipred[i] == 0){
        dv[i] = 0;
      }else{
        real ipred_tmp = ipred[i];
        real sigma_tmp = ipred_tmp*sigma_p;
        dv[i] = normal_lb_rng(ipred_tmp, sigma_tmp, 0.0);
      }
    }
  }
}


  