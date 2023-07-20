// First Order Absorption (oral/subcutaneous)
// One-compartment PK Model with MM elimination
// IIV on VC, VMAX, KM, KA
// exponential error - DV = CP*exp(eps)
// ODE solution using Torsten

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
  int n_subjects_new;
  int n_time_new;
  array[n_time_new] real time;
  array[n_time_new] real amt;
  array[n_time_new] int cmt;
  array[n_time_new] int evid;
  array[n_time_new] real rate;
  array[n_time_new] real ii;
  array[n_time_new] int addl;
  array[n_time_new] int ss;
  array[n_subjects_new] int subj_start;
  array[n_subjects_new] int subj_end;
  
}
transformed data{ 
  
  int n_random = 4;       // Number of random effects
  int n_cmt = 2;          // Number of compartments (depot, central)

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
generated quantities{

  vector[n_time_new] cp; // concentration with no residual error
  vector[n_time_new] dv; // concentration with residual error

  {
    row_vector[n_random] typical_values = 
                                      to_row_vector({TVVC, TVVMAX, TVKM, TVKA});
    
    matrix[n_subjects_new, n_random] eta_new;
    matrix[n_subjects_new, n_random] theta_new;
    matrix[n_time_new, n_cmt] x_cp;

    for(i in 1:n_subjects_new){
      eta_new[i, ] = multi_normal_cholesky_rng(rep_vector(0, n_random),
                                               diag_pre_multiply(omega, L))';
    }
    theta_new = (rep_matrix(typical_values, n_subjects_new) .* exp(eta_new));

    for(j in 1:n_subjects_new){
      
      array[n_random] real theta_params = to_array_1d(theta_new[j]); // access the parameters for subject nn
      
      x_cp[subj_start[j]:subj_end[j],] =
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
                       theta_params)';

      cp[subj_start[j]:subj_end[j]] = x_cp[subj_start[j]:subj_end[j], 2] ./ theta_params[1];
    
    }


    for(i in 1:n_time_new){
      if(cp[i] == 0){
        dv[i] = 0;
      }else{
        dv[i] = lognormal_rng(log(cp[i]), sigma);
      }
    }
  
  }

}


