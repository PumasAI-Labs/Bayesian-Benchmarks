// IV Infusion
// Two-compartment PK Model with MM elimination
// IIV on VC, Q, VP, VMAX, KM
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
  
  vector two_cmt_lin_plus_mm_ode(real t, vector y, array[] real params, 
                                 array[] real x_r, array[] int x_i){
    
    real cl = params[1];
    real vc = params[2];
    real q = params[3];
    real vp = params[4];
    real vmax = params[5];
    real km = params[6];
    real ka = params[7];
    
    real ke = cl/vc;
    real k_cp = q/vc;
    real k_pc = q/vp;
    real conc = y[2]/vc;
    
    vector[3] dydt;

    dydt[1] = -ka*y[1];
    dydt[2] = ka*y[1] - (ke + k_cp)*y[2] - vmax*conc/(km + conc) + k_pc*y[3];
    dydt[3] = k_cp*y[2] - k_pc*y[3];
    
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
  
  int n_random = 5;       // Number of random effects
  int n_cmt = 3;          // Number of compartments (depot, central, peripheral)

}
parameters{ 
  
  // real<lower = 0> TVCL;       
  real<lower = 0> TVVC;
  real<lower = 0> TVQ;
  real<lower = 0> TVVP;
  real<lower = 0> TVVMAX;       
  real<lower = 0> TVKM;
  
  vector<lower = 0>[n_random] omega;
  cholesky_factor_corr[n_random] L;
  
  real<lower = 0> sigma_p;
  
  matrix[n_random, n_subjects] Z;
  
}
generated quantities{

  vector[n_time_new] cp; // concentration with no residual error
  vector[n_time_new] dv; // concentration with residual error

  {
    row_vector[n_random] typical_values = to_row_vector({TVVC, TVQ, TVVP,
                                                        TVVMAX, TVKM});
    
    matrix[n_subjects_new, n_random] eta_new;
    matrix[n_subjects_new, n_random] theta_new;
    matrix[n_time_new, n_cmt] x_cp;

    for(i in 1:n_subjects_new){
      eta_new[i, ] = multi_normal_cholesky_rng(rep_vector(0, n_random),
                                               diag_pre_multiply(omega, L))';
    }
    theta_new = (rep_matrix(typical_values, n_subjects_new) .* exp(eta_new));

    for(j in 1:n_subjects_new){
      
      array[n_random + 2] real theta_params = 
                    to_array_1d(append_col(0, append_col(theta_new[j], 0))); // access the parameters for subject nn
      
      x_cp[subj_start[j]:subj_end[j],] =
        pmx_solve_rk45(two_cmt_lin_plus_mm_ode,
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

      cp[subj_start[j]:subj_end[j]] = x_cp[subj_start[j]:subj_end[j], 2] ./ theta_params[2];
    
    }


    for(i in 1:n_time_new){
      if(cp[i] == 0){
         dv[i] = 0;
      }else{
        real cp_tmp = cp[i];
        real sigma_tmp = cp_tmp*sigma_p;
        dv[i] = normal_lb_rng(cp_tmp, sigma_tmp, 0.0);
      }
    }
  
  }

}


