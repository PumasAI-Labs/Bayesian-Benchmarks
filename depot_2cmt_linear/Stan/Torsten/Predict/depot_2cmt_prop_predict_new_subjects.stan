// First Order Absorption (oral/subcutaneous)
// Two-compartment PK Model
// IIV on CL, VC, Q, VP, and Ka (full covariance matrix)
// proportional error - DV = CP(1 + eps_p)
// Analytical solution using Torsten
// Predictions are generated from a normal that is truncated below at 0

functions{
  
  real normal_lb_rng(real mu, real sigma, real lb){
    
    real p_lb = normal_cdf(lb | mu, sigma);
    real u = uniform_rng(p_lb, 1);
    real y = mu + sigma * inv_Phi(u);
    return y;

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
  
  real<lower = 0> TVCL;       
  real<lower = 0> TVVC; 
  real<lower = 0> TVQ;       
  real<lower = 0> TVVP; 
  // real<lower = 0> TVKA;
  real<lower = 0.5*(TVCL/TVVC + TVQ/TVVC + TVQ/TVVP + 
    sqrt((TVCL/TVVC + TVQ/TVVC + TVQ/TVVP)^2 - 4*TVCL/TVVC*TVQ/TVVP))> TVKA; 
  
  vector<lower = 0>[n_random] omega;
  cholesky_factor_corr[n_random] L;
  
  real<lower = 0> sigma_p;
  
  matrix[n_random, n_subjects] Z;
  
}
generated quantities{

  vector[n_time_new] cp; // concentration with no residual error
  vector[n_time_new] dv; // concentration with residual error

  {
    row_vector[n_random] typical_values = 
      to_row_vector({TVCL, TVVC, TVQ, TVVP, TVKA});
    
    matrix[n_subjects_new, n_random] eta_new;
    matrix[n_subjects_new, n_random] theta_new;
    matrix[n_time_new, n_cmt] x_cp;

    for(i in 1:n_subjects_new){
      eta_new[i, ] = multi_normal_cholesky_rng(rep_vector(0, n_random),
                                               diag_pre_multiply(omega, L))';
    }
    theta_new = (rep_matrix(typical_values, n_subjects_new) .* exp(eta_new));

    for(j in 1:n_subjects_new){
      
      row_vector[n_random] theta_j = theta_new[j]; // access the parameters for subject j
      real cl = theta_j[1];
      real vc = theta_j[2];
      real q = theta_j[3];
      real vp = theta_j[4];
      real ka = theta_j[5];

      array[n_random] real theta_params = {cl, q, vc, vp, ka};
      
      x_cp[subj_start[j]:subj_end[j],] =
        pmx_solve_twocpt(time[subj_start[j]:subj_end[j]],
                         amt[subj_start[j]:subj_end[j]],
                         rate[subj_start[j]:subj_end[j]],
                         ii[subj_start[j]:subj_end[j]],
                         evid[subj_start[j]:subj_end[j]],
                         cmt[subj_start[j]:subj_end[j]],
                         addl[subj_start[j]:subj_end[j]],
                         ss[subj_start[j]:subj_end[j]],
                         theta_params)';

      cp[subj_start[j]:subj_end[j]] = x_cp[subj_start[j]:subj_end[j], 2] ./ vc;
    
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


