// First Order Absorption (oral/subcutaneous)
// One-compartment PK Model
// IIV on CL, VC, KA
// exponential error - DV = CP*exp(eps)
// Analytical solution using Torsten

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
  real<lower = 0> TVKA;
  
  real<lower = 0> omega_cl;
  real<lower = 0> omega_vc;
  real<lower = 0> omega_ka;
  
  corr_matrix[3] R;  // Correlation matrix before transforming to Omega.
                     // Can in theory change this to having inputs for
                     // cor_cl_vc, cor_cl_kja, ... and then construct the 
                     // correlation matrix in transformed data, but it's easy
                     // enough to do in R
  
  real<lower = 0> sigma;
  
}
transformed data{
  
  int n_random = 3;

  vector[n_random] omega = [omega_cl, omega_vc, omega_ka]';
  
  matrix[n_random, n_random] L = cholesky_decompose(R);
  
}
model{
  
}
generated quantities{
  
  vector[n_total] dv; // concentration with residual error
  
  vector[n_subjects] CL;
  vector[n_subjects] VC;
  vector[n_subjects] KA;
  
  {
  
    vector[n_random] typical_values = to_vector({TVCL, TVVC, TVKA});
    
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
    KA = col(theta, 3);
    
    for(j in 1:n_subjects){
      
      array[n_random] real theta_params = {CL[j], VC[j], KA[j]};
      
      x_cp[subj_start[j]:subj_end[j],] =
        pmx_solve_onecpt(time[subj_start[j]:subj_end[j]],
                         amt[subj_start[j]:subj_end[j]],
                         rate[subj_start[j]:subj_end[j]],
                         ii[subj_start[j]:subj_end[j]],
                         evid[subj_start[j]:subj_end[j]],
                         cmt[subj_start[j]:subj_end[j]],
                         addl[subj_start[j]:subj_end[j]],
                         ss[subj_start[j]:subj_end[j]],
                         theta_params)';

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

