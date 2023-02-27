// IV Infusion
// One-compartment PK Model with MM elimination
// IIV on VC, VMAX, KM
// proportional plus additive error - DV = CP(1 + eps_p) + eps_a
// ODE solution using Torsten
// Implements threading for within-chain parallelization 
// Deals with BLOQ values by the "CDF trick" (M4)
// Since we have a normal distribution on the error, but the DV must be > 0, it
//   truncates the likelihood below at 0
// For PPC, it generates values from a normal that is truncated below at 0


functions{

  array[] int sequence(int start, int end) { 
    array[end - start + 1] int seq;
    for (n in 1:num_elements(seq)) {
      seq[n] = n + start - 1;
    }
    return seq; 
  } 
  
  int num_between(int lb, int ub, array[] int y){
    
    int n = 0;
    for(i in 1:num_elements(y)){
      if(y[i] >= lb && y[i] <= ub)
         n = n + 1;
    }
    return n;
    
  }
  
  array[] int find_between(int lb, int ub, array[] int y) {
    // vector[num_between(lb, ub, y)] result;
    array[num_between(lb, ub, y)] int result;
    int n = 1;
    for (i in 1:num_elements(y)) {
      if (y[i] >= lb && y[i] <= ub) {
        result[n] = y[i];
        n = n + 1;
      }
    }
    return result;
  }
  
  vector find_between_vec(int lb, int ub, array[] int idx, vector y) {
    
    vector[num_between(lb, ub, idx)] result;
    int n = 1;
    if(num_elements(idx) != num_elements(y)) reject("illegal input");
    for (i in 1:rows(y)) {
      if (idx[i] >= lb && idx[i] <= ub) {
        result[n] = y[i];
        n = n + 1;
      }
    }
    return result;
  }
  
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
  
  vector iv_1cmt_mm_ode(real t, vector y, array[] real params, 
                        array[] real x_r, array[] int x_i){
    
    real vc = params[1];
    real vmax = params[2];
    real km = params[3];
    
    real conc = y[1]/vc;
    
    vector[1] dydt;

    dydt[1] = -vmax*conc/(km + conc);
    
    return dydt;
  }
  
  real partial_sum_lpmf(array[] int seq_subj, int start, int end,
                        vector dv_obs, array[] int dv_obs_id, array[] int i_obs,
                        array[] real amt, array[] int cmt, array[] int evid, 
                        array[] real time, array[] real rate, array[] real ii, 
                        array[] int addl, array[] int ss,
                        array[] int subj_start, array[] int subj_end, 
                        real TVVC, real TVVMAX, real TVKM,
                        vector omega, matrix L, matrix Z, 
                        real sigma_p, 
                        vector lloq, array[] int bloq,
                        int n_random, int n_subjects, int n_total, int n_cmt){
                           
    real ptarget = 0;
    row_vector[n_random] typical_values = 
                                      to_row_vector({TVVC, TVVMAX, TVKM});
    
    matrix[n_subjects, n_random] eta = diag_pre_multiply(omega, L * Z)';

    matrix[n_subjects, n_random] theta =
                          (rep_matrix(typical_values, n_subjects) .* exp(eta));
    
                              
    int N = end - start + 1;    // number of subjects in this slice  
    vector[n_total] dv_ipred;   
    matrix[n_total, 1] x_ipred;
  
    int n_obs_slice = num_between(subj_start[start], subj_end[end], i_obs);
    array[n_obs_slice] int i_obs_slice = find_between(subj_start[start], 
                                                      subj_end[end], i_obs);
                                                
    vector[n_obs_slice] dv_obs_slice = find_between_vec(start, end, 
                                                        dv_obs_id, dv_obs);
    
    vector[n_obs_slice] ipred_slice;
    
    vector[n_obs_slice] lloq_slice = lloq[i_obs_slice];
    array[n_obs_slice] int bloq_slice = bloq[i_obs_slice];
    
    
    for(n in 1:N){            // loop over subjects in this slice
    
      int nn = n + start - 1; // nn is the ID of the current subject
      
      // row_vector[n_random] theta_nn = theta[nn]; // access the parameters for subject nn
      // real vc = theta_nn[1];
      // real vmax = theta_nn[2];
      // real km = theta_nn[3];
      // 
      // array[n_random] real theta_params = {vc, vmax, km}; 
      
      array[n_random] real theta_params = to_array_1d(theta[nn]); // access the parameters for subject nn
      
      x_ipred[subj_start[nn]:subj_end[nn],] =
        pmx_solve_rk45(iv_1cmt_mm_ode,
                       n_cmt,
                       time[subj_start[nn]:subj_end[nn]],
                       amt[subj_start[nn]:subj_end[nn]],
                       rate[subj_start[nn]:subj_end[nn]],
                       ii[subj_start[nn]:subj_end[nn]],
                       evid[subj_start[nn]:subj_end[nn]],
                       cmt[subj_start[nn]:subj_end[nn]],
                       addl[subj_start[nn]:subj_end[nn]],
                       ss[subj_start[nn]:subj_end[nn]],
                       theta_params)';
                      
      dv_ipred[subj_start[nn]:subj_end[nn]] = 
        x_ipred[subj_start[nn]:subj_end[nn], 1] ./ theta_params[1];
    
    }
  
    ipred_slice = dv_ipred[i_obs_slice];
    
    for(i in 1:n_obs_slice){
      real sigma_tmp = ipred_slice[i]*sigma_p;
      if(bloq_slice[i] == 1){
        ptarget += log_diff_exp(normal_lcdf(lloq_slice[i] | ipred_slice[i], 
                                                            sigma_tmp),
                                normal_lcdf(0.0 | ipred_slice[i], sigma_tmp)) -
                   normal_lccdf(0.0 | ipred_slice[i], sigma_tmp); 
      }else{
        ptarget += normal_lpdf(dv_obs_slice[i] | ipred_slice[i], sigma_tmp) -
                   normal_lccdf(0.0 | ipred_slice[i], sigma_tmp);
      }
    }                                         
                              
    return ptarget;
                           
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
  
  // real<lower = 0> scale_tvcl;       // Prior Scale parameter for CL
  real<lower = 0> scale_tvvc;       // Prior Scale parameter for VC
  real<lower = 0> scale_tvvmax;     // Prior Scale parameter for VMAX
  real<lower = 0> scale_tvkm;       // Prior Scale parameter for KM
  
  // real<lower = 0> scale_omega_cl;   // Prior scale parameter for omega_cl
  real<lower = 0> scale_omega_vc;   // Prior scale parameter for omega_vc
  real<lower = 0> scale_omega_vmax; // Prior scale parameter for omega_vmax
  real<lower = 0> scale_omega_km;   // Prior scale parameter for omega_km
  
  real<lower = 0> lkj_df_omega;     // Prior degrees of freedom for omega cor mat
  
  real<lower = 0> scale_sigma_p;    // Prior Scale parameter for proportional error
 
}
transformed data{ 
  
  int grainsize = 1;
  
  vector<lower = 0>[n_obs] dv_obs = dv[i_obs];
  array[n_obs] int dv_obs_id = ID[i_obs];
  
  vector[n_obs] lloq_obs = lloq[i_obs];
  array[n_obs] int bloq_obs = bloq[i_obs];
  
  int n_random = 3;                    // Number of random effects
  int n_cmt = 1;                       // Number of compartments in ODEs
  
  array[n_random] real scale_omega = {scale_omega_vc,
                                      scale_omega_vmax, scale_omega_km}; 
  
  array[n_subjects] int seq_subj = sequence(1, n_subjects); // reduce_sum over subjects
  
}
parameters{ 
  
  // real<lower = 0> TVCL;       
  real<lower = 0> TVVC;
  real<lower = 0> TVVMAX;       
  real<lower = 0> TVKM;
  
  vector<lower = 0>[n_random] omega;
  cholesky_factor_corr[n_random] L;
  
  real<lower = 0> sigma_p;
  
  matrix[n_random, n_subjects] Z;
  
}

model{ 
  
  // Priors
  // TVCL ~ lognormal(log(location_tvcl), scale_tvcl);
  TVVC ~ lognormal(log(location_tvvc), scale_tvvc);
  TVVMAX ~ lognormal(log(location_tvvmax), scale_tvvmax);
  TVKM ~ lognormal(log(location_tvkm), scale_tvkm);

  omega ~ normal(0, scale_omega);
  L ~ lkj_corr_cholesky(lkj_df_omega);
  
  sigma_p ~ normal(0, scale_sigma_p);
  
  to_vector(Z) ~ std_normal();
  
  // Likelihood
  target += reduce_sum(partial_sum_lupmf, seq_subj, grainsize,
                       dv_obs, dv_obs_id, i_obs,
                       amt, cmt, evid, time, 
                       rate, ii, addl, ss, subj_start, subj_end, 
                       TVVC, TVVMAX, TVKM, omega, L, Z,
                       sigma_p,
                       lloq, bloq,
                       n_random, n_subjects, n_total, n_cmt);
}
generated quantities{
  
  // real<lower = 0> sigma_sq_p = square(sigma_p);

  // real<lower = 0> omega_cl = omega[1];
  real<lower = 0> omega_vc = omega[1];
  real<lower = 0> omega_vmax = omega[2];
  real<lower = 0> omega_km = omega[3];

  // // real<lower = 0> omega_sq_cl = square(omega_cl);
  // real<lower = 0> omega_sq_vc = square(omega_vc);
  // real<lower = 0> omega_sq_vmax = square(omega_vmax);
  // real<lower = 0> omega_sq_km = square(omega_km);
  
  // // real cor_cl_vc;
  // // real cor_cl_vmax;
  // // real cor_cl_km;
  // real cor_vc_vmax;
  // real cor_vc_km;
  // real cor_vmax_km;
  
  // // real omega_cl_vc;
  // // real omega_cl_vmax;
  // // real omega_cl_km;
  // real omega_vc_vmax;
  // real omega_vc_km;
  // real omega_vmax_km;
  
  // // vector[n_subjects] eta_cl;
  // vector[n_subjects] eta_vc;
  // vector[n_subjects] eta_vmax;
  // vector[n_subjects] eta_km;
  // vector[n_subjects] CL;
  vector[n_subjects] VC;
  vector[n_subjects] VMAX;
  vector[n_subjects] KM;

  // vector[n_obs] ipred;
  // vector[n_obs] pred;
  // vector[n_obs] dv_ppc;
  // vector[n_obs] log_lik;
  // vector[n_obs] res;
  // vector[n_obs] wres;
  // vector[n_obs] ires;
  // vector[n_obs] iwres;
 
  {
    row_vector[n_random] typical_values = 
                                      to_row_vector({TVVC, TVVMAX, TVKM});

    matrix[n_random, n_random] R = multiply_lower_tri_self_transpose(L);
    matrix[n_random, n_random] Omega = quad_form_diag(R, omega);

    matrix[n_subjects, n_random] eta = diag_pre_multiply(omega, L * Z)';

    matrix[n_subjects, n_random] theta =
                          (rep_matrix(typical_values, n_subjects) .* exp(eta));

    // vector[n_total] dv_pred;
    // matrix[n_total, n_cmt] x_pred;
    // vector[n_total] dv_ipred;
    // matrix[n_total, n_cmt] x_ipred;
    // 
    // array[n_random] real theta_params_tv = {TVVC, TVVMAX, TVKM};
    // 
    // // eta_cl = col(eta, 1);
    // eta_vc = col(eta, 1);
    // eta_vmax = col(eta, 2);
    // eta_km = col(eta, 3);

    // CL = col(theta, 1);
    VC = col(theta, 1);
    VMAX = col(theta, 2);
    KM = col(theta, 3);

    // // cor_cl_vc = R[1, 2];
    // // cor_cl_vmax = R[1, 3];
    // // cor_cl_km = R[1, 4];
    // cor_vc_vmax = R[1, 1];
    // cor_vc_km = R[1, 3];
    // cor_vmax_km = R[2, 3];
    // 
    // // omega_cl_vc = Omega[1, 2];
    // // omega_cl_vmax = Omega[1, 3];
    // // omega_cl_km = Omega[1, 4];
    // omega_vc_vmax = Omega[1, 1];
    // omega_vc_km = Omega[1, 3];
    // omega_vmax_km = Omega[2, 3];

    // for(j in 1:n_subjects){
    // 
    //   array[n_random] real theta_params = {VC[j], VMAX[j], KM[j]};
    //   
    
    
    // x_ipred[subj_start[j]:subj_end[j],] =
    //   pmx_solve_rk45(iv_1cmt_mm_ode,
    //                  n_cmt,
    //                  time[subj_start[j]:subj_end[j]],
    //                  amt[subj_start[j]:subj_end[j]],
    //                  rate[subj_start[j]:subj_end[j]],
    //                  ii[subj_start[j]:subj_end[j]],
    //                  evid[subj_start[j]:subj_end[j]],
    //                  cmt[subj_start[j]:subj_end[j]],
    //                  addl[subj_start[j]:subj_end[j]],
    //                  ss[subj_start[j]:subj_end[j]],
    //                  theta_params)';

    //         
    //                   
    //  dv_ipred[subj_start[j]:subj_end[j]] = 
    //    x_ipred[subj_start[j]:subj_end[j], 1] ./ VC[j];
    // 
    // x_pred[subj_start[j]:subj_end[j],] =
    //   pmx_solve_rk45(iv_1cmt_mm_ode,
    //                  n_cmt,
    //                  time[subj_start[j]:subj_end[j]],
    //                  amt[subj_start[j]:subj_end[j]],
    //                  rate[subj_start[j]:subj_end[j]],
    //                  ii[subj_start[j]:subj_end[j]],
    //                  evid[subj_start[j]:subj_end[j]],
    //                  cmt[subj_start[j]:subj_end[j]],
    //                  addl[subj_start[j]:subj_end[j]],
    //                  ss[subj_start[j]:subj_end[j]],
    //                  theta_params_tv)';
    // 
    //   dv_pred[subj_start[j]:subj_end[j]] = 
    //     x_pred[subj_start[j]:subj_end[j], 1] ./ TVVC;
    // }

    // pred = dv_pred[i_obs];
    // ipred = dv_ipred[i_obs];

  }

  // res = dv_obs - pred;
  // ires = dv_obs - ipred;
  // 
  // for(i in 1:n_obs){
  //   real ipred_tmp = ipred[i];
  //   real sigma_tmp = ipred_tmp*sigma_p;
  //   dv_ppc[i] = normal_lb_rng(ipred_tmp, sigma_p, 0.0);
  //   if(bloq_obs[i] == 1){
  //     // log_lik[i] = log(normal_cdf(lloq_obs[i] | ipred_tmp, sigma_tmp) -
  //     //                  normal_cdf(0.0 | ipred_tmp, sigma_tmp)) -
  //     //              normal_lccdf(0.0 | ipred_tmp, sigma_tmp);
  //     log_lik[i] = log_diff_exp(normal_lcdf(lloq_obs[i] | ipred_tmp, sigma_tmp),
  //                               normal_lcdf(0.0 | ipred_tmp, sigma_tmp)) -
  //                  normal_lccdf(0.0 | ipred_tmp, sigma_tmp);
  //   }else{
  //     log_lik[i] = normal_lpdf(dv_obs[i] | ipred_tmp, sigma_tmp) -
  //                  normal_lccdf(0.0 | ipred_tmp, sigma_tmp);
  //   }
  //   wres[i] = res[i]/sigma_tmp;
  //   iwres[i] = ires[i]/sigma_tmp;
  // }
  
}

