// Depot dosing with first-order absorption
// Two-compartment PK Model
// Macro parameterization (CL, VC, Q, VP, KA)
// IIV on CL, VC, Q, VP, KA (full covariance matrix)
// Exponential error - DV = CP*exp(eps)
// Closed form solution using Torsten
// Implements threading for within-chain parallelization
// Deals with BLOQ values by the "CDF trick"
// Lognormal prior parameterization on population parameters to match Metrum


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
  
  real partial_sum_lpmf(array[] int seq_subj, int start, int end,
                        vector dv_obs, array[] int dv_obs_id, array[] int i_obs,
                        array[] real amt, array[] int cmt, array[] int evid, 
                        array[] real time, array[] real rate, array[] real ii, 
                        array[] int addl, array[] int ss,
                        array[] int subj_start, array[] int subj_end, 
                        real TVCL, real TVVC, real TVQ, real TVVP, real TVKA,
                        vector omega, matrix L, matrix Z, 
                        real sigma,
                        vector lloq, array[] int bloq,
                        int n_random, int n_subjects, int n_total){
                           
    real ptarget = 0;
    row_vector[n_random] typical_values = 
      to_row_vector({TVCL, TVVC, TVQ, TVVP, TVKA});
    
    matrix[n_subjects, n_random] eta = diag_pre_multiply(omega, L * Z)';

    matrix[n_subjects, n_random] theta =
                          (rep_matrix(typical_values, n_subjects) .* exp(eta));
    
                              
    int N = end - start + 1;    // number of subjects in this slice  
    vector[n_total] dv_ipred;   // = rep_vector(0, n_total);
    matrix[n_total, 3] x_ipred; // = rep_matrix(0, n_total, 2);
  

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
      
      row_vector[n_random] theta_nn = theta[nn]; // access the parameters for subject nn
      real cl = theta_nn[1];
      real vc = theta_nn[2];
      real q = theta_nn[3];
      real vp = theta_nn[4];
      real ka = theta_nn[5];
      array[5] real theta_params = {cl, q, vc, vp, ka};
      
      x_ipred[subj_start[nn]:subj_end[nn], ] =
        pmx_solve_twocpt(time[subj_start[nn]:subj_end[nn]],
                         amt[subj_start[nn]:subj_end[nn]],
                         rate[subj_start[nn]:subj_end[nn]],
                         ii[subj_start[nn]:subj_end[nn]],
                         evid[subj_start[nn]:subj_end[nn]],
                         cmt[subj_start[nn]:subj_end[nn]],
                         addl[subj_start[nn]:subj_end[nn]],
                         ss[subj_start[nn]:subj_end[nn]],
                         theta_params)';
                      
      dv_ipred[subj_start[nn]:subj_end[nn]] = 
        x_ipred[subj_start[nn]:subj_end[nn], 2] ./ vc;
    
    }
  
    
    ipred_slice = dv_ipred[i_obs_slice];
    
    for(i in 1:n_obs_slice){
      if(bloq_slice[i] == 1){
        ptarget += lognormal_lcdf(lloq_slice[i] | log(ipred_slice[i]), sigma);
      }else{
        ptarget += lognormal_lpdf(dv_obs_slice[i] | log(ipred_slice[i]), sigma);
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
  
  real<lower = 0> location_tvcl;   // Prior Location parameter for CL
  real<lower = 0> location_tvvc;   // Prior Location parameter for VC
  real<lower = 0> location_tvq;    // Prior Location parameter for Q
  real<lower = 0> location_tvvp;   // Prior Location parameter for VP
  real<lower = 0> location_tvka;   // Prior Location parameter for KA
  
  real<lower = 0> scale_tvcl;      // Prior Scale parameter for CL
  real<lower = 0> scale_tvvc;      // Prior Scale parameter for VC
  real<lower = 0> scale_tvq;       // Prior Scale parameter for Q
  real<lower = 0> scale_tvvp;      // Prior Scale parameter for VP
  real<lower = 0> scale_tvka;      // Prior Scale parameter for KA
  
  // real<lower = 0> scale_omega_cl; // Prior scale parameter for omega_cl
  // real<lower = 0> scale_omega_vc; // Prior scale parameter for omega_vc
  // real<lower = 0> scale_omega_q;  // Prior scale parameter for omega_q
  // real<lower = 0> scale_omega_vp; // Prior scale parameter for omega_vp
  // real<lower = 0> scale_omega_ka; // Prior scale parameter for omega_ka
  
  real<lower = 0> lkj_df_omega;   // Prior degrees of freedom for omega cor mat
  
  real<lower = 0> scale_sigma;    // Prior Scale parameter for lognormal error
 
}
transformed data{ 
  
  int grainsize = 1;
  
  vector<lower = 0>[n_obs] dv_obs = dv[i_obs];
  array[n_obs] int dv_obs_id = ID[i_obs];
  
  vector[n_obs] lloq_obs = lloq[i_obs];
  array[n_obs] int bloq_obs = bloq[i_obs];
  
  int n_random = 5;                    // Number of random effects
  
  // array[n_random] real scale_omega = {scale_omega_cl, scale_omega_vc, 
  //                                     scale_omega_q, scale_omega_vp,
  //                                     scale_omega_ka}; 
  
  array[n_subjects] int seq_subj = sequence(1, n_subjects); // reduce_sum over subjects
  
}
parameters{ 
  
  real<lower = 0> TVCL;       
  real<lower = 0> TVVC; 
  real<lower = 0> TVQ;
  real<lower = 0> TVVP;
  real<lower = 0> TVKA;
  
  vector<lower = 0.01, upper = 2>[n_random] omega;
  // vector<lower = 0>[n_random] omega;
  cholesky_factor_corr[n_random] L;
  
  real<lower = 0> sigma;
  
  matrix[n_random, n_subjects] Z;
  
}

model{ 
  
  // Priors
  TVCL ~ lognormal(log(location_tvcl), scale_tvcl);
  TVVC ~ lognormal(log(location_tvvc), scale_tvvc);
  TVQ ~ lognormal(log(location_tvq), scale_tvq);
  TVVP ~ lognormal(log(location_tvvp), scale_tvvp);
  TVKA ~ lognormal(log(location_tvka), scale_tvka);

  // omega ~ normal(0, scale_omega);
  L ~ lkj_corr_cholesky(lkj_df_omega);
  
  sigma ~ cauchy(0, scale_sigma);
  
  to_vector(Z) ~ std_normal();
  
  // Likelihood
  target += reduce_sum(partial_sum_lupmf, seq_subj, grainsize,
                       dv_obs, dv_obs_id, i_obs,
                       amt, cmt, evid, time, 
                       rate, ii, addl, ss, subj_start, subj_end, 
                       TVCL, TVVC, TVQ, TVVP, TVKA, omega, L, Z,
                       sigma,
                       lloq, bloq,
                       n_random, n_subjects, n_total);
}
generated quantities{
  
  real<lower = 0> sigma_sq = square(sigma);
  
  real<lower = 0> omega_cl = omega[1];
  real<lower = 0> omega_vc = omega[2];
  real<lower = 0> omega_q = omega[3];
  real<lower = 0> omega_vp = omega[4];
  real<lower = 0> omega_ka = omega[5];

  real<lower = 0> omega_sq_cl = square(omega_cl);
  real<lower = 0> omega_sq_vc = square(omega_vc);
  real<lower = 0> omega_sq_q = square(omega_q);
  real<lower = 0> omega_sq_vp = square(omega_vp);
  real<lower = 0> omega_sq_ka = square(omega_ka);

  real cor_cl_vc;
  real cor_cl_q;
  real cor_cl_vp;
  real cor_cl_ka;
  real cor_vc_q;
  real cor_vc_vp;
  real cor_vc_ka;
  real cor_q_vp;
  real cor_q_ka;
  real cor_vp_ka;
  
  real omega_cl_vc;
  real omega_cl_q;
  real omega_cl_vp;
  real omega_cl_ka;
  real omega_vc_q;
  real omega_vc_vp;
  real omega_vc_ka;
  real omega_q_vp;
  real omega_q_ka;
  real omega_vp_ka;

  vector[n_subjects] eta_cl;
  vector[n_subjects] eta_vc;
  vector[n_subjects] eta_q;
  vector[n_subjects] eta_vp;
  vector[n_subjects] eta_ka;
  vector[n_subjects] CL;
  vector[n_subjects] VC;
  vector[n_subjects] Q;
  vector[n_subjects] VP;
  vector[n_subjects] KA;
  
  vector[n_obs] ipred;
  vector[n_obs] pred;
  vector[n_obs] dv_ppc;
  vector[n_obs] log_lik;
  vector[n_obs] res;
  vector[n_obs] wres;
  vector[n_obs] ires;
  vector[n_obs] iwres;
 
  {
    row_vector[n_random] typical_values = 
      to_row_vector({TVCL, TVVC, TVQ, TVVP, TVKA});

    matrix[n_random, n_random] R = multiply_lower_tri_self_transpose(L);
    matrix[n_random, n_random] Omega = quad_form_diag(R, omega);

    matrix[n_subjects, n_random] eta = diag_pre_multiply(omega, L * Z)';

    matrix[n_subjects, n_random] theta =
                          (rep_matrix(typical_values, n_subjects) .* exp(eta));

    vector[n_total] dv_pred;
    matrix[n_total, 3] x_pred;
    vector[n_total] dv_ipred;
    matrix[n_total, 3] x_ipred;
    
    array[5] real theta_params_tv = {TVCL, TVQ, TVVC, TVVP, TVKA};
    
    eta_cl = col(eta, 1);
    eta_vc = col(eta, 2);
    eta_q = col(eta, 3);
    eta_vp = col(eta, 4);
    eta_ka = col(eta, 5);
    
    CL = col(theta, 1);
    VC = col(theta, 2);
    Q = col(theta, 3);
    VP = col(theta, 4);
    KA = col(theta, 5);

    cor_cl_vc = R[1, 2];
    cor_cl_q = R[1, 3];
    cor_cl_vp = R[1, 4];
    cor_cl_ka = R[1, 5];
    cor_vc_q = R[2, 3];
    cor_vc_vp = R[2, 4];
    cor_vc_ka = R[2, 5];
    cor_q_vp = R[3, 4];
    cor_q_ka = R[3, 5];
    cor_vp_ka = R[4, 5];
    
    omega_cl_vc = Omega[1, 2];
    omega_cl_q = Omega[1, 3];
    omega_cl_vp = Omega[1, 4];
    omega_cl_ka = Omega[1, 5];
    omega_vc_q = Omega[2, 3];
    omega_vc_vp = Omega[2, 4];
    omega_vc_ka = Omega[2, 5];
    omega_q_vp = Omega[3, 4];
    omega_q_ka = Omega[3, 5];
    omega_vp_ka = Omega[4, 5];

    for(j in 1:n_subjects){

      array[5] real theta_params = {CL[j], Q[j], VC[j], VP[j], KA[j]};
      
      x_ipred[subj_start[j]:subj_end[j],] =
        pmx_solve_twocpt(time[subj_start[j]:subj_end[j]],
                         amt[subj_start[j]:subj_end[j]],
                         rate[subj_start[j]:subj_end[j]],
                         ii[subj_start[j]:subj_end[j]],
                         evid[subj_start[j]:subj_end[j]],
                         cmt[subj_start[j]:subj_end[j]],
                         addl[subj_start[j]:subj_end[j]],
                         ss[subj_start[j]:subj_end[j]],
                         theta_params)';
                      
                      
      dv_ipred[subj_start[j]:subj_end[j]] = 
        x_ipred[subj_start[j]:subj_end[j], 2] ./ VC[j];

      x_pred[subj_start[j]:subj_end[j],] =
        pmx_solve_twocpt(time[subj_start[j]:subj_end[j]],
                         amt[subj_start[j]:subj_end[j]],
                         rate[subj_start[j]:subj_end[j]],
                         ii[subj_start[j]:subj_end[j]],
                         evid[subj_start[j]:subj_end[j]],
                         cmt[subj_start[j]:subj_end[j]],
                         addl[subj_start[j]:subj_end[j]],
                         ss[subj_start[j]:subj_end[j]],
                         theta_params_tv)';

      dv_pred[subj_start[j]:subj_end[j]] = 
        x_pred[subj_start[j]:subj_end[j], 2] ./ TVVC;
    }

    pred = dv_pred[i_obs];
    ipred = dv_ipred[i_obs];

  }

  res = dv_obs - pred;                                                                             
  ires = dv_obs - ipred;

  for(i in 1:n_obs){
    dv_ppc[i] = lognormal_rng(log(ipred[i]), sigma);
    if(bloq_obs[i] == 1){
      log_lik[i] = lognormal_lcdf(lloq_obs[i] | log(ipred[i]), sigma);
    }else{
      log_lik[i] = lognormal_lpdf(dv_obs[i] | log(ipred[i]), sigma);
    }
    wres[i] = res[i]/sigma;
    iwres[i] = ires[i]/sigma;
  }
  
}

