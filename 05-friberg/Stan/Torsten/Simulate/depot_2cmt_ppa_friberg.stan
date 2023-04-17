// First Order Absorption (oral/subcutaneous)
// Two-compartment PK Model with linear elimination
// Friberg-Karlsson model for myelosuppression
// IIV on CL, VC, Q, VP, KA, MTT, CIRC0, GAMMA, ALPHA
// proportional error - DV = CP(1 + eps_p) for both PK and PD
// ODE solution using Torsten - full ODE system
// Since we have a normal distribution on the error, but the DV must be > 0, it
//   generates values from a normal that is truncated below at 0

functions{
  
  real normal_lb_rng(real mu, real sigma, real lb){
    
    real p_lb = normal_cdf(lb | mu, sigma);
    real u = uniform_rng(p_lb, 1);
    real y = mu + sigma * inv_Phi(u);
    return y;

  }
  
  vector depot_2cmt_friberg_ode(real t, vector y, array[] real params, 
                                array[] real x_r, array[] int x_i){
    
    real cl = params[1];
    real vc = params[2];
    real q = params[3];
    real vp = params[4];
    real ka = params[5];
    real mtt = params[6];
    real circ_0 = params[7];
    real gamma = params[8];
    real alpha = params[9];
    
    real ke = cl/vc;
    real k_cp = q/vc;
    real k_pc = q/vp;
    real k_tr = 4/mtt; // k_tr = (n_tr + 1)/mtt    
    
    real conc = y[2]/vc;
    
    real e_drug = fmin(1.0, alpha*conc); // Maybe reparameterize this so no more fmin?
    real prol = y[4] + circ_0;
    real transit_1 = y[5] + circ_0; 
    real transit_2 = y[6] + circ_0;
    real transit_3 = y[7] + circ_0;
    real circ = y[8] + circ_0; // fmax(machine_precision(), y[8] + circ_0)
    
    vector[8] dydt;
    
    dydt[1] = -ka*y[1];                               // depot
    dydt[2] = ka*y[1] - (ke + k_cp)*y[2] + k_pc*y[3]; // central
    dydt[3] = k_cp*y[2] - k_pc*y[3];                  // peripheral
    dydt[4] = k_tr*prol*((1 - e_drug)*(circ_0/circ)^gamma - 1);  // proliferative cells
    dydt[5] = k_tr*(prol - transit_1);                // transit 1
    dydt[6] = k_tr*(transit_1 - transit_2);           // transit 2
    dydt[7] = k_tr*(transit_2 - transit_3);           // transit 3
    dydt[8] = k_tr*(transit_3 - circ);                // circulating blood cells
    
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
  real<lower = 0> TVQ;
  real<lower = 0> TVVP;
  real<lower = 0> TVKA;
  real<lower = 0> TVMTT;
  real<lower = 0> TVCIRC0;
  real<lower = 0> TVGAMMA;
  real<lower = 0> TVALPHA;
  
  real<lower = 0> omega_cl;
  real<lower = 0> omega_vc;
  real<lower = 0> omega_q;
  real<lower = 0> omega_vp;
  real<lower = 0> omega_ka;
  real<lower = 0> omega_mtt;
  real<lower = 0> omega_circ_0;
  real<lower = 0> omega_gamma;
  real<lower = 0> omega_alpha;
  
  corr_matrix[9] R;  // Correlation matrix before transforming to Omega.
                     // Can in theory change this to having inputs for
                     // cor_cl_vc, ... and then construct the 
                     // correlation matrix in transformed data, but it's easy
                     // enough to do in R
  
  real<lower = 0> sigma_p;
  real<lower = 0> sigma_a;
  real<lower = -1, upper = 1> cor_p_a;
  
  real<lower = 0> sigma_p_pd;
  real<lower = 0> sigma_a_pd;
  real<lower = -1, upper = 1> cor_p_a_pd;
  
}
transformed data{
  
  int n_random = 9;
  int n_cmt = 8;

  vector[n_random] omega = [omega_cl, omega_vc, omega_q, omega_vp, omega_ka,
                            omega_mtt, omega_circ_0, omega_gamma, omega_alpha]';
  
  matrix[n_random, n_random] L = cholesky_decompose(R);

  vector[2] sigma = [sigma_p, sigma_a]';
  matrix[2, 2] R_Sigma = rep_matrix(1, 2, 2);
  R_Sigma[1, 2] = cor_p_a;
  R_Sigma[2, 1] = cor_p_a;
  
  vector[2] sigma_pd = [sigma_p_pd, sigma_a_pd]';
  matrix[2, 2] R_Sigma_pd = rep_matrix(1, 2, 2);
  R_Sigma_pd[1, 2] = cor_p_a_pd;
  R_Sigma_pd[2, 1] = cor_p_a_pd;
  
  matrix[2, 2] Sigma = quad_form_diag(R_Sigma, sigma);
  matrix[2, 2] Sigma_pd = quad_form_diag(R_Sigma_pd, sigma_pd);
  
  array[n_subjects] real bioav = rep_array(1.0, n_subjects); // Hardcoding, but could be data or a parameter in another situation
  array[n_subjects] real t_lag = rep_array(1.0, n_subjects); // Hardcoding, but could be data or a parameter in another situation
  
}
model{
  
}
generated quantities{
  
  vector[n_total] dv; // dependent variable with residual error
  vector[n_total] cp; // dependent variable with no residual error
  
  vector[n_subjects] CL;
  vector[n_subjects] VC;
  vector[n_subjects] Q;
  vector[n_subjects] VP;
  vector[n_subjects] KA;
  vector[n_subjects] MTT;
  vector[n_subjects] CIRC0;
  vector[n_subjects] GAMMA;
  vector[n_subjects] ALPHA;
  
  {
  
    vector[n_random] typical_values = 
      to_vector({TVCL, TVVC, TVQ, TVVP, TVKA,
                 TVMTT, TVCIRC0, TVGAMMA, TVALPHA});
    
    matrix[n_random, n_subjects] eta;   
    matrix[n_subjects, n_random] theta; 
    matrix[n_total, n_cmt] x_cp;
    
    for(i in 1:n_subjects){
      eta[, i] = multi_normal_cholesky_rng(rep_vector(0, n_random),
                                           diag_pre_multiply(omega, L));
    }
    theta = (rep_matrix(typical_values, n_subjects) .* exp(eta))';

    CL = col(theta, 1);
    VC = col(theta, 2);
    Q = col(theta, 3);
    VP = col(theta, 4);
    KA = col(theta, 5);
    MTT = col(theta, 6);
    CIRC0 = col(theta, 7);
    GAMMA = col(theta, 8);
    ALPHA = col(theta, 9);
    
    for(j in 1:n_subjects){
      
      array[n_random] real params = {CL[j], VC[j], Q[j], VP[j], KA[j],
                                     MTT[j], CIRC0[j], GAMMA[j], ALPHA[j]};
                              
      x_cp[subj_start[j]:subj_end[j],] =
        pmx_solve_rk45(depot_2cmt_friberg_ode,
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
      
      for(k in subj_start[j]:subj_end[j]){
        if(cmt[k] == 2){                     // PK observation => cmt = 2  
          cp[k] = x_cp[k, 2] / VC[j];
        }else if(cmt[k] == 3){               // PD observation => cmt = 3
          cp[k] = x_cp[k, 8] + CIRC0[j];
        }
      }
    
    }

    for(i in 1:n_total){
      if(cp[i] == 0){
         dv[i] = 0;
      }else{
        if(cmt[i] == 2){
          real cp_tmp = cp[i];
          real sigma_tmp = sqrt(square(cp_tmp) * Sigma[1, 1] + Sigma[2, 2] + 
                           2*cp_tmp*Sigma[2, 1]);  
          dv[i] = normal_lb_rng(cp_tmp, sigma_tmp, 0.0);
        }else if(cmt[i] == 3){
          real cp_tmp = cp[i];
          real sigma_tmp = sqrt(square(cp_tmp) * Sigma_pd[1, 1] + Sigma_pd[2, 2] + 
                           2*cp_tmp*Sigma_pd[2, 1]);
          dv[i] = normal_lb_rng(cp_tmp, sigma_tmp, 0.0);
        }
      }
    }
  }
}

