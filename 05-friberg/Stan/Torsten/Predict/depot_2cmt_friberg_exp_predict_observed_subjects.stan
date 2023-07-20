// First Order Absorption (oral/subcutaneous)
// Two-compartment PK Model with linear elimination
// Friberg-Karlsson model for myelosuppression
// IIV on CL, VC, Q, VP, KA, MTT, CIRC0, GAMMA, ALPHA
// exponential error - DV = CP*exp(eps) for both PK and PD
// ODE solution using Torsten - full ODE system

functions{
  
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
  
  int n_random = 9;                    // Number of random effects
  int n_cmt = 8;                       // Number of compartments in ODEs

}
parameters{ 
  
  real<lower = 0> TVCL;       
  real<lower = 0> TVVC; 
  real<lower = 0> TVQ;
  real<lower = 0> TVVP;
  // real<lower = 0> TVKA;
  real<lower = 0.5*(TVCL/TVVC + TVQ/TVVC + TVQ/TVVP +
    sqrt((TVCL/TVVC + TVQ/TVVC + TVQ/TVVP)^2 - 4*TVCL/TVVC*TVQ/TVVP))> TVKA;
  real<lower = 0> TVMTT;
  real<lower = 0> TVCIRC0;
  real<lower = 0> TVGAMMA;
  real<lower = 0> TVALPHA;
  
  vector<lower = 0>[n_random] omega;
  cholesky_factor_corr[n_random] L;
  
  real<lower = 0> sigma;
  real<lower = 0> sigma_pd;
  
  matrix[n_random, n_subjects] Z;
  
}
generated quantities{
  
  vector[n_time_new] ipred; // ipred for the observed individuals at the new timepoints
  vector[n_time_new] pred;  // pred for the observed individuals at the new timepoints
  vector[n_time_new] dv;    // dv for the observed individuals at the new timepoints
 
  {
    row_vector[n_random] typical_values = 
      to_row_vector({TVCL, TVVC, TVQ, TVVP, TVKA,
                     TVMTT, TVCIRC0, TVGAMMA, TVALPHA});

    matrix[n_random, n_random] R = multiply_lower_tri_self_transpose(L);
    matrix[n_random, n_random] Omega = quad_form_diag(R, omega);

    matrix[n_subjects, n_random] eta = diag_pre_multiply(omega, L * Z)';

    matrix[n_subjects, n_random] theta =
                          (rep_matrix(typical_values, n_subjects) .* exp(eta));
                          
    matrix[n_time_new, n_cmt] x_pred;
    matrix[n_time_new, n_cmt] x_ipred;
    
    array[n_random] real theta_params_tv = {TVCL, TVVC, TVQ, TVVP, TVKA,
                                            TVMTT, TVCIRC0, TVGAMMA, TVALPHA};
    
    vector[n_subjects] CL = col(theta, 1);
    vector[n_subjects] VC = col(theta, 2);
    vector[n_subjects] Q = col(theta, 3);
    vector[n_subjects] VP = col(theta, 4);
    vector[n_subjects] KA = col(theta, 5);
    vector[n_subjects] MTT = col(theta, 6);
    vector[n_subjects] CIRC0 = col(theta, 7);
    vector[n_subjects] GAMMA = col(theta, 8);
    vector[n_subjects] ALPHA = col(theta, 9);
    
    for(j in 1:n_subjects){

      array[n_random] real theta_params = {CL[j], VC[j], Q[j], VP[j], KA[j],
                                           MTT[j], CIRC0[j], GAMMA[j], ALPHA[j]};

      x_ipred[subj_start[j]:subj_end[j],] =
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
                       theta_params)';

     x_pred[subj_start[j]:subj_end[j],] =
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
                      theta_params_tv)';
        
      for(k in subj_start[j]:subj_end[j]){
        if(cmt[k] == 2){                     // PK observation => cmt = 2  
          ipred[k] = x_ipred[k, 2] / VC[j];
          pred[k] = x_pred[k, 2] / VC[j];
        }else if(cmt[k] == 3){               // PD observation => cmt = 3
          ipred[k] = x_ipred[k, 8] + CIRC0[j];
          pred[k] = x_pred[k, 8] + TVCIRC0;
        }
      }
    }
    
    for(i in 1:n_time_new){
      if(ipred[i] == 0){
        dv[i] = 0;
      }else{
        if(cmt[i] == 2){
          dv[i] = lognormal_rng(log(ipred[i]), sigma);
        }else if(cmt[i] == 3){
          dv[i] = lognormal_rng(log(ipred[i]), sigma_pd);
        }
      }
    }
  }
}


  