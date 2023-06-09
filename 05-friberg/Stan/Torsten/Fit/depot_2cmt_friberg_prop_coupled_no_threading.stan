// First Order Absorption (oral/subcutaneous)
// Two-compartment PK Model with linear elimination
// Friberg-Karlsson model for myelosuppression
// IIV on CL, VC, Q, VP, KA, MTT, CIRC0, GAMMA, ALPHA
// proportional error - DV = IPRED*(1 + eps_p) for both PK and PD
// ODE solution using Torsten - full ODE system
// Deals with BLOQ values by the "CDF trick" (M4)
// Since we have a normal distribution on the error, but the DV must be > 0, it
//   truncates the likelihood below at 0

functions{

  real normal_lb_rng(real mu, real sigma, real lb){
    
    real p_lb = normal_cdf(lb | mu, sigma);
    real u = uniform_rng(p_lb, 1);
    real y = mu + sigma * inv_Phi(u);
    return y;

  }

  vector depot_2cmt_friberg_ode(real t, vector y, vector y_pk,
                                array[] real params, array[] real x_r, 
                                array[] int x_i){
    
    real vc = params[3];
    real mtt = params[6];
    real circ_0 = params[7];
    real gamma = params[8];
    real alpha = params[9];
    
    real conc = y_pk[2]/vc;
    
    real k_tr = 4/mtt; // k_tr = (n_tr + 1)/mtt    
    real e_drug = fmin(1.0, alpha*conc); // Maybe reparameterize this so no more fmin?
    real prol = y[1] + circ_0;
    real transit_1 = y[2] + circ_0; 
    real transit_2 = y[3] + circ_0;
    real transit_3 = y[4] + circ_0;
    real circ = y[5] + circ_0; // fmax(machine_precision(), y[5] + circ_0)
    
    vector[5] dydt;
    
    dydt[1] = k_tr*prol*((1 - e_drug)*(circ_0/circ)^gamma - 1);  // proliferative cells
    dydt[2] = k_tr*(prol - transit_1);                // transit 1
    dydt[3] = k_tr*(transit_1 - transit_2);           // transit 2
    dydt[4] = k_tr*(transit_2 - transit_3);           // transit 3
    dydt[5] = k_tr*(transit_3 - circ);                // circulating blood cells
    
    return dydt;
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
  
  real<lower = 0> location_tvcl;    // Prior Location parameter for CL
  real<lower = 0> location_tvvc;    // Prior Location parameter for VC
  real<lower = 0> location_tvq;     // Prior Location parameter for Q
  real<lower = 0> location_tvvp;    // Prior Location parameter for VP
  real<lower = 0> location_tvka;    // Prior Location parameter for KA
  real<lower = 0> location_tvmtt;   // Prior Location parameter for MTT
  real<lower = 0> location_tvcirc0; // Prior Location parameter for CIRC0
  real<lower = 0> location_tvgamma; // Prior Location parameter for GAMMA
  real<lower = 0> location_tvalpha; // Prior Location parameter for ALPHA
  
  real<lower = 0> scale_tvcl;      // Prior Scale parameter for CL
  real<lower = 0> scale_tvvc;      // Prior Scale parameter for VC
  real<lower = 0> scale_tvq;       // Prior Scale parameter for Q
  real<lower = 0> scale_tvvp;      // Prior Scale parameter for VP
  real<lower = 0> scale_tvka;      // Prior Scale parameter for KA
  real<lower = 0> scale_tvmtt;     // Prior Scale parameter for MTT
  real<lower = 0> scale_tvcirc0;   // Prior Scale parameter for CIRC0
  real<lower = 0> scale_tvgamma;   // Prior Scale parameter for GAMMA
  real<lower = 0> scale_tvalpha;   // Prior Scale parameter for ALPHA
  
  real<lower = 0> scale_omega_cl;    // Prior scale parameter for omega_cl
  real<lower = 0> scale_omega_vc;    // Prior scale parameter for omega_vc
  real<lower = 0> scale_omega_q;     // Prior scale parameter for omega_q
  real<lower = 0> scale_omega_vp;    // Prior scale parameter for omega_vp
  real<lower = 0> scale_omega_ka;    // Prior scale parameter for omega_ka
  real<lower = 0> scale_omega_mtt;   // Prior scale parameter for omega_mtt
  real<lower = 0> scale_omega_circ0; // Prior scale parameter for omega_circ0
  real<lower = 0> scale_omega_gamma; // Prior scale parameter for omega_gamma
  real<lower = 0> scale_omega_alpha; // Prior scale parameter for omega_alpha
  
  real<lower = 0> lkj_df_omega;   // Prior degrees of freedom for omega cor mat
  
  real<lower = 0> scale_sigma_p;     // Prior Scale parameter for proportional error for PK
  real<lower = 0> scale_sigma_p_pd;  // Prior Scale parameter for proportional error for PD
  
}
transformed data{ 
  
  int grainsize = 1;
  
  vector<lower = 0>[n_obs] dv_obs = dv[i_obs];
  array[n_obs] int dv_obs_id = ID[i_obs];
  
  vector[n_obs] lloq_obs = lloq[i_obs];
  array[n_obs] int bloq_obs = bloq[i_obs];
  
  int n_random = 9;                    // Number of random effects
  int n_cmt = 3;                       // Number of compartments in PK system
  int n_ode = 5;                       // Number of compartments in PD system
  
  array[n_random] real scale_omega = 
    {scale_omega_cl, scale_omega_vc, scale_omega_q, scale_omega_vp, scale_omega_ka,
     scale_omega_mtt, scale_omega_circ0, scale_omega_gamma, scale_omega_alpha}; 
  
  array[n_cmt + n_ode] real bioav = rep_array(1.0, n_cmt + n_ode);	
  array[n_cmt + n_ode] real tlag = rep_array(0.0, n_cmt + n_ode);
  
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
  
  real<lower = 0> sigma_p;
  real<lower = 0> sigma_p_pd;
  
  matrix[n_random, n_subjects] Z;
  
}
transformed parameters{
  
  vector[n_subjects] CL;
  vector[n_subjects] VC;
  vector[n_subjects] Q;
  vector[n_subjects] VP;
  vector[n_subjects] KA;
  vector[n_subjects] MTT;
  vector[n_subjects] CIRC0;
  vector[n_subjects] GAMMA;
  vector[n_subjects] ALPHA;
  
  vector[n_obs] ipred;

  {

    row_vector[n_random] typical_values = 
      to_row_vector({TVCL, TVVC, TVQ, TVVP, TVKA,
                     TVMTT, TVCIRC0, TVGAMMA, TVALPHA});
    
    matrix[n_subjects, n_random] eta = diag_pre_multiply(omega, L * Z)';

    matrix[n_subjects, n_random] theta =
                          (rep_matrix(typical_values, n_subjects) .* exp(eta));
                          
    vector[n_total] dv_ipred;   
    matrix[n_total, n_cmt + n_ode] x_ipred;
                          
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
      
      x_ipred[subj_start[j]:subj_end[j],] =
        pmx_solve_twocpt_rk45(depot_2cmt_friberg_ode,
                              n_ode,
                              time[subj_start[j]:subj_end[j]],
                              amt[subj_start[j]:subj_end[j]],
                              rate[subj_start[j]:subj_end[j]],
                              ii[subj_start[j]:subj_end[j]],
                              evid[subj_start[j]:subj_end[j]],
                              cmt[subj_start[j]:subj_end[j]],
                              addl[subj_start[j]:subj_end[j]],
                              ss[subj_start[j]:subj_end[j]],
                              {CL[j], Q[j], VC[j], VP[j], KA[j],
                                      MTT[j], CIRC0[j], GAMMA[j], ALPHA[j]}, 
                              bioav, tlag)';
      
      for(k in subj_start[j]:subj_end[j]){
        if(cmt[k] == 2){                     // PK observation => cmt = 2  
          dv_ipred[k] = x_ipred[k, 2] / VC[j];
        }else if(cmt[k] == 3){               // PD observation => cmt = 3
          dv_ipred[k] = x_ipred[k, 8] + CIRC0[j];
        }
      }
    
    }
    
    ipred = dv_ipred[i_obs];
    
  }
  
}
model{ 
  
  // Priors
  TVCL ~ lognormal(log(location_tvcl), scale_tvcl);
  TVVC ~ lognormal(log(location_tvvc), scale_tvvc);
  TVQ ~ lognormal(log(location_tvq), scale_tvq);
  TVVP ~ lognormal(log(location_tvvp), scale_tvvp);
  // TVKA ~ lognormal(log(location_tvka), scale_tvka);
  TVKA ~ lognormal(log(location_tvka), scale_tvka) T[0.5*(TVCL/TVVC + TVQ/TVVC + TVQ/TVVP +
          sqrt((TVCL/TVVC + TVQ/TVVC + TVQ/TVVP)^2 - 4*TVCL/TVVC*TVQ/TVVP)), ];
  TVMTT ~ lognormal(log(location_tvmtt), scale_tvmtt);
  TVCIRC0 ~ lognormal(log(location_tvcirc0), scale_tvcirc0);
  TVGAMMA ~ lognormal(log(location_tvgamma), scale_tvgamma);
  TVALPHA ~ lognormal(log(location_tvalpha), scale_tvalpha);

  omega ~ normal(0, scale_omega);
  L ~ lkj_corr_cholesky(lkj_df_omega);
  
  sigma_p ~ normal(0, scale_sigma_p);
  sigma_p_pd ~ normal(0, scale_sigma_p_pd);
  
  to_vector(Z) ~ std_normal();
  
  // Likelihood
  for(i in 1:n_obs){
    if(cmt[i] == 2){
      real sigma_tmp = ipred[i]*sigma_p;
      if(bloq[i] == 1){
        target += log_diff_exp(normal_lcdf(lloq[i] | ipred[i], sigma_tmp),
                               normal_lcdf(0.0 | ipred[i], sigma_tmp)) -
                  normal_lccdf(0.0 | ipred[i], sigma_tmp); 
      }else{
        target += normal_lpdf(dv_obs[i] | ipred[i], sigma_tmp) -
                  normal_lccdf(0.0 | ipred[i], sigma_tmp);
      }
    }else if(cmt[i] == 3){
      real sigma_tmp = ipred[i]*sigma_p_pd;
      if(bloq[i] == 1){
        target += log_diff_exp(normal_lcdf(lloq[i] | ipred[i], sigma_tmp),
                               normal_lcdf(0.0 | ipred[i], sigma_tmp)) -
                  normal_lccdf(0.0 | ipred[i], sigma_tmp); 
      }else{
         target += normal_lpdf(dv_obs[i] | ipred[i], sigma_tmp) -
                   normal_lccdf(0.0 | ipred[i], sigma_tmp);
      }
    }                                         
  }  
}
