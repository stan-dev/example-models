functions {
#include "utils.stan"
#include "model_lib.stan"
  
  // we fit a 1-cmt oral dosing situation
  matrix pk_system(vector lref, vector Dt, vector theta) {
    // as we fitting a 1-cmt oral dosing situation such that k1=k12 (all
    // mass out of 1 goes to 2)
    return(pk_1cmt_metabolite(lref, Dt, theta[1], theta[1], theta[2], 0, 0));
  }

  // note that n are the additional doses to be added such that in total
  // n+1 are added
  matrix pk_system_addl(vector lref, vector Dt, int cmt, real lamt, real tau, int n, vector theta) {
    matrix[num_elements(Dt), num_elements(lref)] lstate;
    matrix[num_elements(Dt), num_elements(lref)] lstate_mdose;
    vector[num_elements(lref)] lref_mdose;
    int S;
    
    // evolve reference state freely...
    lstate <- pk_system(lref, Dt, theta);
    
    // ... and add the extra doses correctly time-shifted
    S <- num_elements(lref);
    lref_mdose <- rep_vector(-35, S);
    lref_mdose[cmt] <- lamt;
    //if(prod(Dt - tau * n) < 0) reject("All requested times must be past the last dosing addl event.");
    /*
      for(i in 1:num_elements(Dt))
      if((Dt[i] - tau * n) < 0)
    reject();
    ("All requested times must be past the last dosing addl event.");
    */
    lstate_mdose <- pk_1cmt_metabolite(lref_mdose, Dt - tau * n, theta[1], theta[1], theta[2], tau, n+1);
    for(s in 1:S)
      for(t in 1:num_elements(Dt))
        lstate[t,s] <- log_sum_exp(lstate_mdose[t,s], lstate[t,s]);
    return(lstate);
  }
  
}
data {
  int<lower = 1> N; // number of lines of nm data set
  vector<lower=0>[N] time;
  vector<lower=0>[N] amt;
  int cmt[N];
  int<lower=0, upper=1> mdv[N];
  int<lower=0, upper=2> evid[N];
  int<lower=1, upper=N> id[N];
  int<lower=0> addl[N];
  vector<lower=0>[N] tau;

  vector<lower=0>[N] dv; // observations

  vector[3] prior_theta_mean;
  vector<lower=0>[3] prior_theta_sd;
}
transformed data {
  int dose_ind[count_elem(evid, 1)];
  int obs_ind[count_elem(mdv, 0)];
  int obs_M[rle_elem_count(id)];
  int dose_M[rle_elem_count(id)];
  int obs_time_rank[count_elem(mdv, 0)];
  vector[count_elem(mdv, 0)] obs_time;
  vector[count_elem(mdv, 0)] obs_ldv;
  int obs_dose_given[count_elem(mdv, 0)];
  int obs_cmt[count_elem(mdv, 0)];
  vector[count_elem(evid, 1)] dose_time;
  vector[count_elem(evid, 1)] dose_tau;
  int dose_addl[count_elem(evid, 1)];
  vector[count_elem(evid, 1)] dose_lamt;
  int dose_cmt[count_elem(evid, 1)];
  int J;
  int O;
  row_vector[rle_elem_count(id)] zero;
  matrix[rle_elem_count(id),2] Init_lstate;
  vector[rle_elem_count(id)] init_time;

  dose_ind <- which_elem(evid, 1);
  obs_ind  <- which_elem(mdv , 0);
  
  // note: We implicitly assume here that every patient has at least
  // one dose. If not, this line breaks, but in any PK problem this
  // will be given.
  dose_M <- rle_int(id[dose_ind]);
  obs_M  <- rle_int(id[obs_ind] );

  obs_time <- time[obs_ind];
  obs_ldv  <- log(dv[obs_ind]);
  obs_cmt  <- cmt[obs_ind];

  dose_time <- time[dose_ind];
  dose_tau <- tau[dose_ind];
  dose_lamt <- log(amt[dose_ind]);
  dose_addl <- addl[dose_ind];
  dose_cmt <- cmt[dose_ind];
  
  obs_time_rank <- find_interval_blocked(obs_M, obs_time, dose_M, dose_time);
  obs_dose_given <- count_dose_given_blocked(obs_M, obs_time, dose_M, dose_time, dose_tau, dose_addl);
  
  J <- rle_elem_count(id);
  O <- count_elem(mdv, 0);

  zero <- rep_row_vector(0, J);

  Init_lstate <- rep_matrix(-25, J, 2);
  init_time   <- rep_vector(0, J);
}
parameters {
  ordered[2] theta_lelim;
  real theta_lV;
  vector<lower=0>[2] omega;
  matrix[2,J] xi;
  real<lower=0> sigma_y;
}
transformed parameters {
  vector[3] theta;
  matrix[3,J] Theta;
  row_vector<lower=0>[J] kDelta;

  // log(ka)
  theta[1] <- theta_lelim[2]; // absorbtion must be faster (hence larger) than elimination
  // log(ke)
  theta[2] <- theta_lelim[1];
  // log(V)
  theta[3] <- theta_lV;
  
  Theta[1] <- rep_row_vector(theta[1], J);
  // ncp parametrization
  //Theta[2:3] <- rep_matrix(theta[2:3], J) + diag_pre_multiply(omega, xi);
  // cp parametrization
  Theta[2:3] <- xi;

  // kDelta is only defined to ensure that we have a faster absorption
  // than elimination (avoid "flip-flop") for each patient
  kDelta <- Theta[1] - Theta[2];
}
model {
  vector[O] ipred;

  theta_lelim[2] ~ normal(prior_theta_mean[1], prior_theta_sd[1]);
  theta_lelim[1] ~ normal(prior_theta_mean[2], prior_theta_sd[2]);
  theta_lV       ~ normal(prior_theta_mean[3], prior_theta_sd[3]);

  // ncp parametrization
  //to_vector(xi) ~ normal(0, 1);
  // cp parametrization
  xi[1] ~ normal(theta[2], omega[1]);
  xi[2] ~ normal(theta[3], omega[2]);
  
  omega ~ normal(0, 1);
  sigma_y ~ normal(0, 1);

  {
    matrix[2,J] Lscale;
    matrix[O,2] ly;

    Lscale[1] <- zero;
    Lscale[2] <- Theta[3];
    
    ly <- evaluate_model_fast(dose_M, dose_lamt, dose_cmt, dose_time, dose_tau, dose_addl,
                              Init_lstate, init_time,
                              obs_M, obs_time, obs_time_rank, obs_dose_given,
                              Theta[1:2]',
                              Lscale');
    for (i in 1:O)
      ipred[i] <- ly[i, obs_cmt[i]];
  }

  obs_ldv ~ normal(ipred, sigma_y);
}
generated quantities {
}
