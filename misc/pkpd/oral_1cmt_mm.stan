functions {
#include "utils.stan"
#include "model_lib.stan"

  // returns dlog(y)/dt (=1/y dy/dt); first order absorbtion
  // calculcated analytically; parametrization is optimized for
  // efficient computation
  real[] pk_1cmt_mm_lode(real t, real[] ly, real[] theta, real[] x_r, int[] x_i) {
    real dly_dt[1];
    real ka;
    real k0;
    real lAm;
    real ct;

    ka <- theta[1];
    k0 <- theta[2];
    lAm <- theta[3];
    ct <- exp(theta[4] - t * ka - ly[1]);
    
    dly_dt[1] <- ka * ct - k0 * inv_logit(lAm - ly[1]);

    return(dly_dt);
  }

  
  // avoid to make the simple exponential part of the ODE which
  // decreases size of sensitivity system
  matrix pk_system(vector lref, vector Dt, vector theta, real[] x_r, int[] x_i) {
    matrix[num_elements(Dt),num_elements(lref)] sol;
    real int_sol[num_elements(Dt), 1];
    int P;
    real theta_tilde[num_elements(theta)+1];
    real lref_tilde[1];
    real ka;
    P <- num_elements(theta);
    theta_tilde[1:P] <- to_array_1d(theta);
    theta_tilde[P+1] <- lref[1];
    lref_tilde[1] <- lref[2];
    int_sol <- integrate_ode_rk45(pk_1cmt_mm_lode, lref_tilde, 0, to_array_1d(Dt), theta_tilde, x_r, x_i, 1e-4, 1e-4, 1000);
    ka <- theta[1];
    for(i in 1:num_elements(Dt)) {
      sol[i,1] <- lref[1] - ka * Dt[i];
      sol[i,2] <- int_sol[i,1];
    }
    return(sol);
  }

  // note that n are the additional doses to be added such that in total
  // n+1 are added
  matrix pk_system_addl(vector lref, vector Dt, int cmt, real lamt, real tau, int n, vector theta, real[] x_r, int[] x_i) {
    matrix[num_elements(Dt), num_elements(lref)] lstate;
    reject("ADDL dose coding not supported with ODEs!");
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

  vector[4] prior_theta_mean;
  vector<lower=0>[4] prior_theta_sd;
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
  int dose_next_obs[count_elem(mdv, 1)];
  int J;
  int O;
  row_vector[rle_elem_count(id)] zero;
  matrix[rle_elem_count(id),2] Init_lstate;
  vector[rle_elem_count(id)] init_time;
  real x_r[0];
  int x_i[0];

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

  dose_next_obs <- count_obs_event_free_blocked(obs_M, obs_time_rank, dose_M);
  
  J <- rle_elem_count(id);
  O <- count_elem(mdv, 0);

  zero <- rep_row_vector(0, J);

  // we initialize the main cmt to 0.1 in order to avoid too steep
  // derivatives when the first dose is injected
  Init_lstate <- append_col(rep_vector(-25, J), rep_vector(log(0.1), J));
  init_time   <- rep_vector(-1E-3, J);
}
parameters {
  ordered[2] theta_lelim;
  real theta_lAm;
  real theta_lV;
  vector<lower=0>[2] omega;
  matrix[2,J] xi;
  real<lower=0> sigma_y;
}
transformed parameters {
  vector[4] theta;
  matrix[3,J] Theta;
  row_vector<lower=0>[J] kDelta;

  // theta is on log-scale
  theta[1] <- theta_lelim[2];  // ka is larger than k0 (faster elimination than absorbtion)
  theta[2] <- theta_lelim[1];  // k0 = Vm/(V * Km)
  theta[3] <- theta_lAm;       // Am = Km*V
  theta[4] <- theta_lV;        // V

  // prepare parameters to pass into model function
  Theta[1] <- rep_row_vector(exp(theta[1]), J); // ka
  Theta[2] <- exp(xi[1]);                       // k0
  Theta[3] <- rep_row_vector(theta[3], J);      // log(Am)

  // kDelta is only defined to ensure that we have a faster absorption
  // than elimination (avoid "flip-flop") for each patient
  kDelta <- Theta[1] - Theta[2];
}
model {
  vector[O] ipred;

  theta_lelim[2] ~ normal(prior_theta_mean[1], prior_theta_sd[1]);
  theta_lelim[1] ~ normal(prior_theta_mean[2], prior_theta_sd[2]);
  theta_lAm      ~ normal(prior_theta_mean[3], prior_theta_sd[3]);
  theta_lV       ~ normal(prior_theta_mean[4], prior_theta_sd[4]);

  // cp parametrization
  xi[1] ~ normal(theta[2], omega[1]);
  xi[2] ~ normal(theta[4], omega[2]);
  
  omega ~ normal(0, 1);
  sigma_y ~ normal(0, 1);

  {
    matrix[O,2] ly;
    matrix[2,J] Lscale;

    Lscale[1] <- zero;
    Lscale[2] <- xi[2];

    ly <- evaluate_model_fast(dose_M, dose_lamt, dose_cmt, dose_time, dose_tau, dose_addl, dose_next_obs,
                              Init_lstate, init_time,
                              obs_M, obs_time, obs_time_rank, obs_dose_given, 
                              Theta',
                              Lscale',
                              x_r, x_i);
    for (i in 1:O)
      ipred[i] <- ly[i, obs_cmt[i]];
  }

  obs_ldv ~ normal(ipred, sigma_y);
}
generated quantities {
}
