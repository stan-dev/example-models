functions {
  array[] real two_pool_feedback(real t, array[] real C, array[] real theta,
                                 array[] real x_r, array[] int x_i) {
    real k1;
    real k2;
    real alpha21;
    real alpha12;
    
    array[2] real dC_dt;
    
    k1 = theta[1];
    k2 = theta[2];
    alpha21 = theta[3];
    alpha12 = theta[4];
    
    dC_dt[1] = -k1 * C[1] + alpha12 * k2 * C[2];
    dC_dt[2] = alpha21 * k1 * C[1] - k2 * C[2];
    
    return dC_dt;
  }
  
  array[] real evolved_CO2(int N_t, real t0, array[] real ts, real gamma,
                           real totalC_t0, real k1, real k2, real alpha21,
                           real alpha12, data array[] real x_r,
                           data array[] int x_i) {
    array[2] real C_t0; // initial state
    array[4] real theta; // ODE parameters
    array[N_t, 2] real C_hat; // predicted pool content
    
    array[N_t] real eCO2_hat;
    
    C_t0[1] = gamma * totalC_t0;
    C_t0[2] = (1 - gamma) * totalC_t0;
    
    theta[1] = k1;
    theta[2] = k2;
    theta[3] = alpha21;
    theta[4] = alpha12;
    
    C_hat = integrate_ode(two_pool_feedback, C_t0, t0, ts, theta, x_r, x_i);
    
    for (t in 1 : N_t) {
      eCO2_hat[t] = totalC_t0 - sum(C_hat[t]);
    }
    return eCO2_hat;
  }
}
data {
  real<lower=0> totalC_t0; // initial total carbon
  
  real t0; // initial time
  int<lower=0> N_t; // number of measurement times
  array[N_t] real<lower=t0> ts; // measurement times
  
  vector<lower=0>[N_t] eCO2mean; // measured cumulative evolved carbon
  array[N_t] real<lower=0> eCO2sd; // measured cumulative evolved carbon
}
transformed data {
  array[0] real x_r; // no real data for ODE system
  array[0] int x_i; // no integer data for ODE system
}
parameters {
  real<lower=0> k1; // pool 1 decomposition rate
  real<lower=0> k2; // pool 2 decomposition rate
  
  real<lower=0> alpha21; // transfer coeff from pool 2 to 1
  real<lower=0> alpha12; // transfer coeff from pool 1 to 2
  
  real<lower=0, upper=1> gamma; // partitioning coefficient
  
  real<lower=0> sigma; // observation std dev
  
  vector<lower=0>[N_t] eCO2; // evolved CO2
}
transformed parameters {
  array[N_t] real eCO2_hat;
  eCO2_hat = evolved_CO2(N_t, t0, ts, gamma, totalC_t0, k1, k2, alpha21,
                         alpha12, x_r, x_i);
}
model {
  gamma ~ beta(10, 1); // identifies pools
  k1 ~ normal(0, 1); // weakly informative
  k2 ~ normal(0, 1);
  alpha21 ~ normal(0, 1);
  alpha12 ~ normal(0, 1);
  sigma ~ normal(0, 1);
  
  // likelihood marginalizes out eCO2
  
  // true eCO2 (unknown) generated from system plus noise
  eCO2 ~ normal(eCO2_hat, sigma);
  
  // measurement eCO2mean (observed) of true eCO2 value with
  // measurement noise eCO2sd (observed)
  eCO2mean ~ normal(eCO2, eCO2sd);
}
