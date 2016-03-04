functions {

  /**
   * ODE system for two pool model with feedback and no inputs.
   *
   * This is the version that does not deal with measurement error.
   *
   * System State C is two dimensional with C[1] and C[2] 
   * being carbon in pools 1 and 2.
   *
   * The system has parameters
   *
   *   theta = (k1, k2, alpha21, alpha12)
   *
   * where
   *
   *   k1:       pool 1 decomposition rate
   *   k2:       pool 2 decomposition rate
   *   alpha21:  transfer coefficient from pool 2 to pool 1
   *   alpha12:  transfer coefficient from pool 1 to pool 2
   *
   * The system time derivatives are
   *
   *   d.C[1] / d.t  =  -k1 * C[1]  +  alpha12 * k2 * C[2]
   *
   *   d.C[2] / d.t  =  alpha21 * k1 * C[1]  -  k2 * C[2]
   *
   * @param t time at which derivatives are evaluated.
   * @param C system state at which derivatives are evaluated.
   * @param theta parameters for system.
   * @param x_r real constants for system (empty).
   * @param x_i integer constants for system (empty).
   */
  real[] two_pool_feedback(real t, real[] C, real[] theta,
                           real[] x_r, int[] x_i) {
    real k1;
    real k2;
    real alpha21;
    real alpha12;

    real dC_dt[2];

    k1 <- theta[1];
    k2 <- theta[2];
    alpha21 <- theta[3];
    alpha12 <- theta[4];
    
    dC_dt[1] <- -k1 * C[1] + alpha12 * k2 * C[2];
    dC_dt[2] <- - k2 * C[2] + alpha21 * k1 * C[1] ;

    return dC_dt;
  }

  /**
   * Compute total evolved CO2 from the system given the specified
   * parameters and times.  This is done by simulating the system
   * defined by the ODE function two_pool_feedback and then
   * subtracting the sum of the CO2 estimated in each pool from the
   * initial CO2.
   *
   * @param T number of times.
   * @param t0 initial time.
   * @param ts observation times.
   * @param gamma partitioning coefficient.
   * @param k1 decomposition rate for pool 1
   * @param k2 decomposition rate for pool 2
   * @param alpha21 transfer coefficient from pool 2 to 1
   * @param alpha12 transfer coefficient from pool 1 to 2
   * @param x_r real data (empty)
   * @param x_i integer data (empty)
   * @return evolved CO2 for times ts
   */
  real[] evolved_CO2(int N_t, real t0, real[] ts,
                     real gamma, real totalC_t0,
                     real k1, real k2, 
                     real alpha21, real alpha12,
                     real[] x_r, int[] x_i) {

    real C_t0[2];               // initial state
    real theta[4];              // ODE parameters
    real C_hat[N_t,2];          // predicted pool content

    real eCO2_hat[N_t];

    C_t0[1] <- gamma * totalC_t0;
    C_t0[2] <- (1 - gamma) * totalC_t0;

    theta[1] <- k1;
    theta[2] <- k2;
    theta[3] <- alpha21;
    theta[4] <- alpha12;

    C_hat <- integrate_ode(two_pool_feedback, 
                           C_t0, t0, ts, theta, x_r, x_i);

    for (t in 1:N_t)
      eCO2_hat[t] <- totalC_t0 - sum(C_hat[t]);
    return eCO2_hat;
  }

}
data {
  real<lower=0> totalC_t0;     // initial total carbon

  real t0;                     // initial time
  int<lower=0> N_t;            // number of measurement times
  real<lower=t0> ts[N_t];      // measurement times


  real<lower=0> eCO2mean[N_t]; // measured cumulative evolved carbon
}
transformed data {
  real x_r[0];                 // no real data for ODE system
  int x_i[0];                  // no integer data for ODE system
}
parameters {
  real<lower=0> k1;            // pool 1 decomposition rate
  real<lower=0> k2;            // pool 2 decomposition rate

  real<lower=0> alpha21;       // transfer coeff from pool 2 to 1
  real<lower=0> alpha12;       // transfer coeff from pool 1 to 2

  real<lower=0,upper=1> gamma; // partitioning coefficient

  real<lower=0> sigma;         // observation std dev
}
transformed parameters {
  real eCO2_hat[N_t];
  eCO2_hat <- evolved_CO2(N_t, t0, ts, gamma, totalC_t0,
                          k1, k2, alpha21, alpha12, x_r, x_i);
}
model {
  // priors
  gamma ~ beta(10,1);         // identifies pools
  k1 ~ normal(0,1);           // weakly informative
  k2 ~ normal(0,1);
  alpha21 ~ normal(0,1);
  alpha12 ~ normal(0,1);
  sigma ~ cauchy(0,1);

  // likelihood
  for (t in 1:N_t)
    eCO2mean[t] ~ normal(eCO2_hat[t], sigma);   // normal error
}
