functions {
  real[] dz_dt(real t,       // time
               real[] z,     // system state
               real[] theta, // parameters
               real[] x_r,   // unused data
               int[] x_i) {
    real u = z[1];
    real v = z[2];

    real alpha = theta[1];
    real beta = theta[2];
    real gamma = theta[3];
    real delta = theta[4];

    real du_dt = (alpha - beta * v) * u;
    real dv_dt = (-gamma + delta * u) * v;

    return { du_dt, dv_dt };
  }
}
data {
  int<lower = 0> N;         // num measurements
  real ts[N];               // measurement times > 0
  real y0[2];               // initial measured population
  real<lower = 0> y[N, 2];  // measured population at measurement times
}
parameters {
  real<lower = 0> theta[4];  // theta = { alpha, beta, gamma, delta }
  real<lower = 0> z0[2];     // initial population
  real<lower = 0> sigma[2];  // measurement errors
}
transformed parameters {
  real z[N, 2]               // population for remaining years
    = integrate_ode_rk45(dz_dt, z0, 0, ts, theta,
                         rep_array(0.0, 0), rep_array(0, 0), // no data
                         1e-6, 1e-5, 1e3);  // rel tol, abs tol, max steps
  // 1e-9, 1e-, 2e4);  // control params
}
model {
  // priors, weakly informative
  sigma ~ normal(0, 0.5);
  theta[1:2] ~ normal(0, 1);
  theta[3:4] ~ normal(0, 0.2);
  z0[1] ~ normal(10, 10);
  z0[2] ~ normal(50, 50);

  // likelihood, lognormal measurement error
  y0 ~ lognormal(log(z0), sigma);
  for (k in 1:2)
    y[ , k] ~ lognormal(log(z[, k]), sigma[k]);
}
