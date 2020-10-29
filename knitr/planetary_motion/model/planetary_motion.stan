
functions {
  real[] ode (real t,
              real[] y,
              real[] theta,
              real[] x_r, int[] x_i) {
    vector[2] q = to_vector({y[1], y[2]});
    real r_cube = pow(dot_self(q), 1.5);
    vector[2] p = to_vector({y[3], y[4]});
    real m = x_r[1];
    real k = theta[1];

    vector[4] dydt;
    dydt[1:2] = p ./ m;
    dydt[3:4] = - k * q ./ r_cube;

    return to_array_1d(dydt);
  }
}

data {
  int n;
  real q_obs[n, 2];
}

transformed data {
  real t0 = 0;
  int n_coord = 2;
  real q0[n_coord] = {1.0, 0.0};
  real p0[n_coord] = {0.0, 1.0};
  real y0[n_coord * 2] = append_array(q0, p0);

  real m = 1.0;

  real t[n];
  for (i in 1:n) t[i] = i * 1.0 / 10;

  int x_i[0];
  
  real<lower = 0> sigma_x = 0.01;
  real<lower = 0> sigma_y = 0.01;
  
  // ODE tuning parameters
  real rel_tol = 1e-6;  // 1e-12;
  real abs_tol = 1e-6;  // 1e-12;
  int max_steps = 1000;  // 1e6
  
}

parameters {
  real<lower = 0> k;
  // real<lower = 0> sigma_x;
  // real<lower = 0> sigma_y;
}

transformed parameters {
  real y[n, n_coord * 2]
    = integrate_ode_bdf(ode, y0, t0, t, {k}, {m}, x_i,
                         rel_tol, abs_tol, max_steps);
}

model {
  // sigma_x ~ normal(0, 1);
  // sigma_y ~ normal(0, 1);
  k ~ normal(0, 1);
  // k ~ normal(1, 0.1);

  q_obs[, 1] ~ normal(y[, 1], sigma_x);
  q_obs[, 2] ~ normal(y[, 2], sigma_y);
}

generated quantities {
  // real q_pred[n, 2];
  real qx_pred[n];
  real qy_pred[n];

  qx_pred = normal_rng(y[, 1], sigma_x);
  qy_pred = normal_rng(y[, 2], sigma_y);
}
