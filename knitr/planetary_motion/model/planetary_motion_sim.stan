
// It's assumed that the planet orbits a star at position (0, 0),
// and that the mass of the planet is negligable next to that of
// the star (meaning the star doesn't move).

functions {
  real[] ode (real t,
              real[] y,
              real[] theta,
              real[] x_r, int[] x_i) {
    vector[2] q = to_vector({y[1], y[2]});
    real r_cube = pow(dot_self(q), 1.5);
    vector[2] p = to_vector({y[3], y[4]});
    real m = x_r[1];
    real k = x_r[2];

    vector[4] dydt;
    dydt[1:2] = p ./ m;
    dydt[3:4] = - k * q ./ r_cube;

    return to_array_1d(dydt);
  }
}

data {
  int n;
  real<lower = 0> sigma_x;
  real<lower = 0> sigma_y;
}

transformed data {
  // intial state at t = 0
  real t0 = 0;
  int n_coord = 2;
  real q0[n_coord] = {1.0, 0.0};
  real p0[n_coord] = {0.0, 1.0};
  real y0[n_coord * 2] = append_array(q0, p0);

  // model parameters
  real m = 1.0;
  real k = 1.0;

  // exact motion
  real t[n];
  for (i in 1:n) t[i] = (i * 1.0) / 10;
  real x_r[2] = {m, k};

  real theta[0];
  int x_i[0];

  real y[n, n_coord * 2]
    = integrate_ode_rk45(ode, y0, t0, t, theta, x_r, x_i);
}

parameters {
  real phi;
}

model {
  phi ~ normal(0, 1);
}

generated quantities {
  real q_obs[n, 2];

  q_obs[, 1] = normal_rng(y[, 1], sigma_x);
  q_obs[, 2] = normal_rng(y[, 2], sigma_y);
}
