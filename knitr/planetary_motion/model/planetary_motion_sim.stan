// It's assumed that the planet orbits a star at position (0, 0),
// and that the mass of the planet is negligable next to that of
// the star (meaning the star doesn't move).

functions {
  vector ode(real t, vector y, real m, real k) {
    vector[2] q = y[1:2];
    vector[2] p = y[3:4];
    real r_cube = pow(dot_self(q), 1.5);
    
    vector[4] dydt;
    dydt[1:2] = p / m;
    dydt[3:4] = -k * q / r_cube;
    
    return dydt;
  }
}
data {
  int n;
  real<lower=0> sigma_x;
  real<lower=0> sigma_y;
}
transformed data {
  // intial state at t = 0
  real t0 = 0;
  int n_coord = 2;
  vector[n_coord] q0 = to_vector({1.0, 0.0});
  vector[n_coord] p0 = to_vector({0.0, 1.0});
  vector[n_coord * 2] y0 = append_row(q0, p0);
  
  // model parameters
  real m = 1.0;
  real k = 1.0;
  
  // exact motion
  array[n] real t;
  for (i in 1 : n) {
    t[i] = (i * 1.0) / 10;
  }
  array[2] real x_r = {m, k};
  array[n] vector[n_coord * 2] y = ode_rk45(ode, y0, t0, t, m, k);
}
parameters {
  real phi;
}
model {
  phi ~ normal(0, 1);
}
generated quantities {
  array[n, 2] real q_obs;
  
  q_obs[ : , 1] = normal_rng(y[ : , 1], sigma_x);
  q_obs[ : , 2] = normal_rng(y[ : , 2], sigma_y);
}
