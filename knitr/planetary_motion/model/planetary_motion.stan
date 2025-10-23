functions {
  vector ode(real t, vector y, real k, real m) {
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
  array[n, 2] real q_obs;
}
transformed data {
  real t0 = 0;
  int n_coord = 2;
  vector[n_coord] q0 = to_vector({1.0, 0.0});
  vector[n_coord] p0 = to_vector({0.0, 1.0});
  vector[n_coord * 2] y0 = append_row(q0, p0);
  
  real m = 1.0;
  
  array[n] real t;
  for (i in 1 : n) {
    t[i] = i * 1.0 / 10;
  }
  
  real<lower=0> sigma_x = 0.01;
  real<lower=0> sigma_y = 0.01;
  
  // ODE tuning parameters
  real rel_tol = 1e-6;
  real abs_tol = 1e-6;
  int max_steps = 1000;
}
parameters {
  real<lower=0> k;
}
transformed parameters {
  array[n] vector[n_coord * 2] y = ode_bdf_tol(ode, y0, t0, t, rel_tol, abs_tol, max_steps, k, m);
}
model {
  k ~ normal(0, 1);
  
  q_obs[ : , 1] ~ normal(y[ : , 1], sigma_x);
  q_obs[ : , 2] ~ normal(y[ : , 2], sigma_y);
}
generated quantities {
  array[n] real qx_pred;
  array[n] real qy_pred;
  
  qx_pred = normal_rng(y[ : , 1], sigma_x);
  qy_pred = normal_rng(y[ : , 2], sigma_y);
}
