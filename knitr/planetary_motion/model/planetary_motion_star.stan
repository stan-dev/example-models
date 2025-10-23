functions {
  vector ode(real t, vector y, real k, vector star, real m) {
    vector[2] q = y[1:2];
    vector[2] p = y[3:4];
    
    real r_cube = pow(dot_self(q - star), 1.5);
    
    vector[4] dydt;
    dydt[1:2] = p / m;
    dydt[3:4] = -k * (q - star) / r_cube;
    
    return dydt;
  }
}
data {
  int n;
  array[n, 2] real q_obs;
  array[n] real time;
  
  real<lower=0> sigma;
}
transformed data {
  real t0 = 0;
  int n_coord = 2;
  real m = 1.0;
  
  real<lower=0> sigma_x = sigma;
  real<lower=0> sigma_y = sigma;
  
  // ODE tuning parameters
  real rel_tol = 1e-6;
  real abs_tol = 1e-6;
  int max_steps = 1000;
}
parameters {
  real<lower=0> k;
  vector[n_coord] q0;
  vector[n_coord] p0;
  vector[n_coord] star;
}
transformed parameters {
  vector[n_coord * 2] y0 = append_row(q0, p0);
  array[n] vector[n_coord * 2] y = ode_rk45_tol(ode, y0, t0, time, rel_tol, abs_tol, max_steps, k, star, m);
}
model {
  k ~ normal(1, 0.001); // prior derive based on solar system 
  // (still pretty uninformative)
  
  p0[1] ~ normal(0, 1);
  p0[2] ~ lognormal(0, 1); // impose p0 to be positive.
  q0 ~ normal(0, 1);
  star ~ normal(0, 0.5);
  
  q_obs[ : , 1] ~ normal(y[ : , 1], sigma_x);
  q_obs[ : , 2] ~ normal(y[ : , 2], sigma_y);
}
generated quantities {
  array[n] real qx_pred;
  array[n] real qy_pred;
  
  qx_pred = normal_rng(y[ : , 1], sigma_x);
  qy_pred = normal_rng(y[ : , 2], sigma_y);
}
