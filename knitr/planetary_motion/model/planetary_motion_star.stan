// fixes a calculation error in planetary_motion_star.stan.

functions {
  array[] real ode(real t, array[] real y, array[] real theta,
                   array[] real x_r, array[] int x_i) {
    vector[2] q = to_vector({y[1], y[2]});
    vector[2] s = to_vector({theta[2], theta[3]});
    
    real r_cube = pow(dot_self(q - s), 1.5);
    vector[2] p = to_vector({y[3], y[4]});
    real m = x_r[1];
    real k = theta[1];
    
    vector[4] dydt;
    dydt[1 : 2] = p ./ m;
    dydt[3 : 4] = -k * (q - s) ./ r_cube;
    
    return to_array_1d(dydt);
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
  array[0] int x_i;
  
  real<lower=0> sigma_x = sigma;
  real<lower=0> sigma_y = sigma;
  
  // ODE tuning parameters
  real rel_tol = 1e-6;
  real abs_tol = 1e-6;
  int max_steps = 1000;
}
parameters {
  real<lower=0> k;
  array[n_coord] real q0;
  array[n_coord] real p0;
  array[n_coord] real star;
}
transformed parameters {
  array[n_coord * 2] real y0 = append_array(q0, p0);
  array[n_coord + 1] real theta = append_array({k}, star);
  
  array[n, n_coord * 2] real y = integrate_ode_rk45(ode, y0, t0, time, theta,
                                                    {m}, x_i, rel_tol,
                                                    abs_tol, max_steps);
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
