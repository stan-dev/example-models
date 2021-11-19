functions {
  array[] real sir(real t, array[] real y, array[] real theta,
                   array[] real x_r, array[] int x_i) {
    real S = y[1];
    real I = y[2];
    real R = y[3];
    real N = x_i[1];
    
    real beta = theta[1];
    real gamma = theta[2];
    
    real dS_dt = -beta * I * S / N;
    real dI_dt = beta * I * S / N - gamma * I;
    real dR_dt = gamma * I;
    
    return {dS_dt, dI_dt, dR_dt};
  }
}
data {
  int<lower=1> n_days;
  array[3] real y0;
  real t0;
  array[n_days] real ts;
  int N;
  real<lower=0> beta;
  real<lower=0> gamma;
  real<lower=0> phi;
}
transformed data {
  array[0] real x_r;
  array[1] int x_i;
  x_i[1] = N;
}
generated quantities {
  array[2] real theta;
  array[n_days] real pred_cases;
  array[n_days, 3] real y;
  theta[1] = beta;
  theta[2] = gamma;
  y = integrate_ode_rk45(sir, y0, t0, ts, theta, x_r, x_i);
  pred_cases = neg_binomial_2_rng(col(to_matrix(y), 2) + 1e-9, phi);
}
