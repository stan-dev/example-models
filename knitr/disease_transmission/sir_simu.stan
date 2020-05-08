functions {
  real[] sir(real t, real[] y, real[] theta, 
             real[] x_r, int[] x_i) {

      real S = y[1];
      real I = y[2];
      real R = y[3];
      real N = x_i[1];
      
      real beta = theta[1];
      real gamma = theta[2];
      
      real dS_dt = -beta * I * S / N;
      real dI_dt =  beta * I * S / N - gamma * I;
      real dR_dt =  gamma * I;
      
      return {dS_dt, dI_dt, dR_dt};
  }
}
data {
  int<lower=1> n_days;
  real y0[3];
  real t0;
  real ts[n_days];
  int N;
  real<lower=0> beta;
  real<lower=0> gamma;
  real<lower=0> phi;
}
transformed data {
  real x_r[0];
  int x_i[1];
  x_i[1]=N;
}
generated quantities {
    real theta[2];
    real pred_cases[n_days];
    real y[n_days,3];
    theta[1] = beta;
    theta[2] = gamma;
    y = integrate_ode_rk45(sir, y0, t0, ts, theta, x_r, x_i);
    pred_cases = neg_binomial_2_rng(col(to_matrix(y),2) + 1e-9, phi);
}
