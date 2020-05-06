functions {
  //theta[1] = beta
  //theta[2] = gamma
  //x_i[1] = N, total population
  // y = S, I, R
  real[] sir(real t,
             real[] y,
             real[] theta,
             real[] x_r,
             int[] x_i) {
    real dydt[3];
    dydt[1] =  -  theta[1] * y[1] * y[2] / x_i[1];
    dydt[2] =    theta[1] * y[1] * y[2] / x_i[1] - theta[2] * y[2];
    dydt[3] = theta[2] * y[2];
    return dydt;

  }
}

data {
  int<lower=1> T;
  real y0[3];
  real t0;
  real ts[T];
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
    real pred_cases[T];
    real y[T,3];
    theta[1] = beta;
    theta[2] = gamma;
    y = integrate_ode_rk45(sir, y0, t0, ts, theta, x_r, x_i);
    pred_cases = neg_binomial_2_rng(col(to_matrix(y),2) + 1e-9, phi);
}
