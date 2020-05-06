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
  int cases[T];
}

transformed data {
  real x_r[0];
  int x_i[1];
  x_i[1]=N;
}

parameters {
  real<lower=0> gamma;
  real<lower=0> beta;
  real<lower=0> phi_inv;
}

transformed parameters{
  real y[T,3];
  real phi = 1. / phi_inv;
  {
    real theta[2];
    theta[1] = beta;
    theta[2] = gamma;

    y = integrate_ode_rk45(sir, y0, t0, ts, theta, x_r, x_i);
  }
}


model {
  //priors
  beta ~ normal(2, 1);
  gamma ~ normal(0.4, 0.5);
  phi_inv ~ exponential(5);
  
  //sampling distribution
  //col(matrix x, int n) - The n-th column of matrix x. Here the number of infected people
  cases ~ neg_binomial_2(col(to_matrix(y),2), phi);
}

generated quantities {
  real R0 = beta/gamma;
  real recovery_time = 1/gamma;
  real pred_cases[T];
  pred_cases = neg_binomial_2_rng(col(to_matrix(y),2), phi);
}
