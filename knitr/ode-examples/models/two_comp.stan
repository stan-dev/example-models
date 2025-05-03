// two compartment (sim)
functions {
  real[] two_comp (real t, real[] y, real[] theta, real[] x_r, int[] x_i) {
    real dydt[2];
    real a = theta[1];
    real b = theta[2];
    dydt[1] = -a*y[1];
    dydt[2] = a*y[1] - b*y[2];
    return dydt;
  }
}
data {
  int<lower=1> T;
  int<lower=1> S;
  int<lower=1> K;
  real y0[S];
  real t0;
  real ts[T];
  real theta[K];
}
transformed data {
  real x[0];
  int x_int[0];
}
model {
}
generated quantities {
  real y_hat[T,S];
  y_hat = integrate_ode_bdf(two_comp, y0, t0, ts, theta, x, x_int);
}
