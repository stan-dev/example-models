// two compartment forcing (sim)
functions {
  int max_indx (real t, real[] x) {
    int N = num_elements(x);
    int indx = 1;
    while (indx <= N && x[indx] <= t) {
      indx = indx + 1;
    }
    indx = indx - 1;
    return indx;
  }
  
  real force_constant (real t, real[] x_step, real[] y_step) {
    int indx;
    real value;
    indx = max_indx(t, x_step);
    if (indx == 0)
      value = 0;
    else
      value = y_step[indx];
    return value;
  }
  
  real[] two_comp (real t, real[] y, real[] theta, real[] x_r, int[] x_i) {
    real dydt[2];
    real a = theta[1];
    real b = theta[2];
    real c = force_constant(t, {0.0,3,5,7,9}, {0.0,3,0,3,0});
    dydt[1] = -a*y[1] + c;
    dydt[2] = a*y[1] - b*y[2];
    return dydt;
  }
}
data {
  int<lower=1> T;
  int<lower=1> K;
  int<lower=1> S;
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
  real y_hat[T,2];
  y_hat = integrate_ode_bdf(two_comp, y0, t0, ts, theta, x, x_int);
}
