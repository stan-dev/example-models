// two compartment event and forcing (sim)
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
    real c = force_constant(t, {0.0,3,5,7,9}, {0.0,3,0,3,0});
    dydt[1] = -theta[1] * y[1] + c;
    dydt[2] = theta[1] * y[1] - theta[2] * y[2];
    return dydt;
  }
}
data {
  int<lower=1> T;
  int<lower=1> E;
  int<lower=1> S;
  int<lower=1> K;
  real y0[E,S];
  real t0[E];
  real ts[T];
  real theta[K];
  int start[E];
  int end[E];
}
transformed data {
  real x[0];
  int x_int[0];
}
model {
}
generated quantities {
  real y_hat[T,S];
  // 1. stop the integrator at t
  // 2. take the last value of y_hat (value at t-1) as the new y0 (initial state)
  // 3. combine y0 with the event at t
  // 4. restart the integrator
  for (i in 1:E) {
    real y_state[S];
    if (i == 1)
      y_state = y0[1,];
    else
      y_state = to_array_1d(to_vector(y0[i,]) + to_vector(y_hat[end[i-1],]));
    y_hat[start[i]:end[i],] = integrate_ode_bdf(two_comp, y_state, t0[i], ts[start[i]:end[i]], theta, x, x_int);
  }
}
