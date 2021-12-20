/*
 * Leuk: Cox regression 
 * URL of OpenBugs' implementation: 
 *   http://www.openbugs.net/Examples/Leuk.html
 */

data {
  int<lower=0> N;
  int<lower=0> NT;
  array[N] int<lower=0> obs_t;
  array[NT + 1] int<lower=0> t;
  array[N] int<lower=0> fail;
  array[N] real Z;
}
transformed data {
  array[N, NT] int Y;
  array[N, NT] int dN;
  real c;
  real r;
  for (i in 1 : N) {
    for (j in 1 : NT) {
      Y[i, j] = int_step(obs_t[i] - t[j] + .000000001);
      dN[i, j] = Y[i, j] * fail[i]
                 * int_step(t[j + 1] - obs_t[i] - .000000001);
    }
  }
  c = 0.001;
  r = 0.1;
}
parameters {
  real beta;
  array[NT] real<lower=0> dL0;
}
model {
  beta ~ normal(0, 1000);
  for (j in 1 : NT) {
    dL0[j] ~ gamma(r * (t[j + 1] - t[j]) * c, c);
    for (i in 1 : N) {
      if (Y[i, j] != 0) {
        target += poisson_lpmf(dN[i, j] | Y[i, j] * exp(beta * Z[i]) * dL0[j]);
      }
    }
  }
}
generated quantities {
  array[NT] real S_placebo;
  array[NT] real S_treat;
  
  for (j in 1 : NT) {
    // Survivor function = exp(-Integral{l0(u)du})^exp(beta*z)
    real s;
    s = 0;
    for (i in 1 : j) {
      s = s + dL0[i];
    }
    S_treat[j] = pow(exp(-s), exp(beta * -0.5));
    S_placebo[j] = pow(exp(-s), exp(beta * 0.5));
  }
}
