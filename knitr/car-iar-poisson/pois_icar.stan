data {
  int<lower=0> N;
  int<lower=0> N_edges;
  int<lower=1, upper=N> node1[N_edges];  // node1[i] adjacent to node2[i]
  int<lower=1, upper=N> node2[N_edges];  // and node1[i] < node2[i]

  int<lower=0> y[N];              // count outcomes
  vector<lower=0>[N] E;           // exposure
}
transformed data {
  vector[N] log_E = log(E);
}
parameters {
  real beta0;             // intercept
  real<lower=0> sigma;    // overall standard deviation
  vector[N] phi;         // spatial effects
}
model {
  y ~ poisson_log(log_E + beta0 + phi * sigma);
  target += -0.5 * dot_self(phi[node1] - phi[node2]);
  beta0 ~ normal(0.0, 2.5);
  sigma ~ normal(0.0, 5.0);
  // soft sum-to-zero constraint on phi)
  sum(phi) ~ normal(0, 0.001 * N);  // equivalent to mean(phi) ~ normal(0,0.001)
}
generated quantities {
  vector[N] eta = log_E + beta0 + phi * sigma;
  vector[N] mu = exp(eta);
}
