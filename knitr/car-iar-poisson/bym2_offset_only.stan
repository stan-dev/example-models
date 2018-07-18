// use for case study NYC data subset
data {
  int<lower=0> N;
  int<lower=0> N_edges;
  int<lower=1, upper=N> node1[N_edges];  // node1[i] adjacent to node2[i]
  int<lower=1, upper=N> node2[N_edges];  // and node1[i] < node2[i]

  int<lower=0> y[N];              // count outcomes
  vector<lower=0>[N] E;           // exposure

  real<lower=0> scaling_factor; // scales the variance of the spatial effects
}
transformed data {
  vector[N] log_E = log(E);
}
parameters {
  real beta0;                // intercept

  real<lower=0> sigma;        // overall standard deviation
  real<lower=0, upper=1> rho; // proportion unstructured vs. spatially structured variance

  vector[N] theta;       // heterogeneous effects
  vector[N] phi;  // spatial effects
}
transformed parameters {
  vector[N] convolved_re;
  // variance of each component should be approximately equal to 1
  convolved_re =  sqrt(1 - rho) * theta + sqrt(rho / scaling_factor) * phi;
}
model {
  y ~ poisson_log(log_E + beta0 + convolved_re * sigma);

  target += -0.5 * dot_self(phi[node1] - phi[node2]);

  beta0 ~ normal(0, 1);
  theta ~ normal(0, 1);
  sigma ~ normal(0, 1);
  rho ~ beta(0.5, 0.5);
  // soft sum-to-zero constraint on phi)
  sum(phi) ~ normal(0, 0.001 * N);  // equivalent to mean(phi) ~ normal(0,0.001)
}
generated quantities {
  real log_precision = -2.0 * log(sigma);
  real logit_rho = log(rho / (1.0 - rho));
  vector[N] eta = log_E + beta0 + convolved_re * sigma;
  vector[N] mu = exp(eta);
}
