data {
  int<lower=0> N;
  int<lower=0> N_edges;
  int<lower=1, upper=N> node1[N_edges];  // node1[i] adjacent to node2[i]
  int<lower=1, upper=N> node2[N_edges];  // and node1[i] < node2[i]

  int<lower=0> y[N];              // count outcomes
  vector[N] x;                    // predictor
  vector<lower=0>[N] E;           // exposure

  real<lower=0> scaling_factor; //the scaling factor to make the ICAR variances approxiamtely one
}
transformed data {
  vector[N] log_E = log(E);
}
parameters {
  real beta0;                // intercept
  real beta1;                // slope

  real<lower=0> sigma;        // overall standard deviation
  real<lower=0, upper=1> rho; // proportion unstructured vs. spatially structured variance

  vector[N] theta;       // heterogeneous effects
  vector[N - 1] phi_raw; // raw spatial effects
}
transformed parameters {
  vector[N] phi;
  vector[N] bym2_re;

  phi[1:(N - 1)] = phi_raw;
  phi[N] = -sum(phi_raw);

  // NB: scaling_factor scales the spatial effect so the variance is approxiamtely 1.
  // This is NOT a magic number, and comes as data.
  // Divide by sqrt of scaling factor to properly scale precision matrix phi.
  bym2_re =  sqrt(1 - rho) * theta + sqrt(rho / scaling_factor) * phi;
}
model {
  y ~ poisson_log(log_E + beta0 + beta1 * x + bym2_re * sigma);

  // This is the prior for phi! (up to proportionality)
  target += -0.5 * dot_self(phi[node1] - phi[node2]);

  beta0 ~ normal(0.0, 5.0);
  beta1 ~ normal(0.0, 5.0);
  theta ~ normal(0.0, 1.0);
  sigma ~ normal(0,5);
  rho ~ beta(0.5, 0.5);
}
generated quantities {
  vector[N] mu = exp(log_E + beta0 + beta1 * x + bym2_re);
  real log_precision = -2.0 * log(sigma);
  real logit_rho = log(rho / (1.0 - rho));
}

