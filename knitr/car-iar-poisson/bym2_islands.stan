data {
  int<lower=0> N;
  int<lower=0> N_edges;
  int<lower=1, upper=N> node1[N_edges];  // node1[i] adjacent to node2[i]
  int<lower=1, upper=N> node2[N_edges];  // and node1[i] < node2[i]

  int<lower=0> y[N];              // count outcomes
  vector<lower=0>[N] E;           // exposure
  //  int<lower=1> K;                 // num covariates
  //  matrix[N, K] x;                 // design matrix

  int<lower=0, upper=N> N_singletons;
  int<lower=0, upper=N> N_components;
  int<lower=1, upper=N> nodes_per_component[N_components];

  vector[N_components] scales; // per-component scaling factor
                               // makes the ICAR variances approx == 1
}
transformed data {
  vector[N] log_E = log(E);
  int N_connected = N - N_singletons;
  int N_con_comp = N_components - N_singletons;
  vector<lower=0>[N_connected] scaling_factor;   // per-node scaling factor
  int component_starts[N_components];
  int component_ends[N_components];
  int c_offset = 1;
  // calculate component offsets, set up scaling factor
  for (i in 1:N_components) {
    component_starts[i] = c_offset;
    c_offset = c_offset + nodes_per_component[i];
    component_ends[i] = c_offset - 1;
  }
  for (i in 1:N_con_comp) {
    for (j in component_starts[i]:component_ends[i]) {
      scaling_factor[j] = scales[i];
    }
  }
}
parameters {
  real beta0;                // intercept
  //  vector[K] betas;       // covariates

  real<lower=0> sigma;        // random effects scale
  real<lower=0, upper=1> rho; // proportion unstructured vs. spatially structured variance

  vector[N_connected] theta;       // heterogeneous effects
  vector[N_connected] phi;         // raw spatial effects
  vector[N_singletons] singletons_re; // random effects for areas with no neighbours
}
transformed parameters {
  vector[N] re;

  // Divide by sqrt of scaling factor to properly scale precision matrix phi.
  re[1:N_connected] = sqrt(1 - rho) * theta + sqrt(rho * inv(scaling_factor)) .* phi;
  re[(N_connected+1):N] = singletons_re;
}
model {
  y ~ poisson_log(log_E + beta0 + re * sigma);
  //  y ~ poisson_log(log_E + beta0 + x * betas + convolved_re * sigma);  // co-variates

  // This is the prior for phi! (up to proportionality)
  target += -0.5 * dot_self(phi[node1] - phi[node2]);

  beta0 ~ normal(0.0, 2.5);
  //  betas ~ normal(0.0, 2.5);
  theta ~ normal(0.0, 1.0);
  sigma ~ normal(0,5);
  rho ~ beta(0.5, 0.5);
  singletons_re ~ normal(0.0, 1.0);

 for (i in 1:N_components) {
   sum(phi[component_starts[i]:(component_ends[i])]) ~ normal(0, 0.001 * nodes_per_component[i]);
 }

}
generated quantities {
  real log_precision = -2.0 * log(sigma);
  real logit_rho = logit(rho);

  //  vector[N] eta = log_E + beta0 + x * betas + convolved_re * sigma; // co-variates
  vector[N] eta = log_E + beta0 + re * sigma;
  vector[N] mu = exp(eta);
}
