data {
  int<lower=0> N;
  array[N] int<lower=0> y; // count outcomes
  vector<lower=0>[N] E; // exposure
  int<lower=1> K; // num covariates
  matrix[N, K] xs; // design matrix

  // neighbor graph structure
  int<lower=0, upper=N> N_components;
  array[N_components] int<lower=1, upper=N> component_sizes;
  int<lower = 0> N_edges;  // number of neighbor pairs
  array[2, N_edges] int<lower = 1, upper = (sum(component_sizes))> neighbors;  // columnwise adjacent
  vector<lower=0>[N_components] scaling_factors;
}

transformed data {
  int N_connected = sum(component_sizes);
  int N_singletons = N - N_connected;
  if (N_singletons < 0) {
    reject("Inconsistent inputs: sum(component_sizes) > N");
  }
  vector<lower=0>[N_connected] taus;
  array[N_components, 2] int node_idxs;
  array[N_components, 2] int edge_idxs;
  int node_idx = 1;
  int edge_idx = 1;
  for (n in 1:N_components) {
    taus[node_idx: component_sizes[n] + node_idx - 1]
      = rep_vector(scaling_factors[n], component_sizes[n]);
    node_idxs[n, 1] = node_idx;
    node_idxs[n, 2] = component_sizes[n] + node_idx - 1;
    node_idx += component_sizes[n];
  }
  vector[N] log_E = log(E);
  // center continuous predictors 
  vector[K] means_xs;  // column means of xs before centering
  matrix[N, K] xs_centered;  // centered version of xs
  for (k in 1:K) {
    means_xs[k] = mean(xs[, k]);
    xs_centered[, k] = xs[, k] - means_xs[k];
  }
}
parameters {
  real beta0;  // intercept
  vector[K] betas;  // covariates
  
  real<lower=0> sigma;  // random effects scale
  real<lower=0, upper=1> rho;  // proportion unstructured vs. spatially structured variance
  
  vector[N_connected] theta; // heterogeneous effects
  vector[N_singletons] singletons_re; // random effects for areas with no neighbours

  vector[N_connected] phi;
}
transformed parameters {
  vector[N] gamma;  // BYM2, per Freni-Sterrantino
  gamma[1 : N_connected] = sqrt(1 - rho) * theta + sqrt(rho * inv(taus)) .* phi;
  gamma[N_connected + 1 : N] = singletons_re;
}

model {
  y ~ poisson_log(log_E + beta0 + xs_centered * betas + gamma * sigma);
  beta0 ~ std_normal();
  betas ~ std_normal();
  theta ~ std_normal();
  singletons_re ~ std_normal();
  sigma ~ std_normal();
  rho ~ beta(0.5, 0.5);
  target += -0.5 * dot_self(phi[neighbors[1]] - phi[neighbors[2]]);  // ICAR
  for (n in 1:N_components) {   // component-wise sum-to-zero constraint
    sum(phi[node_idxs[n, 1] : node_idxs[n, 2]]) ~ normal(0,
							 0.001 * component_sizes[n]);
  }
}
generated quantities {
  real beta_intercept = beta0 - dot_product(means_xs, betas);  // adjust intercept
  array[N] int y_rep;
  {
    vector[N] eta = log_E + beta0 + xs_centered * betas + gamma * sigma;
    y_rep = max(eta) < 26 ? poisson_log_rng(eta) : rep_array(-1, N);
  }
}
