data {
  int<lower=0> N;
  int<lower=0> N_edges;
  int<lower=1, upper=N> node1[N_edges];  // node1[i] adjacent to node2[i]
  int<lower=1, upper=N> node2[N_edges];  // and node1[i] < node2[i]

  int<lower=0> y[N];              // count outcomes
  vector<lower=0>[N] E;           // exposure

  vector[N] scaling_factor; // the scaling factor to make the ICAR variances approxiamtely one
  int<lower=0> N_singletons;
  int<lower=1> singletons[N_singletons];
}
transformed data {
  vector[N] log_E = log(E);

  int N_connected = N - N_singletons;
  vector<lower=0>[N_connected] scaling_factor_connected;
  int sfc_idx = 1;
  int node_map[N];
  for (n in 1:N) {
    node_map[n] = 0;
  }
  for (n in 1:N_singletons) {
    int idx2 = singletons[n];
    node_map[idx2] = n;
  }
  for (n in 1:N) {
    if (node_map[n] == 0) {
      scaling_factor_connected[sfc_idx] = scaling_factor[n];
      sfc_idx = sfc_idx + 1;
    }
  }
}
parameters {
  real beta0;                // intercept

  real<lower=0> sigma;        // overall standard deviation
  real<lower=0, upper=1> rho; // proportion unstructured vs. spatially structured variance

  vector[N_connected] theta;       // heterogeneous effects
  vector[N_connected - 1] phi_raw; // raw spatial effects
  vector[N_singletons] singletons_re; // random effects for areas with no neighbors
}
transformed parameters {
  vector[N_connected] phi;
  vector[N_connected] convolved_re;
  vector[N] re;

  // need to sum-to-zero on a per-component basis.
  phi[1:(N_connected - 1)] = phi_raw;
  phi[N_connected] = -sum(phi_raw);
  
  // Divide by sqrt of scaling factor to properly scale precision matrix phi.
  convolved_re =  sqrt(1 - rho) * theta + sqrt(rho * inv(scaling_factor_connected)) .* phi;
  // construct RE component vector
  // nodes with neighbors have convolved RE, islands have std normal
  {
    int idx = 1;
    for (n in 1:N) {
      if (node_map[n] == 0) {
        re[n] = convolved_re[idx];
        idx = idx + 1;
      } else {
        int idx2 = node_map[n];
        re[n] = singletons_re[idx2];
      }
    }
  }
}
model {
  y ~ poisson_log(log_E + beta0 + re * sigma);

  // This is the prior for phi! (up to proportionality)
  target += -0.5 * dot_self(phi[node1] - phi[node2]);

  beta0 ~ normal(0.0, 2.5);
  theta ~ normal(0.0, 1.0);
  sigma ~ normal(0,5);
  rho ~ beta(0.5, 0.5);
  singletons_re ~ normal(0.0, 1.0);
}
generated quantities {
  real log_precision = -2.0 * log(sigma);
  real logit_rho = log(rho / (1.0 - rho));
  vector[N] eta = log_E + beta0 + re * sigma;
  vector[N] mu = exp(eta);
  vector[N] log_lik;
  int y_rep[N];
  if (max(eta) > 20) {
    print("max eta too big: ", max(eta));
    for (n in 1:N) {
      y_rep[n] = -1;
      log_lik[n] = not_a_number();
    }
  } else {
      for (n in 1:N) {
        y_rep[n] = poisson_log_rng(eta[n]);
        log_lik[n] = poisson_log_lpmf(y[n] | eta[n]);
      }
  }
}
