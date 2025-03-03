functions {
  /**
   * Component-wise constrain sum-to-zero vectors
   *
   * @param phi unconstrained vector of zero-sum slices
   * @param idxs component start and end indices
   * @param sizes component sizes
   * @return vector phi_ozs, the vector whose slices sum to zero
   */
  vector zero_sum_components_lp(vector phi, array[ , ] int idxs, array[] int sizes) {
    vector[sum(sizes)] phi_ozs;
    int idx_phi = 1;
    int idx_ozs = 1;
    for (i in 1:size(sizes)) {
      int n = sizes[i];
      phi_ozs[idx_ozs : idx_ozs + n - 1] = 
        zero_sum_constrain_lp(segment(phi, idx_phi, n - 1));
      idx_phi += n - 1;
      idx_ozs += n;
    }
    return phi_ozs;
  }

  /**
   * Constrain sum-to-zero vector
   *
   * @param y unconstrained zero-sum parameters
   * @return vector z, the vector whose slices sum to zero
   */
  vector zero_sum_constrain_lp(vector y) {
    int N = num_elements(y);
    vector[N + 1] z = zeros_vector(N + 1);
    real sum_w = 0;
    for (ii in 1:N) {
      int i = N - ii + 1; 
      real n = i;
      real w = y[i] * inv_sqrt(n * (n + 1));
      sum_w += w;
      z[i] += sum_w;     
      z[i + 1] -= w * n;    
    }
    return z;
  }
}

data {
  int<lower=0> N;
  array[N] int<lower=0> y; // count outcomes
  vector<lower=0>[N] E; // exposure
  int<lower=1> K; // num covariates
  matrix[N, K] xs; // design matrix

  int<lower=0, upper=N> N_components;
  array[N_components] int<lower=1, upper=N> component_sizes;
  vector<lower=0>[N_components] scaling_factors;

  // neighbor graph structure
  int<lower = 0> N_edges;  // number of neighbor pairs
  array[2, N_edges] int<lower = 1, upper = (sum(component_sizes))> neighbors;  // columnwise adjacent
}

transformed data {
  int N_connected = sum(component_sizes);
  int N_singletons = N - N_connected;
  if (N_singletons < 0) {
    reject("Inconsistent inputs: sum(component_sizes) > N");
  }
  vector<lower=0>[N_connected] taus;
  array[N_components, 2] int component_idxs;
  int idx = 1;
  for (n in 1:N_components) {
    taus[idx: component_sizes[n] + idx - 1]
      = rep_vector(scaling_factors[n], component_sizes[n]);
    component_idxs[n, 1] = idx;
    component_idxs[n, 2] = component_sizes[n] + idx - 1;
    idx += component_sizes[n];
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

  // concatenation of unconstrained vectors, each size N - 1
  vector[N_connected - N_components] phi_raw;
}

transformed parameters {
  vector[N_connected] phi = zero_sum_components_lp(phi_raw, component_idxs, component_sizes);
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
}

generated quantities {
  real beta_intercept = beta0 - dot_product(means_xs, betas);  // adjust intercept
  array[N] int y_rep;
  {
    vector[N] eta = log_E + beta0 + xs_centered * betas + gamma * sigma;
    y_rep = max(eta) < 26 ? poisson_log_rng(eta) : rep_array(-1, N);
  }
}
