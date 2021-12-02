data {
  int<lower=0> N; // number of data points
  int<lower=1> D; // number of dimensions
  int<lower=1> K; // number of clusters
  array[N] vector[D] y; // observations
}
transformed data {
  real<upper=0> neg_log_K;
  neg_log_K = -log(K);
}
parameters {
  array[K] vector[D] mu; // cluster means
}
transformed parameters {
  array[N, K] real<upper=0> soft_z; // log unnormalized cluster assigns
  for (n in 1 : N) {
    for (k in 1 : K) {
      soft_z[n, k] = neg_log_K - 0.5 * dot_self(mu[k] - y[n]);
    }
  }
}
model {
  for (k in 1 : K) {
    mu[k] ~ normal(0, 1);
  } // prior
  for (n in 1 : N) {
    target += log_sum_exp(soft_z[n]);
  } // likelihood
}
