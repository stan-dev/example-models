data {
  int<lower=2> K; // num topics
  int<lower=2> V; // num words
  int<lower=1> M; // num docs
  int<lower=1> N; // total word instances
  array[N] int<lower=1, upper=V> w; // word n
  array[N] int<lower=1, upper=M> doc; // doc ID for word n
  vector<lower=0>[K] alpha; // topic prior
  vector<lower=0>[V] beta; // word prior
}
parameters {
  simplex[K] theta; // topic prevalence
  array[K] simplex[V] phi; // word dist for topic k
}
model {
  array[M, K] real gamma;
  
  theta ~ dirichlet(alpha);
  for (k in 1 : K) {
    phi[k] ~ dirichlet(beta);
  }
  
  for (m in 1 : M) {
    for (k in 1 : K) {
      gamma[m, k] = categorical_lpmf(k | theta);
    }
  }
  for (n in 1 : N) {
    for (k in 1 : K) {
      gamma[doc[n], k] = gamma[doc[n], k] + categorical_lpmf(w[n] | phi[k]);
    }
  }
  for (m in 1 : M) {
    target += log_sum_exp(gamma[m]);
  }
  
  // to normalize s.t. gamma[m,k] = log Pr[Z2[m] = k|data]
  // gamma[m] <- gamma[m] - log_sum_exp(gamma[m]);
}
