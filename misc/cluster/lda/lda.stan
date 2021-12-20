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
  array[M] simplex[K] theta; // topic dist for doc m
  array[K] simplex[V] phi; // word dist for topic k
}
model {
  for (m in 1 : M) {
    theta[m] ~ dirichlet(alpha);
  } // prior
  for (k in 1 : K) {
    phi[k] ~ dirichlet(beta);
  } // prior
  for (n in 1 : N) {
    array[K] real gamma;
    for (k in 1 : K) {
      gamma[k] = log(theta[doc[n], k]) + log(phi[k, w[n]]);
    }
    target += log_sum_exp(gamma); // likelihood
  }
}
