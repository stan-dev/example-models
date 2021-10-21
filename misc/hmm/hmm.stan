data {
  int<lower=1> K; // num categories
  int<lower=1> V; // num words
  int<lower=0> T; // num instances
  array[T] int<lower=1, upper=V> w; // words
  array[T] int<lower=1, upper=K> z; // categories
  vector<lower=0>[K] alpha; // transit prior
  vector<lower=0>[V] beta; // emit prior
}
parameters {
  array[K] simplex[K] theta; // transit probs
  array[K] simplex[V] phi; // emit probs
}
model {
  for (k in 1 : K) {
    theta[k] ~ dirichlet(alpha);
  }
  for (k in 1 : K) {
    phi[k] ~ dirichlet(beta);
  }
  for (t in 1 : T) {
    w[t] ~ categorical(phi[z[t]]);
  }
  for (t in 2 : T) {
    z[t] ~ categorical(theta[z[t - 1]]);
  }
}
