// LSAT: item response
// http://www.openbugs.net/Examples/Lsat.html

data {
  int<lower=0> N; // 1000, number of students
  int<lower=0> R; // 32, number of patterns of results: 2^T
  int<lower=0> T; // 5, number of questions
  array[R] int<lower=0> culm;
  array[R, T] int<lower=0> response;
}
transformed data {
  array[T, N] int r;
  vector[N] ones;
  
  for (j in 1 : culm[1]) {
    for (k in 1 : T) {
      r[k, j] = response[1, k];
    }
  }
  for (i in 2 : R) {
    for (j in (culm[i - 1] + 1) : culm[i]) {
      for (k in 1 : T) {
        r[k, j] = response[i, k];
      }
    }
  }
  for (i in 1 : N) {
    ones[i] = 1.0;
  }
}
parameters {
  array[T] real alpha;
  vector[N] theta;
  real<lower=0> beta;
}
model {
  alpha ~ normal(0, 100.);
  theta ~ normal(0, 1);
  beta ~ normal(0.0, 100.);
  for (k in 1 : T) {
    r[k] ~ bernoulli_logit(beta * theta - alpha[k] * ones);
  }
}
generated quantities {
  real mean_alpha;
  array[T] real a;
  mean_alpha = mean(alpha);
  for (t in 1 : T) {
    a[t] = alpha[t] - mean_alpha;
  }
}
