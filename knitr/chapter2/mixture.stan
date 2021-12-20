data {
  int K;
  int N;
  array[N] real y;
  array[K] real mu;
}
parameters {
  simplex[K] theta;
  real sigma;
}
model {
  array[K] real ps;
  sigma ~ cauchy(0, 2.5);
  mu ~ normal(0, 10);
  for (n in 1 : N) {
    for (k in 1 : K) {
      ps[k] = log(theta[k]) + normal_lpdf(y[n] | mu[k], sigma);
    }
    target += log_sum_exp(ps);
  }
}
