data {
  int<lower=1> J;
  int<lower=1> K;
  int<lower=1> N;
  array[N] int<lower=1, upper=J> jj;
  array[N] int<lower=1, upper=K> kk;
  array[N] int<lower=0, upper=1> y;
}
parameters {
  real delta;
  array[J] real alpha;
  array[K] real beta;
  array[K] real log_gamma; // log of discrimination of question k
}
model {
  alpha ~ normal(0, 1);
  beta ~ normal(0, 1);
  delta ~ normal(0.75, 1);
  log_gamma ~ normal(0, 1);
  for (n in 1 : N) {
    y[n] ~ bernoulli_logit(exp(log_gamma[kk[n]])
                           * (alpha[jj[n]] - beta[kk[n]] + delta));
  }
}
