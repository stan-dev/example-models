// Fit a Gaussian process's hyperparameters
// for squared exponential prior

data {
  int<lower=1> N1;
  real x1[N1];
  int<lower=0, upper=1> z1[N1];
  int<lower=1> N2;
  real x2[N2];
}
transformed data {
  real delta = 1e-9;
  int<lower=1> N = N1 + N2;
  real x[N];
  for (n in 1:N1) x[n] = x1[n];
  for (n in 1:N2) x[N1 + n] = x2[n];
}
parameters {
  real<lower=0> alpha;
  real<lower=0> rho;
  real a;
  vector[N] eta;
}
transformed parameters {
  vector[N] f;
  {
    matrix[N, N] L_K;
    matrix[N, N] K = cov_exp_quad(x, alpha, rho);
  
    // diagonal elements
    for (n in 1:N)
      K[n, n] = K[n, n] + delta;
    
    L_K = cholesky_decompose(K);
    f = L_K * eta;
  }
}
model {
  alpha ~ normal(0, 1);
  rho ~ gamma(4, 4);
  a ~ normal(0, 1);
  eta ~ normal(0, 1);

  z1 ~ bernoulli_logit(a + f[1:N1]);
}
generated quantities {
  int z2[N2];
  for (n in 1:N2)
    z2[n] = bernoulli_logit_rng(a + f[N1 + n]);
}
