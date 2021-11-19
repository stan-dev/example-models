// Fit the hyperparameters of a latent-variable Gaussian process with an
// exponentiated quadratic kernel and a Bernoulli likelihood and predict
// out-of-sample observations

data {
  int<lower=1> N1;
  array[N1] real x1;
  array[N1] int<lower=0, upper=1> z1;
  int<lower=1> N2;
  array[N2] real x2;
}
transformed data {
  real delta = 1e-9;
  int<lower=1> N = N1 + N2;
  array[N] real x;
  for (n1 in 1 : N1) {
    x[n1] = x1[n1];
  }
  for (n2 in 1 : N2) {
    x[N1 + n2] = x2[n2];
  }
}
parameters {
  real<lower=0> rho;
  real<lower=0> alpha;
  real a;
  vector[N] eta;
}
transformed parameters {
  vector[N] f;
  {
    matrix[N, N] L_K;
    matrix[N, N] K = gp_exp_quad_cov(x, alpha, rho);
    
    // diagonal elements
    for (n in 1 : N) {
      K[n, n] = K[n, n] + delta;
    }
    
    L_K = cholesky_decompose(K);
    f = L_K * eta;
  }
}
model {
  rho ~ inv_gamma(5, 5);
  alpha ~ normal(0, 1);
  a ~ normal(0, 1);
  eta ~ normal(0, 1);
  
  z1 ~ bernoulli_logit(a + f[1 : N1]);
}
generated quantities {
  array[N2] int z2;
  for (n2 in 1 : N2) {
    z2[n2] = bernoulli_logit_rng(a + f[N1 + n2]);
  }
}
