data {
  int<lower=1> I;               // # items
  int<lower=1> J;               // # persons
  int<lower=1> N;               // # observations
  int<lower=1, upper=I> ii[N];  // item for n
  int<lower=1, upper=J> jj[N];  // person for n
  int<lower=0, upper=1> y[N];   // correctness for n
}
parameters {
  vector[J] theta;              // abilities
  vector[2] xi[I];              // alpha/beta pair vectors
  vector[2] mu;                 // vector for alpha/beta means
  vector<lower=0>[2] tau;       // vector for alpha/beta residual sds
  cholesky_factor_corr[2] L_Omega;
}
transformed parameters {
  vector[I] alpha;
  vector[I] beta;
  for (i in 1:I) {
    alpha[i] <- exp(xi[i,1]);
    beta[i] <- xi[i,2];
  }
}
model {
  matrix[2,2] L_Sigma;
  L_Sigma <- diag_pre_multiply(tau, L_Omega);
  for (i in 1:I)
    xi[i] ~ multi_normal_cholesky(mu, L_Sigma);
  theta ~ normal(0, 1);
  L_Omega ~ lkj_corr_cholesky(4);
  mu[1] ~ normal(0,1);
  tau[1] ~ exponential(.1);
  mu[2] ~ normal(0,5);
  tau[2] ~ exponential(.1);
  y ~ bernoulli_logit(alpha[ii] .* (theta[jj] - beta[ii]));
}
generated quantities {
  corr_matrix[2] Omega;
  Omega <- multiply_lower_tri_self_transpose(L_Omega);
}
