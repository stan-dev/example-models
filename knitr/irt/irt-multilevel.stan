data {
  int<lower=1> I; // # items
  int<lower=1> J; // # students
  int<lower=1> K; // # schools
  int<lower=1> N; // data items
  array[N] int<lower=1, upper=I> ii; // item for n
  array[N] int<lower=1, upper=J> jj; // student for n
  array[N] int<lower=1, upper=J> kk; // school for n
  array[N] int<lower=0, upper=1> y; // correctness for n
}
parameters {
  vector[J] alpha; // ability for student j
  vector[I] beta; // difficulty for item i
  vector[K] gamma; // ability for school k
  vector<lower=0>[I] delta; // discriminativeness for item i
  real zeta; // intercept
  real<lower=0> mu_delta;
  real<lower=0, upper=20> sigma_beta;
  real<lower=0, upper=20> sigma_gamma;
  real<lower=0, upper=20> sigma_delta;
}
model {
  vector[N] log_odds;
  
  alpha ~ normal(0, 1); // priors
  beta ~ normal(0, sigma_beta);
  gamma ~ normal(0, sigma_gamma);
  delta ~ lognormal(mu_delta, sigma_delta);
  zeta ~ normal(0, 10);
  mu_delta ~ lognormal(0, 10);
  
  for (n in 1 : N) {
    // likelihood
    log_odds[n] = delta[ii[n]]
                  * (zeta + alpha[jj[n]] + gamma[kk[n]] - beta[ii[n]]);
  }
  y ~ bernoulli_logit(log_odds);
}
