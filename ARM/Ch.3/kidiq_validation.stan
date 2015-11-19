data {
  int<lower=0> N;
  vector<lower=0, upper=200>[N] ppvt;
  vector<lower=0, upper=1>[N] hs;
  vector<lower=0, upper=200>[N] afqt;
}
parameters {
  vector[3] beta;
  real<lower=0> sigma;
}
model {
  sigma ~ cauchy(0, 2.5);
  ppvt ~ normal(beta[1] + beta[2] * hs + beta[3] * afqt, sigma);
}
