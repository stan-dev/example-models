data {
  int<lower=0> N;
  vector<lower=0, upper=200>[N] ppvt;
  vector<lower=0, upper=1>[N] hs;
  vector<lower=0, upper=200>[N] afqt;
}
transformed data {
  matrix[N,2] x = [hs', afqt']';
}
parameters {
  real alpha;
  vector[2] beta;
  real<lower=0> sigma;
}
model {
  sigma ~ cauchy(0, 2.5);
  ppvt ~ normal_id_glm(x, alpha, beta, sigma);
}
