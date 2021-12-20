data {
  int<lower=0> N;
  vector<lower=0, upper=200>[N] kid_score;
  vector<lower=0, upper=1>[N] mom_hs;
}
transformed data {
  matrix[N, 1] x = [mom_hs']';
}
parameters {
  real alpha;
  vector[1] beta;
  real<lower=0> sigma;
}
model {
  sigma ~ cauchy(0, 2.5);
  kid_score ~ normal_id_glm(x, alpha, beta, sigma);
}
