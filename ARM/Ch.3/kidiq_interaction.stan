data {
  int<lower=0> N;
  vector<lower=0, upper=200>[N] kid_score;
  vector<lower=0, upper=200>[N] mom_iq;
  vector<lower=0, upper=1>[N] mom_hs;
}
transformed data {           // interaction
  vector[N] inter = mom_hs .* mom_iq;
  matrix[N, 3] x = [mom_hs', mom_iq', inter']';
}
parameters {
  real alpha;
  vector[3] beta;
  real<lower=0> sigma;
}
model {
  sigma ~ cauchy(0, 2.5);
  kid_score ~ normal_id_glm(x, alpha, beta, sigma);
}
