data {
  int<lower=0> N;
  vector[N] kid_score;
  vector[N] mom_hs;
  vector[N] mom_iq;
}
transformed data {           // interaction
  vector[N] inter = mom_hs .* mom_iq;
  matrix[N,3] x = [mom_hs', mom_iq', inter']';
}
parameters {
  real alpha;
  vector[3] beta;
  real<lower=0> sigma;
}
model {
  kid_score ~ normal_id_glm(x, alpha, beta, sigma);
}
