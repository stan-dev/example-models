data {
  int<lower=0> N; 
  vector[N] earnings;
  vector[N] height;
  vector[N] sex1;
} 
transformed data {
  vector[N] log_earnings = log(earnings);
  vector[N] male = 2 - sex1;
  vector[N] height_male_inter = height .* male;
  matrix[N,3] x = [height', male', height_male_inter']';
}
parameters {
  real alpha;
  vector[3] beta;
  real<lower=0> sigma;
}
model {
  log_earnings ~ normal_id_glm(x, alpha, beta, sigma);
}
