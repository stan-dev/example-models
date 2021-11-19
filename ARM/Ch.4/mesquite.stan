data {
  int<lower=0> N;
  vector[N] weight;
  vector[N] diam1;
  vector[N] diam2;
  vector[N] canopy_height;
  vector[N] total_height;
  vector[N] density;
  vector[N] group;
}
transformed data {
  matrix[N, 6] x = [diam1', diam2', canopy_height', total_height', density',
                    group']';
}
parameters {
  real alpha;
  vector[6] beta;
  real<lower=0> sigma;
}
model {
  weight ~ normal_id_glm(x, alpha, beta, sigma);
}
