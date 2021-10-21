data {
  int<lower=0> N;
  vector[N] pretest;
  vector[N] setting;
  array[N] int site;
  vector[N] watched_hat;
  vector[N] y;
}
transformed data {
  vector[N] site2;
  vector[N] site3;
  vector[N] site4;
  vector[N] site5;
  matrix[N, 7] x;
  
  for (i in 1 : N) {
    site2[i] = site[i] == 2;
    site3[i] = site[i] == 3;
    site4[i] = site[i] == 4;
    site5[i] = site[i] == 5;
  }
  
  x = [watched_hat', pretest', site2', site3', site4', site5', setting']';
}
parameters {
  real alpha;
  vector[7] beta;
  real<lower=0> sigma;
}
model {
  y ~ normal_id_glm(x, alpha, beta, sigma);
}
