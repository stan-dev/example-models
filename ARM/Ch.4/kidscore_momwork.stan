data {
  int<lower=0> N;
  vector[N] kid_score;
  array[N] int mom_work;
}
transformed data {
  vector[N] work2;
  vector[N] work3;
  vector[N] work4;
  matrix[N, 3] x;
  
  for (i in 1 : N) {
    work2[i] = mom_work[i] == 2;
    work3[i] = mom_work[i] == 3;
    work4[i] = mom_work[i] == 4;
  }
  
  x = [work2', work3', work4']';
}
parameters {
  real alpha;
  vector[3] beta;
  real<lower=0> sigma;
}
model {
  kid_score ~ normal_id_glm(x, alpha, beta, sigma);
}
