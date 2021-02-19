data {
  int<lower=0> N;
  int<lower=0,upper=1> switched[N];
  vector[N] dist;
}
transformed data {
  // rescaling
  vector[N] dist100 = dist / 100.0;   
  matrix[N,1] x = [dist100']';
}
parameters {
  real alpha;
  vector[1] beta;
}
model {
  switched ~ bernoulli_logit_glm(x, alpha, beta);
}
