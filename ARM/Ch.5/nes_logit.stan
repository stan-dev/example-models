data {
  int<lower=0> N;
  vector[N] income;
  int<lower=0,upper=1> vote[N];
}
transformed data {
  matrix[N,1] x = [income']';
}
parameters {
  real alpha;
  vector[1] beta;
}
model {
  vote ~ bernoulli_logit_glm(x, alpha, beta);
}
