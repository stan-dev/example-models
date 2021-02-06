data {
  int<lower=0> N;
  int<lower=0,upper=1> switched[N];
  vector[N] dist;
  vector[N] arsenic;
  vector[N] educ;
}
transformed data {
  vector[N] c_dist100 = (dist - mean(dist)) / 100.0;
  vector[N] c_arsenic = log(arsenic) - mean(log(arsenic));
  vector[N] c_educ4 = (educ - mean(educ)) / 4.0;
  vector[N] da_inter = c_dist100 .* c_arsenic;
  vector[N] de_inter = c_dist100 .* c_educ4;
  vector[N] ae_inter = c_arsenic .* c_educ4;
  matrix[N,6] x = [c_dist100', c_arsenic', c_educ4',
                   da_inter', de_inter', ae_inter']';
}
parameters {
  real alpha;
  vector[6] beta;
}
model {
  switched ~ bernoulli_logit_glm(x, alpha, beta);
}
generated quantities {
  vector[N] pred = inv_logit(alpha + x * beta);
}
