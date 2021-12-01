// GARCH(1,1)

data {
  int<lower=2> T;
  array[T] real r;
  real<lower=0> sigma1;
  int<lower=0> T_out;
}
parameters {
  real mu;
  real<lower=0> alpha0;
  real<lower=0, upper=1> alpha1;
  real<lower=0, upper=(1 - alpha1)> beta1;
}
transformed parameters {
  array[T] real<lower=0> sigma;
  sigma[1] = sigma1;
  for (t in 2 : T) {
    sigma[t] = sqrt(alpha0 + alpha1 * square(r[t - 1] - mu)
                    + beta1 * square(sigma[t - 1]));
  }
}
model {
  r ~ normal(mu, sigma);
}
generated quantities {
  array[T_out] real r_out;
  array[T_out] real sigma_out;
  sigma_out[1] = sqrt(alpha0 + alpha1 * square(r[T] - mu)
                      + beta1 * square(sigma[T]));
  r_out[1] = normal_rng(mu, sigma_out[1]);
  for (t in 2 : T_out) {
    sigma_out[t] = sqrt(alpha0 + alpha1 * square(r_out[t - 1] - mu)
                        + beta1 * square(sigma_out[t - 1]));
    r_out[t] = normal_rng(mu, sigma_out[t]);
  }
}
