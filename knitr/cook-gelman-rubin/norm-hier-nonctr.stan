data {
  int<lower=1> J;     // number of groups
  int<lower=0> N[J];  // group size
  real y[sum(N)];    // observations
}
transformed data {
  int<lower=1, upper=J> jj[sum(N)];  // group for observation
  {
    int start = 1;
    for (j in 1:J) {
      jj[start:(start + N[j] - 1)] = rep_array(j, N[j]);
      start = start + N[j];
    }
  }
}
parameters {
  real mu;
  real<lower=0> tau;
  vector[J] alpha_std;
  real<lower=0> sigma;
}
transformed parameters {
  vector[J] alpha = mu + alpha_std * tau;
}
model {
  mu ~ normal(5, 5);
  tau ~ normal(0, 5);
  sigma ~ normal(0, 5);
  alpha_std ~ normal(0, 1);
  y ~ normal(alpha[jj], sigma);
}
