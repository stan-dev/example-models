// Dyes: variance components model 
//  http://www.openbugs.net/Examples/Dyes.html

data {
  int BATCHES;
  int SAMPLES;
  array[BATCHES, SAMPLES] real y;
  // vector[SAMPLES] y[BATCHES]; 
}
parameters {
  real<lower=0> tau_between;
  real<lower=0> tau_within;
  real theta;
  array[BATCHES] real mu;
}
transformed parameters {
  real sigma_between;
  real sigma_within;
  sigma_between = 1 / sqrt(tau_between);
  sigma_within = 1 / sqrt(tau_within);
}
model {
  theta ~ normal(0.0, 1E5);
  tau_between ~ gamma(.001, .001);
  tau_within ~ gamma(.001, .001);
  
  mu ~ normal(theta, sigma_between);
  for (n in 1 : BATCHES) {
    y[n] ~ normal(mu[n], sigma_within);
  }
}
generated quantities {
  real sigmasq_between;
  real sigmasq_within;
  
  sigmasq_between = 1 / tau_between;
  sigmasq_within = 1 / tau_within;
}
