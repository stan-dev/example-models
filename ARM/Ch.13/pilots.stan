data {
  int<lower=0> N; 
  int<lower=0> n_treaments; 
  int<lower=0> n_airports; 
  int<lower=1,upper=n_treatments> treatment_id[N];
  int<lower=1,upper=n_airports> airport_id[N];
  vector[N] y;
} 
parameters {
  vector[n_treatments] gamma_std;
  vector[n_airports] delta_std;
  real mu;
  real<lower=0> sigma_gamma;
  real<lower=0> sigma_delta;
  real<lower=0> sigma_y;
}
transformed parameters {
  vector[n_treatments] gamma;
  vector[n_airports] delta;
  vector[N] y_hat;

  gamma <- sigma_gamma * gamma_std;
  delta <- sigma_delta * delta_std;

  for (i in 1:N)
    y_hat[i] <- mu + gamma[treatment_id[i]] + delta[airport_id[i]];
} 
model {
  sigma_gamma ~ cauchy(1,5);
  sigma_delta ~ cauchy(1,5);

  mu ~ normal(0, 10);
  gamma_std ~ normal(0, 1);
  delta_std ~ normal(0, 1);
  y ~ normal(y_hat, sigma_y);
}