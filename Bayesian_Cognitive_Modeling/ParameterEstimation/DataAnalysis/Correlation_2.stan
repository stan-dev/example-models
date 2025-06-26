// #### Notes to Stan model #######################################################
// ## 1) All notes from previous model Correlation_1 also apply to this model.
// ## 2) If you change sigmaerror to c(0.03, 10) as suggested in excercise 5.2.2
// ##    warning statements will be more frequent and posterior less smooth.
// ################################################################################

// Pearson Correlation With Uncertainty in Measurement
data { 
  int<lower=0> n;
  array[n] vector[2] x;
  vector[2] sigmaerror;
}

parameters {
  vector[2] mu;
  vector<lower=0>[2] lambda;
  real<lower=-1, upper=1> r;
  array[n] vector[2] y;
}

transformed parameters {
  vector<lower=0>[2] sigma;
  cov_matrix[2] T;

  // Reparameterization
  sigma[1] = inv_sqrt(lambda[1]);
  sigma[2] = inv_sqrt(lambda[2]);
  
  T[1,1] = square(sigma[1]);
  T[1,2] = r * sigma[1] * sigma[2];
  T[2,1] = r * sigma[1] * sigma[2];
  T[2,2] = square(sigma[2]);
}

model {
  // Priors
  mu ~ normal(0, sqrt(1000));
  lambda ~ gamma(0.001, 0.001);

  // Data
  y ~ multi_normal(mu, T);
  for (i in 1:n)
    x[i] ~ normal(y[i], sigmaerror);
}
