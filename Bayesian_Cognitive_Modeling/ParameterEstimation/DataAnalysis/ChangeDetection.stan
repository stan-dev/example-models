// #### Notes to Stan model #######################################################
// ## The model is not very effective.  
// ## 1) Don't change seed or lower iterations. This model converges slowly so if
// ##    you change the values, you'll need to increment iterations significantly  
// ## 2) Code is quite dissimilar to original WinBUGS model - using conditionals 
// ##    instead of step function. This will happen in further models more often.
// ##    There is a difference in what functions are efficient in BUGS and Stan.
// ################################################################################

// Change Detection
data { 
  int n;
  array[n] int t;
  array[n] real c;
}

parameters {
  array[2] real mu;
  real<lower=0> lambda;
  real<lower=0, upper=n> tau;
}

transformed parameters {
  real<lower=0> sigma;

  sigma = inv_sqrt(lambda);
}
model { 
  // Group Means
  mu ~ normal(0, sqrt(1000));
  // Common Precision
  lambda ~ gamma(0.001, 0.001);
    
  // Which Side is Time of Change Point?
  // Data Come From A Gaussian
  for (i in 1:n) {
    if ((t[i] - tau) < 0.0)
      c[i] ~ normal(mu[1], sigma);
    else 
      c[i] ~ normal(mu[2], sigma);
  }
}
