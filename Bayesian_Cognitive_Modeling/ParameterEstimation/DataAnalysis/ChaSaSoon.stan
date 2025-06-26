// #### Notes to Stan model #######################################################
// ## Implementation of this model can be difficult to understand for beginners. 
// ## Therefore I suggest either not trying to understand it and look on WinBUGS
// ## version or go deep into Stan manual.
// ################################################################################
 
// ChaSaSoon Censored Data
data { 
  int<lower=0> n_fails;
  int<lower=0> n;
  int<lower=0> z_observed;
}

parameters { 
  real<lower=0.25, upper=1> theta;  // Uniform Prior on Rate Theta
}

model { 
  // Observed Data
  z_observed ~ binomial(n, theta); 
  
  // Unobserved Data
  target += n_fails * log(binomial_cdf(25 | n, theta) 
                        - binomial_cdf(14 | n, theta));
}
