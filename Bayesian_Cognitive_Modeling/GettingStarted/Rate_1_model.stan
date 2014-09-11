// Inferring a Rate
data { 
  int<lower=1> n; 
  int<lower=0> k;
} 
parameters {
  real<lower=0,upper=1> theta;
} 
model {
  // Prior Distribution for Rate Theta
  theta ~ beta(1, 1);
  
  // Observed Counts
  k ~ binomial(n, theta);
}