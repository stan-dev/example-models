// Inferring a Common Rate
data { 
  int<lower=1> n1; 
  int<lower=1> n2; 
  int<lower=0> k1;
  int<lower=0> k2;
} 
parameters {
  real<lower=0,upper=1> theta;
} 
model {
  // Prior on Single Rate Theta
  theta ~ beta(1, 1);
  // Observed Counts
  k1 ~ binomial(n1, theta);
  k2 ~ binomial(n2, theta);
}