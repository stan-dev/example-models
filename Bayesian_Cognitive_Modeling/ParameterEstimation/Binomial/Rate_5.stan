// Inferring a Common Rate, With Posterior Predictive
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
generated quantities {
  int<lower=0,upper=n1> postpredk1;
  int<lower=0,upper=n2> postpredk2;
    
  // Posterior Predictive
  postpredk1 <- binomial_rng(n1, theta);
  postpredk2 <- binomial_rng(n2, theta);
}