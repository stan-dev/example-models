// Prior and Posterior Prediction
data {
  int<lower=1> n; 
  int<lower=0> k;
} 
parameters {
  real<lower=0,upper=1> theta;
  real<lower=0,upper=1> thetaprior;
}
model {
  // Prior on Rate Theta
  theta ~ beta(1, 1);
  thetaprior ~ beta(1, 1);
  // Observed Data
  k ~ binomial(n, theta);
}
generated quantities {
  int<lower=0> postpredk;
  int<lower=0> priorpredk;
    
  // Posterior Predictive
  postpredk <- binomial_rng(n, theta);
  // Prior Predictive
  priorpredk <- binomial_rng(n, thetaprior);
}