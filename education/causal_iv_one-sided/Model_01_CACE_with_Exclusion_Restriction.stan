data {
  int<lower=1> N;              // Sample size N 
  int<lower=0,upper=1> Z[N];   // Treatment assigned Z
  int<lower=0,upper=1> W[N];   // Treatment received W  
  int<lower=0,upper=1> Y[N];   // Outcome Y  
}

parameters {
  // Population probability of being a complier
  real<lower=0,upper=1> pi_c;
  
  // Probabilities for the binomial outcome distributions
  real<lower=0,upper=1> eta_c0;
  real<lower=0,upper=1> eta_c1;
  real<lower=0,upper=1> eta_n;
}  

transformed parameters {
  // Superpopulation complier average causal effect (CACE)
  // in per-1000 units
  real CACE; 
  CACE = (eta_c1 - eta_c0)*10^3;
}

model {
  // Prior for Complier probability
  // implicit prior: pi_c ~ Unif(0, 1)
  
  // Priors for outcome model parameters
  eta_c0 ~ beta(2, 2);  
  eta_c1 ~ beta(2, 2);  
  eta_n ~ beta(2, 2); 

  // Likelihood
  for(n in 1:N){
    
    // Complier (assigned to treatment)
    if(Z[n] == 1 && W[n] == 1){
      target += log(pi_c) + bernoulli_lpmf(Y[n] | eta_c1) ;
    }
    
    // Never-taker (assigned to treatment)
    else if(Z[n] == 1 && W[n] == 0){
      target +=  log(1 - pi_c) + bernoulli_lpmf(Y[n] | eta_n);
    }
    
    // Complier or Never-taker (assigned to control)
    else if(Z[n] == 0 && W[n] == 0){
      target += log_sum_exp(
        log(1 - pi_c) + bernoulli_lpmf(Y[n] | eta_n),  // Never-taker
        log(pi_c) + bernoulli_lpmf(Y[n] | eta_c0));    // Complier
    }
  }
}
