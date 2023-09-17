data {
  int<lower=1> N;               // Sample size N 
  int<lower=0, upper=1> Z[N];   // Treatment assigned Z
  int<lower=0, upper=1> W[N];   // Treatment received W  
  int<lower=0, upper=1> Y[N];   // Outcome Y  
}

parameters {
  // Population probability of being a complier
  real<lower=0, upper=1> pi_c;
  
  // Probabilities for the binomial outcome distributions
  real<lower=0, upper=1> eta_c0;
  real<lower=0, upper=1> eta_c1;
  real<lower=0, upper=1> eta_n;
}  

transformed parameters {
  // Superpopulation complier average causal effect (CACE)
  // in per-1000 units
  real CACE = (eta_c1 - eta_c0) * 10^3;
}

model {
  // Define local variables for efficiency
  real log_pi_c = log(pi_c);
  real log1m_pi_c = log1m(pi_c);
  
  // Prior for Complier probability
  // implicit prior: pi_c ~ Unif(0, 1)
  
  // Priors for outcome model parameters
  eta_c0 ~ beta(2, 2);  
  eta_c1 ~ beta(2, 2);  
  eta_n ~ beta(2, 2); 

  // Likelihood
  for(n in 1:N){
    
    // Complier (assigned to treatment)
    if (Z[n] == 1 && W[n] == 1){
      target += log_pi_c + bernoulli_lpmf(Y[n] | eta_c1) ;
    }
    
    // Never-taker (assigned to treatment)
    else if (Z[n] == 1 && W[n] == 0){
      target +=  log1m_pi_c + bernoulli_lpmf(Y[n] | eta_n);
    }
    
    // Complier or Never-taker (assigned to control)
    else if (Z[n] == 0 && W[n] == 0){
      target += log_mix(
        pi_c,                           // Complier probability
        bernoulli_lpmf(Y[n] | eta_c0),  // Complier
        bernoulli_lpmf(Y[n] | eta_n)    // Never-taker
      );  
    }
  }
}
