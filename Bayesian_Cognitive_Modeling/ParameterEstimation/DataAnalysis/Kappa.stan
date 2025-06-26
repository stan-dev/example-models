// #### Notes to Stan model #######################################################
// ## 1) This is the first time we use simplex data type. Simplex is similar to 
// ##    vector, but with a property that the sum of all its elements is equal to 1.
// ## 2) Sampling statements for parameters alpha, beta and gamma could be removed
// ##    leading to uniform prior on (0, 1) interval which is the same as beta(1, 1)
// ## 3) Variable n was removed here. Stan doesn't need this information as
// ##    an argument for multinomial distribution. Always make sure that you know 
// ##    what arguments are required for a function / sampling statement. In many 
// ##    cases these are different from BUGS. Very useful for this are last pages 
// ##    of Stan manual
// ################################################################################

// Kappa Coefficient of Agreement
data { 
  array[4] int<lower=0> y;
}

parameters {
  // Underlying Rates

  // Rate Objective Method Decides 'one'
  real<lower=0,upper=1> alpha;

  // Rate Surrogate Method Decides 'one' When Objective Method Decides 'one'
  real<lower=0,upper=1> beta;

  // Rate Surrogate Method Decides 'zero' When Objective Method Decides 'zero'
  real<lower=0,upper=1> gamma;
} 

transformed parameters {
  simplex[4] pi;
  real xi;
  real psi;
  real kappa;

  // Probabilities For Each Count
  pi[1] = alpha * beta;
  pi[2] = alpha * (1 - beta);
  pi[3] = (1 - alpha) * (1 - gamma);
  pi[4] = (1 - alpha) * gamma;
    
  // Derived Measures   
  // Observed Rate of Agreement
  xi = alpha * beta + (1 - alpha) * gamma ;
  // Rate of Chance Agreement
  psi = (pi[1] + pi[2]) * (pi[1] + pi[3]) + (pi[2] + pi[4]) * (pi[3] + pi[4]);  
  // Chance-Corrected Agreement
  kappa = (xi - psi) / (1 - psi);
}

model {
  // Count Data     
  y ~ multinomial(pi);
}
