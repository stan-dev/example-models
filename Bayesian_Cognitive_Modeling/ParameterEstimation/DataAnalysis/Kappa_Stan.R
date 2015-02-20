# clears workspace: 
rm(list=ls()) 

library(rstan)

#### Notes to Stan model #######################################################
## 1) This is the first time we use simplex data type. Simplex is similar to 
##    vector, but with a property that sum of all it's elements is equal to 1.
## 2) Sampling statements for parameters alpha, beta and gamma could be removed
##    leading to uniform prior on (0, 1) interval which is the same as beta(1, 1)
## 3) Variable n was removed here. Stan doesn't need this information as
##    an argument for multinomial distribution. Always make sure that you know 
##    what arguments are required for a function / sampling statement. In many 
##    cases these are different from BUGS. Very useful for this are last pages 
##    of Stan manual
################################################################################
model <- "
// Kappa Coefficient of Agreement
data { 
  int<lower=0> y[4];
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
  pi[1] <- alpha * beta;
  pi[2] <- alpha * (1 - beta);
  pi[3] <- (1 - alpha) * (1 - gamma);
  pi[4] <- (1 - alpha) * gamma;
    
  // Derived Measures   
  // Rate Surrogate Method Agrees With the Objective Method
  xi <- alpha * beta + (1 - alpha) * gamma ;
  // Rate of Chance Agreement
  psi <- (pi[1] + pi[2]) * (pi[1] + pi[3]) + (pi[2] + pi[4]) * (pi[3] + pi[4]);  
  // Chance-Corrected Agreement
  kappa <- (xi - psi) / (1 - psi);
}
model {
  alpha ~ beta(1, 1);  // could be removed
  beta ~ beta(1, 1);  // could be removed
  gamma ~ beta(1, 1);  // could be removed

  // Count Data     
  y ~ multinomial(pi);
}"

# CHOOSE a data set:
# Influenza 
y <- c(14, 4, 5, 210)
# Hearing Loss 
# y <- c(20, 7, 103, 417)
# Rare Disease
# y <- c(0, 0, 13, 157)

data <- list(y=y) # to be passed on to WinBUGS
myinits <- list(
  list(alpha=.5, beta=.5, gamma=.5))

# parameters to be monitored: 
parameters <- c("kappa", "xi", "psi", "alpha", "beta", "gamma", "pi")

# The following command calls Stan with specific options.
# For a detailed description type "?rstan".
samples <- stan(model_code=model,   
                data=data, 
                init=myinits,  # If not specified, gives random inits
                pars=parameters,
                iter=4000, 
                chains=1, 
                thin=1,
                # warmup = 100,  # Stands for burn-in; Default = iter/2
                # seed = 123  # Setting seed; Default is random seed
                )
# Now the values for the monitored parameters are in the "samples" object, 
# ready for inspection.

plot(samples)

samples$mean$kappa

# Compare to Cohen's point estimate
n <- sum(y)
p0 <- (y[1]+y[4])/n
pe <- (((y[1]+y[2]) * (y[1]+y[3])) + ((y[2]+y[4]) * (y[3]+y[4]))) / n^2
kappa.Cohen <- (p0-pe) / (1-pe) 
kappa.Cohen
