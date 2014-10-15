# clears workspace: 
rm(list=ls()) 

library(rstan)

model <- "
// Exam Scores
data { 
  int<lower=1> p;
  int<lower=0> k[p];
  int<lower=1> n;
}
transformed data {
  real psi;

  // First Group Guesses
  psi <- .5;
}
parameters {
  // Second Group Has Some Unknown Greater Rate Of Success
  real<lower=.5,upper=1> phi; 
} 
transformed parameters {
  vector[2] lp_parts[p];

  // Data Follow Binomial With Rate Given By Each Person's Group Assignment
  for (i in 1:p) {
    lp_parts[i,1] <- log(.5) + binomial_log(k[i], n, phi);
    lp_parts[i,2] <- log(.5) + binomial_log(k[i], n, psi); 
  }
}
model {
  for (i in 1:p)
    increment_log_prob(log_sum_exp(lp_parts[i]));  
}
generated quantities {
  int<lower=0,upper=1> z[p];
  
  for (i in 1:p) {
    vector[2] prob;
    prob <- softmax(lp_parts[i]);
    z[i] <- bernoulli_rng(prob[1]);
  }
}"

k <- c(21, 17, 21, 18, 22, 31, 31, 34, 34, 35, 35, 36, 39, 36, 35)
p <- length(k)  # number of people
n <- 40  # number of questions

data <- list(p=p, k=k, n=n) # to be passed on to Stan

myinits <- list(
  list(phi=.75))
  
parameters <- c("phi", "z")  # parameters to be monitored

# The following command calls Stan with specific options.
# For a detailed description type "?rstan".
samples <- stan(model_code=model,   
                data=data, 
                init=myinits,  # If not specified, gives random inits
                pars=parameters,
                iter=20000, 
                chains=1, 
                thin=1,
                # warmup = 100,  # Stands for burn-in; Default = iter/2
                # seed = 123  # Setting seed; Default is random seed
)
# Now the values for the monitored parameters are in the "samples" object, 
# ready for inspection.

print(samples, digits=3)
