# clears workspace: 
rm(list=ls()) 

library(rstan)

model <- "
// Malingering
data { 
  int<lower=1> n;
  int<lower=1> p;
  int<lower=0,upper=n> k[p];
}
parameters {
  vector<lower=0,upper=1>[2] psi;
} 
transformed parameters {
  vector[2] lp_parts[p];

  for (i in 1:p) {
    lp_parts[i, 1] <- log(.5) + binomial_log(k[i], n, psi[1]);
    lp_parts[i, 2] <- log(.5) + binomial_log(k[i], n, psi[2]);
  } 
}
model {
  // Bona Fide Group has Unknown Success Rate Above Chance
  psi[1] ~ uniform(.5, 1);
  // Malingering Group has Unknown Success Rate Below Bona Fide
  psi[2] ~ uniform(0, psi[1]);

  for (i in 1:p)
    increment_log_prob(log_sum_exp(lp_parts[i]));    
}
generated quantities {
  int<lower=0,upper=1> z[p];
  
  for (i in 1:p) {
    vector[2] prob;
    prob <- softmax(lp_parts[i]);
    z[i] <- bernoulli_rng(prob[2]);
  }
}"

k <- c(45, 45, 44, 45, 44, 45, 45, 45, 45, 45, 30,
       20, 6, 44, 44, 27, 25, 17, 14, 27, 35, 30)
p <- length(k)  # number of people
n <- 45         # number of questions

data <- list(p=p, k=k, n=n)  # To be passed on to Stan

myinits <- list(
    list(psi=c(.7, .5)))

parameters <- c("psi","z")  # Parameters to be monitored

# The following command calls Stan with specific options.
# For a detailed description type "?rstan".
samples <- stan(model_code=model,   
                data=data, 
                init=myinits,  # If not specified, gives random inits
                pars=parameters,
                iter=2000, 
                chains=1, 
                thin=1,
                # warmup = 100,  # Stands for burn-in; Default = iter/2
                # seed = 123  # Setting seed; Default is random seed
)
# Now the values for the monitored parameters are in the "samples" object, 
# ready for inspection.

print(samples, digits=3)
