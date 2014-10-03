# clears workspace: 
rm(list=ls()) 

library(rstan)

model <- "
// Multinomial Processing Tree
data { 
  int<lower=1> n;
  int<lower=0,upper=n> k[4];
}
parameters {
  real<lower=0,upper=1> c;
  real<lower=0,upper=1> r;
  real<lower=0,upper=1> u;
} 
transformed parameters {
  simplex[4] theta;
  
  // MPT Category Probabilities for Word Pairs
  theta[1] <- c * r;
  theta[2] <- (1 - c) * u ^ 2;
  theta[3] <- (1 - c) * 2 * u * (1 - u);
  theta[4] <- c * (1 - r) + (1 - c) * (1 - u) ^ 2;
}
model {
  // Priors
  c ~ beta(1, 1);  // can be removed
  r ~ beta(1, 1);  // can be removed
  u ~ beta(1, 1);  // can be removed

  // Data
  k ~ multinomial(theta);
}"

k <- c(45, 24, 97, 254)
n <- sum(k)
data <- list(k=k, n=n) # To be passed on to Stan

myinits <- list(
  list(c=.5, r=.5, u=.5))

parameters <- c("c", "r", "u")  # Parameters to be monitored

# The following command calls Stan with specific options.
# For a detailed description type "?stan".
samples_1 <- stan(model_code=model,   
                data=data, 
                init=myinits,  # If not specified, gives random inits
                pars=parameters,
                iter=11000, 
                chains=1, 
                thin=1,
                warmup=1000,  # Stands for burn-in; Default = iter/2
                # seed=123  # Setting seed; Default is random seed
)

k <- c(106, 41, 107, 166)
n <- sum(k)
data <- list(k=k, n=n) # To be passed on to Stan

samples_2 <- stan(fit=samples_1,   
                  data=data, 
                  init=myinits,  # If not specified, gives random inits
                  pars=parameters,
                  iter=11000, 
                  chains=1, 
                  thin=1,
                  warmup=1000,  # Stands for burn-in; Default = iter/2
                  # seed=123  # Setting seed; Default is random seed
)

k <- c(243, 64, 65, 48)
n <- sum(k)
data <- list(k=k, n=n) # To be passed on to Stan

samples_3 <- stan(fit=samples_1,   
                  data=data, 
                  init=myinits,  # If not specified, gives random inits
                  pars=parameters,
                  iter=11000, 
                  chains=1, 
                  thin=1,
                  warmup=1000,  # Stands for burn-in; Default = iter/2
                  # seed=123  # Setting seed; Default is random seed
)
# Now the values for the monitored parameters are in the "samples" object, 
# ready for inspection.

c1 <- extract(samples_1)$c
r1 <- extract(samples_1)$r
u1 <- extract(samples_1)$u
c2 <- extract(samples_2)$c
r2 <- extract(samples_2)$r
u2 <- extract(samples_2)$u
c3 <- extract(samples_3)$c
r3 <- extract(samples_3)$r
u3 <- extract(samples_3)$u

windows(14, 7)
layout(matrix(1:3, 1, 3, byrow=T))
plot(density(c3), xlim=c(0, 1), lty="dotted")
lines(density(c2), lty="dashed")
lines(density(c1))
plot(density(r3), xlim=c(0, 1), lty="dotted")
lines(density(r2), lty="dashed")
lines(density(r1))
plot(density(u3), xlim=c(0, 1), lty="dotted")
lines(density(u2), lty="dashed")
lines(density(u1))

