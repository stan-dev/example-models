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
  list(c=.5, r=.5, u=.5),
  list(c=.2, r=.2, u=.2),
  list(c=.8, r=.8, u=.8))

parameters <- c("c", "r", "u")  # Parameters to be monitored

# The following command calls Stan with specific options.
# For a detailed description type "?stan".
samples_1 <- stan(model_code=model,   
                data=data, 
                init=myinits,  # If not specified, gives random inits
                pars=parameters,
                iter=6000, 
                chains=3, 
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
                  iter=6000, 
                  chains=3, 
                  thin=1,
                  warmup=1000,  # Stands for burn-in; Default = iter/2
                  # seed=123  # Setting seed; Default is random seed
)

k <- c(243, 64, 65, 48)
n <- sum(k)
data <- list(k=k, n=n) # To be passed on to Stan

samples_6 <- stan(fit=samples_1,   
                  data=data, 
                  init=myinits,  # If not specified, gives random inits
                  pars=parameters,
                  iter=6000, 
                  chains=3, 
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
c6 <- extract(samples_6)$c
r6 <- extract(samples_6)$r
u6 <- extract(samples_6)$u

windows(10, 5)
layout(matrix(1:3, 1, 3, byrow=TRUE))
par(cex=1.1, mar=c(2, 2, 1, 1), mgp=c(.8, .1, 0))

plot(density(c6), xlim=c(0, 1), ylim=c(0, 15), lty="dotted",
     ylab="Probability Density", xlab="c", main="", yaxt="n", xaxt="n")
lines(density(c2), lty="dashed")
lines(density(c1))
axis(1, seq(0, 1, by=.2), tick=FALSE)
legend(0, 15.5, c("Trial1", "Trial 2", "Trial 6"),lty = c(1, 2, 3),col=c("black"),
       text.col = "black", bty = "n")

plot(density(r6), xlim=c(0, 1), ylim=c(0, 15), lty="dotted", ylab="", xlab="r",
     yaxt="n", xaxt="n", main="")
lines(density(r2), lty="dashed")
lines(density(r1))
axis(1, seq(0, 1, by=.2), tick=FALSE)

plot(density(u6), xlim=c(0, 1), ylim=c(0, 15), lty="dotted", ylab="", xlab="u",
     yaxt="n", xaxt="n", main="")
lines(density(u2), lty="dashed")
lines(density(u1))
axis(1, seq(0, 1, by=.2), tick=FALSE)

