# clears workspace: 
rm(list=ls()) 

library(rstan)

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
samples <- stan(file="Kappa.stan",   
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
