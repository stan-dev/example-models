# When you work through the code for the first time, 
# execute each command one at a time to better understand
# what it does.

# clears workspace:  
rm(list=ls()) 

library(rstan)

# to be passed on to Stan
data <- read_rdump("Rate_1.data.R") 

myinits <- list(
  list(theta=.1),  # chain 1 starting value
  list(theta=.9))  # chain 2 starting value

# parameters to be monitored:  
parameters <- c("theta")

# The following command calls Stan with specific options.
# For a detailed description type "?stan".
samples <- stan(file="Rate_1.stan",   
                data=data, 
                init=myinits,  # If not specified, gives random inits
                pars=parameters,
                iter=20000, 
                chains=2, 
                thin=1,
                # warmup=100,  # Stands for burn-in; Default = iter/2
                # seed=123  # Setting seed; Default is random seed
)
# Now the values for the monitored parameters are in the "samples" object, 
# ready for inspection.

# The commands below are useful for a quick overview:
print(samples)  # a rough summary
plot(samples)   # a visual representation

chain <- 1
as.array(samples)[1:15,chain,]  # array: element, chain, column (theta/deviance)

# Collect posterior and prior samples across all chains:
theta <- extract(samples)$theta

# Now let's plot a histogram for theta. 
# NB. Some the plots will not look good in RStudio.
# First, some options to make the plot look better:
par(cex.main = 1.5, mar = c(5, 6, 4, 5) + 0.1, mgp = c(3.5, 1, 0), cex.lab = 1.5,
    font.lab = 2, cex.axis = 1.3, bty = "n", las=1)
Nbreaks <- 80
y       <- hist(theta, Nbreaks, plot=F)
plot(c(y$breaks, max(y$breaks)), c(0,y$density,0), type="S", lwd=2, lty=1,
     xlim=c(0,1), xlab="Rate", ylab="Posterior Density") 

