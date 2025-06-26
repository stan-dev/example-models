# clears workspace: 
rm(list=ls()) 

library(rstan)

c <- scan("changepointdata.txt")
n <- length(c)
t <- 1:n

data <- list(c=c, n=n, t=t) # to be passed on to Stan
myinits <- list(
  list(mu=c(1, 1), lambda=1, tau=n / 2))

# parameters to be monitored:  
parameters <- c("mu", "sigma", "tau")

# The following command calls Stan with specific options.
# For a detailed description type "?rstan".
samples <- stan(file="ChangeDetection.stan",   
                data=data, 
                init=myinits,  # If not specified, gives random inits
                pars=parameters,
                iter=250, 
                chains=1, 
                thin=1,
                warmup = 150,  # Stands for burn-in; Default = iter/2
                seed = 1  # Setting seed; Default is random seed
)
# Now the values for the monitored parameters are in the "samples" object, 
# ready for inspection.

plot(samples)
mean.tau <- mean(extract(samples)$tau)
mean.mu1 <- mean(extract(samples)$mu[,1])
mean.mu2 <- mean(extract(samples)$mu[,2])

#some plotting options to make things look better:
par(cex.main = 1.5, mar = c(5, 6, 4, 5) + 0.1, mgp = c(3.5, 1, 0), cex.lab = 1.5,
    font.lab = 2, cex.axis = 1.3, bty = "n", las=1)
# the plotting: 
plot(c, type="l", main="", ylab="Values", xlab="Samples")
lines(c(1, mean.tau), c(mean.mu1,mean.mu1), lwd=2, col="red")
lines(c(mean.tau+1,length(c)), c(mean.mu2,mean.mu2), lwd=2, col="red")
