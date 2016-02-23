# clears workspace: 
rm(list=ls()) 

library(rstan)

# to be passed on to Stan
data <- read_rdump("Rate_5.data.R")

myinits <- list(
  list(theta=.5))

# parameters to be monitored:  
parameters <- c("theta", "postpredk1", "postpredk2")

# The following command calls Stan with specific options.
# For a detailed description type "?stan".
samples <- stan(file="Rate_5.stan",   
                data=data, 
                init=myinits,  # If not specified, gives random inits
                pars=parameters,
                iter=10000, 
                chains=1, 
                thin=1,
                # warmup = 100,  # Stands for burn-in; Default = iter/2
                # seed = 123  # Setting seed; Default is random seed
)
print(samples)

theta      <- extract(samples)$theta
postpredk1 <- extract(samples)$postpredk1
postpredk2 <- extract(samples)$postpredk2

# Two-panel plot. 
layout(matrix(c(1,2),1,2))
layout.show(2)
# First, a histogram for theta.
par(cex.main = 1.5, mar = c(5, 6, 4, 5) + 0.1, mgp = c(3.5, 1, 0), cex.lab = 1.5,
    font.lab = 2, cex.axis = 1.3, bty = "n", las=1)
Nbreaks <- 80
y       <- hist(theta, Nbreaks, plot=F)
plot(c(y$breaks, max(y$breaks)), c(0,y$density,0), type="S", lwd=2, lty=1,
     xlim=c(0,1), ylim=c(0,10), xlab="Theta", ylab="Density") 
# let's plot a density estimate over this:
lines(density(theta), col="red", lwd=2)

# Second plot, the data space (predictives)
plot(data$k1, data$k2, type="p", pch=4, cex=2, lwd=2, xlab="Success Count 1",
     ylab="Success Count 2", xlim=c(-1, data$n1+1), ylim=c(-1, data$n2+1))        
nsamples <- length(theta)
sc <- 10
for (i in 0:data$n1)
{
  for (j in 0:data$n2) 
  {
    match.preds <- sum(postpredk1==i & postpredk2==j)/nsamples
    if (match.preds > 0)
    {
      points(i,j, pch=0, cex=sc*sqrt(match.preds)) 
    }
  }
}
