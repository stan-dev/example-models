# clears workspace: 
rm(list=ls())

library(rstan)

# to be passed on to Stan
data <- read_rdump("Rate_4.data.R")
# Uncomment for Trompetter Data
# data <- read_rdump("Rate_4b.data.R") 

myinits <- list(
  list(theta=.5, thetaprior=.5))

# parameters to be monitored:  
parameters <- c("theta", "thetaprior", "postpredk", "priorpredk")

# The following command calls Stan with specific options.
# For a detailed description type "?stan".
samples <- stan(file="Rate_4.stan",   
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
  
print(samples)  

######################Plots####################################################
theta   <- extract(samples)$theta
thetaprior  <- extract(samples)$thetaprior 
priorpredk  <- extract(samples)$priorpredk
postpredk   <- extract(samples)$postpredk


layout(matrix(c(1,2),2,1))
layout.show(2)
#Prior and posterior of theta
plot(density(theta, from=0, to=1), zero.line=F, axes=F, main="", xlab="",
             ylab="", xlim=c(0,1), ylim=c(0,6))
axis(1, at=c(0,0.2,0.4,0.6,0.8,1), lab=c("0","0.2","0.4","0.6","0.8","1"),
     cex.axis=0.8)
mtext("Rate", side=1, line=2.25, cex=1.2)
axis(2, at=c(0,2,4,6),cex.axis=0.8)
mtext("Density", side=2, line=2.25, cex=1.2)
lines(density(thetaprior, from=0, to=1), lty=3, col="gray")
legend(0.6,5.75, c("Prior", "Posterior"), lty=c(3,1), col=c ("grey", "black"))

#Prior and posterior predictive
mybreaks <- seq(from=-.5,to=data$n+1,by=1)
my.at    <- seq(from=0,to=data$n,by=1)
hist(postpredk,breaks=mybreaks,freq=F, right=F, ylab="", xlab="", ylim=c(0,0.3),
     main="", axes=F )
axis(1, at=my.at,lab=my.at,cex.axis=0.8)
mtext("Success Count", side=1, line=2.25, cex=1.2)
axis(2,at=c(0,0.1,0.2,0.3),lab=c("0","0.1","0.2","0.3"),cex.axis=0.8)
mtext("Mass", side=2, line=2.25, cex=1.2)
hist(priorpredk, breaks=mybreaks,freq=F,right=F,add=T, lty=3,border="grey")
legend(8,0.3, c("Prior", "Posterior"), lty=c(3,1),col=c("grey", "black"))
