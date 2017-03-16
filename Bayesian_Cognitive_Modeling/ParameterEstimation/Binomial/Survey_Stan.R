# clears workspace: 
rm(list=ls()) 

library(rstan)

# to be passed on to Stan
data <- read_rdump("Survey.data.R")

myinits <- list(
  list(theta=.5))

# parameters to be monitored:  
parameters <- c("theta", "n")

# The following command calls Stan with specific options.
# For a detailed description type "?stan".
samples <- stan(file="Survey.stan",   
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

theta <- extract(samples)$theta
n <- extract(samples)$n

## First calculate MLE:
cc <- -Inf
ind <- 0

for (i in 1:length(n)) {
    logL <- 0
    for(j in 1:data$m) {   
        logL <- logL+lgamma(n[i]+1)-lgamma(data$k[j]+1)-lgamma(n[i]-data$k[j]+1)
        logL <- logL+data$k[j]*log(theta[i])+(n[i]-data$k[j])*log(1-theta[i])
    }
    if (logL>cc) {
        ind <- i
        cc <- logL
    }
}
# end MLE

######################Plots#####################################################
layout(matrix(c(2,0,1,3),2,2,byrow=T), width=c(2/3, 1/3), heights=c(1/3, 2/3))
xhist <- hist(n, plot=F)
yhist <- hist(theta, plot=F)
top <- max(c(xhist$counts, yhist$counts))
xrange <- c(0, data$nmax)
yrange <- c(0, 1)

par(mar=c(5, 5, 1, 1))
plot(n, theta, xlim=xrange, ylim=yrange,ylab="", xlab="")
axis(1)
mtext("Number of Surveys", side=1,line=2.25, cex=1.2)
axis(2, cex=1.2)
las=0
mtext("Rate of Return", side=2 ,line=2.25, cex=1.2)
las=1
points(mean(n), mean(theta), col="red", lwd=3, pch=4) #expectation
points(n[ind], theta[ind], col="green", lwd=3, pch=10) #Maximum Likelihood

par(mar=c(0, 4, 1, 1))
barplot(xhist$counts, axes=FALSE, ylim=c(0, top), space=0,col="lightblue")

par(mar=c(4, 0, 1, 3))
barplot(yhist$counts, axes=FALSE, xlim=c(0, top), space=0, horiz=TRUE,
        col="lightblue")
