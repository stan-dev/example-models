# clears workspace: 
rm(list=ls()) 

library(rstan)

########################### CHOOSE MODEL #######################################
modelChoice <- 1  
#> 1) Simple model with a restriction in parameter declaration (it did the trick)
#> 2) Model replicating approach from the book 
################################################################################

if (modelChoice == 1) {
model <- "
  // Pledgers, Order Restricted Rates
  data { 
    int<lower=0> n1;
    int<lower=0> n2;
    int<lower=0,upper=n1> s1;
    int<lower=0,upper=n2> s2;
  }
  parameters {
    real<lower=0,upper=1> theta2;
    real<lower=0,upper=theta2> theta1;
    real<lower=0,upper=1> theta2prior;
    real<lower=0,upper=theta2prior> theta1prior;
  }
  transformed parameters {
    real<lower=-1,upper=1> delta;
    real<lower=-1,upper=1> deltaprior;

    // Computation of delta and deltaprior 
    delta <- theta1 - theta2;
    deltaprior <- theta1prior - theta2prior;
  }
  model {
    // Prior Sampling
    theta1 ~ beta(1, 1);
    theta2 ~ beta(1, 1);
    theta1prior ~ beta(1, 1);
    theta2prior ~ beta(1, 1);

    // Data
    s1 ~ binomial(n1, theta1);
    s2 ~ binomial(n2, theta2);
  }"
} else if (modelChoice == 2) {
model <- "
  // Pledgers, Order Restricted Rates
  data { 
    int<lower=0> n1;
    int<lower=0> n2;
    int<lower=0,upper=n1> s1;
    int<lower=0,upper=n2> s2;
  }
  transformed data {
    cov_matrix[2] TI;
    vector[2] mu;
    real<lower=0> angle;

    // Constants
    angle <- 45 * pi() / 180;
    TI[1,1] <- 1.0;
    TI[1,2] <- 0.0;
    TI[2,1] <- 0.0;
    TI[2,2] <- 1.0;
    mu[1] <- 0.0;
    mu[2] <- 0.0;
  }
  parameters {
    vector[2] thetap;
    vector[2] thetapprior;
  }
  transformed parameters {
    real<lower=-1,upper=1> delta;
    real<lower=0,upper=1> theta1;
    real<lower=0,upper=1> theta2;
    real<lower=-1,upper=1> deltaprior;
    real<lower=0,upper=1> theta1prior;
    real<lower=0,upper=1> theta2prior;
    
    theta1 <- Phi((cos(angle) * thetap[1]) - (sin(angle) * fabs(thetap[2])));
    theta2 <- Phi((sin(angle) * thetap[1]) + (cos(angle) * fabs(thetap[2])));
    theta1prior <- Phi(cos(angle) * thetapprior[1]
                       - sin(angle) * fabs(thetapprior[2]));
    theta2prior <- Phi(sin(angle) * thetapprior[1] 
                       + cos(angle) * fabs(thetapprior[2]));
    // Difference
    delta <- theta1 - theta2;
    deltaprior <- theta1prior - theta2prior;
  }
  model {
    // Order Constrained Rates
    thetap ~ multi_normal(mu, TI);
    // Prior Sampling
    thetapprior ~ multi_normal(mu, TI);
    // Data
    s1 ~ binomial(n1, theta1);
    s2 ~ binomial(n2, theta2);
  }"
}

# Pledger data:
s1 = 424
s2 = 5416
n1 = 777
n2 = 9072

data = list(s1=s1, s2=s2, n1=n1, n2=n2) # to be passed on to Stan

if (modelChoice == 1) { 
  myinits <- list(
    list(theta1=.8, theta2=.9, theta1prior=.8, theta2prior=.9),
    list(theta1=.5, theta2=.7, theta1prior=.5, theta2prior=.7),
    list(theta1=.1, theta2=.3, theta1prior=.1, theta2prior=.3))
} else if (modelChoice == 2) { 
  myinits <- list(
    list(thetap=c(-.8, -.4), thetapprior=c(-.8, -.4)),
    list(thetap=c(-.5, -.25), thetapprior=c(-.5, -.25)),
    list(thetap=c(-.2, -.1), thetapprior=c(-.2, -.1)))
}

parameters <- c("delta", "deltaprior")

# The following command calls Stan with specific options.
# For a detailed description type "?rstan".
samples <- stan(model_code=model,   
                data=data, 
                init=myinits,  # If not specified, gives random inits
                pars=parameters,
                iter=60000, 
                chains=3, 
                thin=1,
                warmup=10000,  # Stands for burn-in; Default = iter/2
                # seed=123  # Setting seed; Default is random seed
)
# Now the values for the monitored parameters are in the "samples" object, 
# ready for inspection.

######################################################
# Order-restriction. H2: delta < 0
######################################################
# Collect posterior samples across all chains:
delta.posterior  <- extract(samples)$delta
delta.prior      <- extract(samples)$deltaprior

#============ BFs based on logspline fit ===========================
library(polspline) # this package can be installed from within R

fit.posterior <- logspline(delta.posterior, lbound=-1, ubound=0)
posterior     <- dlogspline(0, fit.posterior) # this gives the pdf at point delta = 0
BF02          <- posterior/2 # 0.26, BF20 = 3.78

#======= Plot Order-Restricted Prior and Posterior ======================
par(cex.main = 1.5, mar = c(5, 6, 4, 5) + 0.1, mgp = c(3.5, 1, 0), cex.lab = 1.5,
    font.lab = 2, cex.axis = 1.3, bty = "n", las=1)
Nbreaks <- 20
y <- hist(delta.prior, Nbreaks, plot=F)
plot(c(y$breaks, max(y$breaks)), c(0,y$density,0), type="S", lwd=2, lty=1,
     xlim=c(-1,0), ylim=c(0,25), xlab=" ", ylab="Density", axes=F) 
lines(c(0,0), c(0,2), col="white", lwd=2)
axis(1, at = c(-1, -0.5, 0), lab=c("-1", "-0.5", "0"))
axis(2)
mtext(expression(delta), side=1, line = 2.8, cex=2)
par(new=T)
x <- hist(delta.posterior, Nbreaks, plot=F)
plot(c(x$breaks, max(x$breaks)), c(0,x$density,0), type="S", lwd=2, lty=2,
     xlim=c(-1,0), ylim=c(0,25), xlab=" ", ylab="Density", main ="Full Scale",
     axes=F) 
axis(1, at = c(-1, -0.5, 0), lab=c("-1", "-0.5", "0"))
axis(2)
#now bring in log spline density estimation:
par(new=T)
# plot the prior:
lines(c(-1,0),c(0,2), lty=1, lwd=1)
par(new=T)
plot(fit.posterior, ylim=c(0,25), xlim=c(-1,0), lty=1, lwd=1, axes=F)

text(-0.42, 20, labels = "Posterior", cex = 1.5, pos=4)
text(-0.5, 2.5, labels = "Prior", cex=1.5, pos=4)

######## Second plot, zoom in:
par(cex.main = 1.5, mar = c(5, 6, 4, 5) + 0.1, mgp = c(3.5, 1, 0), cex.lab = 1.5,
    font.lab = 2, cex.axis = 1.3, bty = "n", las=1)
xmin <- -0.05
xmax <- 0
ymax <- 5
plot(0,0, ylim=c(0,ymax), xlim=c(xmin,xmax), lwd=2, lty=3, ylab="Density", 
     xlab=" ", main="Zoomed in", axes=F, col="white") 
#white makes this invisible
axis(1, at = c(xmin, -0.025, xmax), lab=c(paste(xmin), -0.025, paste(xmax)))
axis(2)
mtext(expression(delta), side=1, line = 2.8, cex=2)
par(new=T)
plot(fit.posterior, ylim=c(0,ymax), xlim=c(xmin,xmax), lty=2, lwd=2, axes=F)
lines(c(-1,0),c(0,2), lty=1, lwd=2)
points(0, 2, pch=19, cex=2)
points(0, dlogspline(0, fit.posterior),pch=19, cex=2)

text(-0.04, 1.7, labels = "Prior", cex=1.5, pos=4)
text(-0.037, 4, labels = "Posterior", cex = 1.5, pos=4)
