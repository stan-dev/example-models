# clears workspace: 
rm(list=ls()) 

library(rstan)

model <- "
// Two-sample Comparison of Means
data { 
  int<lower=1> n1;
  int<lower=1> n2;
  vector[n1] x;
  vector[n2] y;
}
parameters {
  real mu;
  real sigmatmp;
  real delta;
} 
transformed parameters {
  real<lower=0> sigma;
  real alpha;

  sigma <- fabs(sigmatmp);
  alpha <- delta * sigma;
}
model {
  // Delta, mu, and sigma Come From (Half) Cauchy Distribution
  mu ~ cauchy(0, 1);
  sigmatmp ~ cauchy(0, 1);
  delta ~ cauchy(0, 1);

  // Data
  x ~ normal(mu + alpha / 2, sigma);
  y ~ normal(mu - alpha / 2, sigma);
}"

x <- c(70,80,79,83,77,75,84,78,75,75,78,82,74,81,72,70,75,72,76,77)
y <- c(56,80,63,62,67,71,68,76,79,67,76,74,67,70,62,65,72,72,69,71)

n1 <- length(x)
n2 <- length(y)

# Rescale
y <- y - mean(x)
y <- y / sd(x)
x <- (x - mean(x)) / sd(x); 

data <- list(x=x, y=y, n1=n1, n2=n2) # to be passed on to Stan

myinits <- list(
  list(delta=rnorm(1,0,3), deltaprior=rnorm(1,0,3), mu=rnorm(1,0,1),
    sigmatmp=runif(1,0,5)),
  list(delta=rnorm(1,0,3), deltaprior=rnorm(1,0,3), mu=rnorm(1,0,1),
    sigmatmp=runif(1,0,5)),
  list(delta=rnorm(1,0,3), deltaprior=rnorm(1,0,3), mu=rnorm(1,0,1), 
    sigmatmp=runif(1,0,5)))

# Parameters to be monitored
parameters <- c("delta")

# The following command calls Stan with specific options.
# For a detailed description type "?rstan".
samples <- stan(model_code=model,   
                data=data, 
                init=myinits,  # If not specified, gives random inits
                pars=parameters,
                iter=20000, 
                chains=3, 
                thin=1,
                warmup=5000,  # Stands for burn-in; Default = iter/2
                # seed=123  # Setting seed; Default is random seed
)
# Now the values for the monitored parameters are in the "samples" object, 
# ready for inspection.

plot(samples)

# Collect posterior samples across all chains:
delta.posterior <- extract(samples)$delta  

#============ BFs based on logspline fit ===========================
library(polspline) # this package can be installed from within R
fit.posterior <- logspline(delta.posterior)

# 95% confidence interval:
x0 <- qlogspline(0.025,fit.posterior)
x1 <- qlogspline(0.975,fit.posterior)

posterior <- dlogspline(0, fit.posterior) # this gives the pdf at point delta = 0
prior     <- dcauchy(0)         # height of order--restricted prior at delta = 0
BF01      <- posterior/prior
BF01

#============ Plot Prior and Posterior  ===========================
par(cex.main = 1.5, mar = c(5, 6, 4, 5) + 0.1, mgp = c(3.5, 1, 0), cex.lab = 1.5,
    font.lab = 2, cex.axis = 1.3, bty = "n", las=1)
xlow  <- -3
xhigh <- 3
yhigh <- 2
Nbreaks <- 80
y       <- hist(delta.posterior, Nbreaks, plot=F)
plot(c(y$breaks, max(y$breaks)), c(0,y$density,0), type="S", lwd=2, lty=2,
     xlim=c(xlow,xhigh), ylim=c(0,yhigh), xlab=" ", ylab="Density", axes=F) 
axis(1, at = c(-4,-3,-2,-1,0,1,2,3,4), lab=c("-4","-3","-2","-1","0",
                                             "1", "2", "3", "4"))
axis(2)
mtext(expression(delta), side=1, line = 2.8, cex=2)
#now bring in log spline density estimation:
par(new=T)
plot(fit.posterior, ylim=c(0,yhigh), xlim=c(xlow,xhigh), lty=1, lwd=1, axes=F)
points(0, dlogspline(0, fit.posterior),pch=19, cex=2)
# plot the prior:
par(new=T)
plot ( function( x ) dcauchy( x, 0, 1 ), xlow, xhigh, ylim=c(0,yhigh),
       xlim=c(xlow,xhigh), lwd=2, lty=1, ylab=" ", xlab = " ", axes=F) 
axis(1, at = c(-4,-3,-2,-1,0,1,2,3,4), lab=c("-4","-3","-2","-1","0",
                                             "1", "2", "3", "4"))
axis(2)
points(0, dcauchy(0), pch=19, cex=2)
