# clears workspace: 
rm(list=ls()) 

library(rstan)

model <- "
// Optional Stopping for Bem Experiments
data { 
  int<lower=0> n;
  vector[2] x[n];
}
parameters {
  vector[2] mu;
  vector<lower=0>[2] lambda;
  real<lower=-1,upper=1> r;
} 
transformed parameters {
  vector<lower=0>[2] sigma;
  cov_matrix[2] T;

  // Reparameterization
  sigma[1] <- inv_sqrt(lambda[1]);
  sigma[2] <- inv_sqrt(lambda[2]);
  
  T[1,1] <- 1 / lambda[1];
  T[1,2] <- r * sigma[1] * sigma[2];
  T[2,1] <- r * sigma[1] * sigma[2];
  T[2,2] <- 1 / lambda[2];
}
model {
  // Priors
  mu ~ normal(0, inv_sqrt(.001));
  lambda ~ gamma(.001, .001);
  
  // Data
  x ~ multi_normal(mu, T);
}"

# Sample size N and effect size E in the Bem experiments
N <- c(100, 150, 97, 99, 100, 150, 200, 100, 50)
E <- c(0.25, 0.20, 0.25, 0.20, 0.22, 0.15, 0.09, 0.19, 0.42)

x <- matrix(cbind(N,E),nrow=9) 
n <- nrow(x) # number of experiments

data <- list(x=x, n=n) # to be passed on to Stan
myinits <- list(
    list(r=0, mu=c(0, 0), lambda=c(1, 1)))

# parameters to be monitored:	
parameters <- c("r", "mu", "sigma")

# The following command calls Stan with specific options.
# For a detailed description type "?rstan".
samples <- stan(model_code=model,   
                 data=data, 
                 init=myinits,  # If not specified, gives random inits
                 pars=parameters,
                 iter=10000, 
                 chains=1, 
                 thin=1,
                 # warmup=100,  # Stands for burn-in; Default = iter/2
                 # seed=123  # Setting seed; Default is random seed
)

r <- extract(samples)$r

#============ BFs based on logspline fit ===========================
library(polspline) # this package can be installed from within R
fit.posterior <- logspline(r)

posterior    <- dlogspline(0, fit.posterior) # this gives the pdf at point r = 0
prior        <- .5                           # height of prior at r = 0
BF10         <- prior/posterior
BF10 # 21.99

# NB. Adding the lower and upper bounds worsens the estimates:
fit.posterior <- logspline(r,lbound=-1,ubound=1)
posterior     <- dlogspline(0, fit.posterior) # this gives the pdf at point r = 0
prior         <- .5                           # height of prior at r = 0
BF10          <- prior/posterior
BF10 # 25.64

#The Fisher z-transformation:
r.F           <- atanh(r)
fit.posterior <- logspline(r.F)
posterior     <- dlogspline(atanh(0), fit.posterior) 
r.F.prior     <- atanh(runif(100000,-1,1))
fit.prior     <- logspline(r.F.prior)
prior         <- dlogspline(atanh(0), fit.prior) 
BF10          <- prior/posterior
BF10 # 20.85

#====================================================================

#============ Plot Prior and Posterior  ===========================
Nbreaks <- 80
posterior.hist <- hist(r.F, Nbreaks, plot=F) 
prior.hist <- hist(r.F.prior, Nbreaks, plot=F) 
par(cex.main = 1.5, mar = c(5, 6, 4, 5) + 0.1, mgp = c(3.5, 1, 0), cex.lab = 1.5,
    font.lab = 2, cex.axis = 1.3, bty = "n", las=1)
xlow  <- min(r.F)
xhigh <- -min(r.F)
yhigh <- 3
Nbreaks <- 80
plot(fit.posterior, ylim=c(0,yhigh), xlim=c(xlow,xhigh), lty=1, lwd=1, axes=F)
axis(1)
axis(2)
lines(c(posterior.hist$breaks, max(posterior.hist$breaks)), c(0,posterior.hist$density,0), type="S", lwd=2, lty=1) 
lines(c(prior.hist$breaks, max(prior.hist$breaks)), c(0,prior.hist$density,0), type="S", lwd=2, lty=1) 
mtext("Arctanh(r)", side=1, line = 2.8, cex=2)
mtext("Density", side=2, line = 2.8, cex=2, las=0)
points(0, dlogspline(0, fit.posterior),pch=19, cex=2)
# plot the prior:
par(new=T)
plot(fit.prior, ylim=c(0,yhigh), xlim=c(xlow,xhigh), lty=1, lwd=1, axes=F)
points(0, dlogspline(0, fit.prior), pch=19, cex=2)
###########################################################################
