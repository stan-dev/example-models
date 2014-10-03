# clears workspace: 
rm(list=ls()) 

library(rstan)

model <- "
// Pearson Correlation
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
  
  T[1,1] <- square(sigma[1]);
  T[1,2] <- r * sigma[1] * sigma[2];
  T[2,1] <- r * sigma[1] * sigma[2];
  T[2,2] <- square(sigma[2]);
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
E <- c(.25, .20, .25, .20, .22, .15, .09, .19, .42)

x <- matrix(cbind(N,E), nrow=9) 
n <- nrow(x)  # number of experiments

data <- list(x=x, n=n)  # to be passed on to Stan
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
                iter=5000, 
                chains=1, 
                thin=1,
                # warmup=100,  # Stands for burn-in; Default = iter/2
                # seed=123  # Setting seed; Default is random seed
)
# Now the values for the monitored parameters are in the "samples" object, 
# ready for inspection.

r <- extract(samples)$r

#============ BFs based on logspline fit ===========================
library(polspline) # this package can be installed from within R
fit.posterior <- logspline(r)

# 95% confidence interval:
x0 <- qlogspline(0.025,fit.posterior)
x1 <- qlogspline(0.975,fit.posterior)

posterior    <- dlogspline(0, fit.posterior) # this gives the pdf at point r = 0
prior        <- .5                           # height of prior at r = 0
BF10         <- prior/posterior
BF10  # 21.99

# Compare to approximation Jeffreys (1961), pp. 289-292:
BF10.J.approx <- function(n, r)
{ 
  BF10 <- 1/(((2*n-1)/pi)^.5 * (1-r^2)^(.5*(n-3)))
  return(BF10)
}
BF10.J.approx(n=length(N)-1, r=cor(N,E)) #16.20

# Compare to exact solution Jeffreys (numerical integration):
BF10.J.exact <- function(n, r)
{
  integrand <- function(rho) {((1-rho^2)^(n/2)) / ((1-rho*r)^(n-.5))}
  BF10      <- integrate(integrand, lower=-1, upper=1)$value/2
  return(BF10)
}
BF10.J.exact(n=length(N)-1, r=cor(N,E)) #21.54
#====================================================================

#============ Plot Prior and Posterior  ===========================
par(cex.main = 1.5, mar = c(5, 6, 4, 5) + 0.1, mgp = c(3.5, 1, 0), cex.lab = 1.5,
    font.lab = 2, cex.axis = 1.3, bty = "n", las=1)
xlow  <- -1
xhigh <- 1
yhigh <- 5
Nbreaks <- 80

y <- hist(r, Nbreaks, plot=F)
plot(c(y$breaks, max(y$breaks)), c(0,y$density,0), type="S", lwd=2, lty=1,
     xlim=c(xlow,xhigh), ylim=c(0,yhigh), xlab=" ", ylab=" ", main = " ", axes=F) 
axis(1)
axis(2)
mtext("Correlation r", side=1, line = 2.8, cex=2)
mtext("Density", side=2, line = 2.8, cex=2, las=0)
#now bring in log spline density estimation:
par(new=T)
plot(fit.posterior, ylim=c(0,yhigh), xlim=c(xlow,xhigh), lty=1, lwd=1, axes=F)
points(0, dlogspline(0, fit.posterior),pch=19, cex=2)
# plot the prior:
lines(c(-1,1),c(0.5,0.5),lwd=2)
points(0, .5, pch=19, cex=2)
###########################################################################
