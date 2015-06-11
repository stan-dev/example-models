# clears workspace: 
rm(list=ls()) 

library(rstan)

model <- "
// Geurts
data { 
  int<lower=1> nsc;
  int<lower=1> nsa;
  int<lower=0> nc[nsc];
  int<lower=0> kc[nsc];
  int<lower=0> na[nsa];
  int<lower=0> ka[nsa];
}
parameters {
  vector[nsa] phia;
  vector[nsc] phic;
  real<lower=0,upper=10> sigma;
  real mu;
  real delta;
  real deltaprior;
} 
transformed parameters {
  vector<lower=0,upper=1>[nsa] thetaa;
  vector<lower=0,upper=1>[nsc] thetac;
  real alpha;

  alpha <- sigma * delta;
  
  for (i in 1:nsc)
    thetac[i] <- Phi(phic[i]);
  for (j in 1:nsa)
    thetaa[j] <- Phi(phia[j]);
}
model {
  // Priors
  mu ~ normal(0, 1);
  delta ~ normal(0, 1);
  // Sampling from Prior Distribution for Delta
  deltaprior ~ normal(0, 1);

  phic ~ normal(mu + alpha / 2, sigma);
  phia ~ normal(mu - alpha / 2, sigma);
  kc ~ binomial(nc, thetac);
  ka ~ binomial(na, thetaa); 
}"

### Geurts data:
# Normal Controls:         
num.errors <- c(15,10,61,11,60,44,63,70,57,11,67,21,89,12,63,11,96,10,37,19,44,
                18,78,27,60,14)
nc         <- c(89,74,128,87,128,121,128,128,128,78,128,106,128,83,128,100,128,
                73,128,86,128,86,128,100,128,79)
kc         <- nc - num.errors
nsc        <- length(kc)
# ADHD:
num.errors <- c(88,50,58,17,40,18,21,50,21,69,19,29,11,76,46,36,37,72,27,92,13,
                39,53,31,49,57,17,10,12,21,39,43,49,17,39,13,68,24,21,27,48,54,
                41,75,38,76,21,41,61,24,28,21)
na         <- c(128,128,128,86,128,117,89,128,110,128,93,107,87,128,128,113,128,
                128,98,128,93,116,128,116,128,128,93,86,86,96,128,128,128,86,128,
                78,128,111,100,95,128,128,128,128,128,128,98,127,128,93,110,96)
ka         <- na - num.errors
nsa        <- length(ka)

# two-sided p-value = .72
t.test(kc / nc, ka / na, alternative=c("two.sided"), paired=F)

data <- list(nc=nc,kc=kc,nsc=nsc,na=na,ka=ka,nsa=nsa) # to be passed on to Stan

myinits <- list(
  list(mu=0, sigma=1, delta=0, phia=rnorm(nsa), phic=rnorm(nsc), deltaprior=0),
  list(mu=-.8, sigma=2, delta=-.5, phia=rnorm(nsa), phic=rnorm(nsc),
       deltaprior=-.5),
  list(mu=.8, sigma=1.5, delta=.5, phia=rnorm(nsa), phic=rnorm(nsc),
       deltaprior=.5))

parameters <- c("delta", "mu", "sigma")  # parameters to be monitored

# The following command calls Stan with specific options.
# For a detailed description type "?rstan".
samples <- stan(model_code=model,   
                data=data, 
                init=myinits,  # If not specified, gives random inits
                pars=parameters,
                iter=20000, 
                chains=3, 
                thin=1,
                # warmup=100,  # Stands for burn-in; Default = iter/2
                # seed=123  # Setting seed; Default is random seed
)
# Now the values for the monitored parameters are in the "samples" object,
# ready for inspection, e.g.;
traceplot(samples, inc_warmup=FALSE)

# Collect posterior samples across all chains:
delta.posterior <- extract(samples)$delta      

#============ Unrestricted analysis      ===========================
#============ BFs based on logspline fit ===========================
library(polspline) # this package can be installed from within R
fit.posterior <- logspline(delta.posterior)

# 95% confidence interval:
(x0 <- qlogspline(0.025,fit.posterior))
(x1 <- qlogspline(0.975,fit.posterior))

posterior <- dlogspline(0, fit.posterior) # this gives the pdf at point delta = 0
prior     <- dnorm(0)            # height of order-restricted prior at delta = 0
(BF01      <- posterior/prior)
# BF01 = 3.96, BF10 = 0.25

#####################################################################
### Plot Prior and Posterior
#####################################################################
par(cex.main = 1.5, mar = c(5, 6, 4, 5) + 0.1, mgp = c(3.5, 1, 0), cex.lab = 1.5,
    font.lab = 2, cex.axis = 1.3, bty = "n", las=1)
Nbreaks <- 80
y <- hist(delta.posterior, Nbreaks, plot=F)
plot(c(y$breaks, max(y$breaks)), c(0,y$density,0), type="S", lwd=2, lty=2,
     xlim=c(-3,3), ylim=c(0,2), xlab=" ", ylab="Density", axes=F) 
axis(1, at = c(-3,-2,-1,0,1,2,3), lab=c("-3","-2", "-1", "0", "1", "2", "3"))
axis(2)
mtext(expression(delta), side=1, line = 2.8, cex=2)
#now bring in log spline density estimation:
par(new=T)
plot(fit.posterior, ylim=c(0,2), xlim=c(-3,3), lty=1, lwd=1, axes=F)
points(0, dlogspline(0, fit.posterior),pch=19, cex=2)
# plot the prior:
curve(dnorm( x, 0, 1 ), -3, 3, xlim=c(-3,3), ylim=c(0,2), lwd=2,
      lty=1, ylab=" ", xlab = " ", add = TRUE) 
points(0, dnorm(0), pch=19, cex=2)

text(0.3, 1.5, labels = "Posterior", cex = 1.5, pos=4)
text(1, 0.3, labels = "Prior", cex=1.5, pos=4)

