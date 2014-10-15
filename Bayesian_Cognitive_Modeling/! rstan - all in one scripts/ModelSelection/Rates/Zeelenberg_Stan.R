# clears workspace: 
rm(list=ls()) 

library(rstan)


########################### CHOOSE MODEL #######################################
modelChoice <- 1  
#> 1) Simple model, direct translation from BUGS - poor performance
#> 2) Model with group distribution reparametrization using "Matt trick". This
#     model performs much better.
################################################################################

if (modelChoice == 1) {
model <- "
// Zeelenberg
data { 
  int<lower=0> ns;
  int<lower=0> nb;
  int<lower=0> nn;
  int<lower=0> sb[ns];
  int<lower=0> sn[ns];
}
parameters {
  real<lower=0> delta;
  real<lower=0> deltaprior;
  real<lower=0> mu;
  real<lower=0,upper=10> sigma;
  real<lower=0,upper=10> sigmaalpha;
  vector[ns] alpha;
  vector[ns] phin;
} 
transformed parameters {
  real<lower=0> mualpha;
  vector<lower=0,upper=1>[ns] thetab;
  vector<lower=0,upper=1>[ns] thetan;
  vector[ns] phib;

  mualpha <- delta * sigmaalpha;
  phib <- phin + alpha;

  // Probit transformation
  for (i in 1:ns) {
    thetab[i] <- Phi(phib[i]);
    thetan[i] <- Phi(phin[i]);
  }
}
model{
  // Priors
  mu ~ normal(0, 1)T[0,];
  // Priming Effect
  delta ~ normal(0, 1)T[0,];
  // Sampling from Prior Distribution for Delta
  deltaprior ~ normal(0, 1)T[0,];
  
  // Individual Parameters
  alpha ~ normal(mualpha, sigmaalpha);
  phin ~ normal(mu, sigma);
  // Data
  sb ~ binomial(nb, thetab);
  sn ~ binomial(nn, thetan);
}"
} else if (modelChoice == 2) {
model <- "
// Zeelenberg
data { 
  int<lower=0> ns;
  int<lower=0> nb;
  int<lower=0> nn;
  int<lower=0> sb[ns];
  int<lower=0> sn[ns];
}
parameters {
  real<lower=0> delta;
  real<lower=0> deltaprior;
  real<lower=0> mu;
  real<lower=0,upper=10> sigma;
  real<lower=0,upper=10> sigmaalpha;
  vector[ns] raw_p;  // Matt trick
  vector[ns] raw_a;  // Matt trick
} 
transformed parameters {
  real<lower=0> mualpha;
  vector<lower=0,upper=1>[ns] thetab;
  vector<lower=0,upper=1>[ns] thetan;
  vector[ns] phib;
  vector[ns] alpha;
  vector[ns] phin;

  mualpha <- delta * sigmaalpha;
  
  // Individual Parameters
  alpha <- mualpha + sigmaalpha * raw_a;  // Matt trick
  phin <- mu + sigma * raw_p;  // Matt trick
  
  phib <- phin + alpha;

  // Probit transformation
  for (i in 1:ns) {
    thetab[i] <- Phi(phib[i]);
    thetan[i] <- Phi(phin[i]);
  }
}
model{
  // Priors
  mu ~ normal(0, 1)T[0,];
  // Priming Effect
  delta ~ normal(0, 1)T[0,];
  // Sampling from Prior Distribution for Delta
  deltaprior ~ normal(0, 1)T[0,];
  
  raw_a ~ normal(0, 1);  // Matt trick
  raw_p ~ normal(0, 1);  // Matt trick

  // Data
  sb ~ binomial(nb, thetab);
  sn ~ binomial(nn, thetan);
}"
}

### Zeelenberg data:
# Study Both:         
sb <- c(15,11,15,14,15,18,16,16,18,16,15,13,18,12,11,13,17,18,16,11,17,18,
        12,18,18,14,21,18,17,10,11,12,16,18,17,15,19,12,21,15,16,20,15,19,
        16,16,14,18,16,19,17,11,19,18,16,16,11,19,18,12,15,18,20, 8,12,19,
        16,16,16,12,18,17,11,20)
nb <- 21

# Study Neither: 
sn <- c(15,12,14,15,13,14,10,17,13,16,16,10,15,15,10,14,17,18,19,12,19,18,
        10,18,16,13,15,20,13,15,13,14,19,19,19,18,13,12,19,16,14,17,15,16,
        15,16,13,15,14,19,12,11,17,13,18,13,13,19,18,13,13,16,18,14,14,17,
        12,12,16,14,16,18,13,13)
nn <- 21
ns <- length(sb)
         
# two-sided p-value = .03
t.test(sb, sn, alternative = c("two.sided"), paired=T)

data = list(sb=sb, sn=sn, nb=nb, nn=nn, ns=ns) # to be passed on to Stan

if (modelChoice == 1) {
  mychains <- 3
  myiterations <- 12000
  mywarmup <- 2000
  
  myinits <- list(
    list(mu=.3, sigma=.5, delta=1, sigmaalpha=1, alpha=rnorm(ns),
         phin=rnorm(ns), deltaprior=abs(rnorm(1))),
    list(mu=.5, sigma=1, delta=.01, sigmaalpha=.5, alpha=rnorm(ns),
         phin=rnorm(ns), deltaprior=abs(rnorm(1))),
    list(mu=.8, sigma=1.5, delta=.5, sigmaalpha=1.5, alpha=rnorm(ns),
         phin=rnorm(ns), deltaprior=abs(rnorm(1))))
} else if (modelChoice == 2) {
  mychains <- 2
  myiterations <- 8000
  mywarmup <- 1000
  
  myinits <- list(
    list(mu=.3, sigma=.5, delta=1, sigmaalpha=1, deltaprior=abs(rnorm(1)),
         raw_a=rnorm(ns), raw_p=rnorm(ns)),
    list(mu=.8, sigma=1.5, delta=.5, sigmaalpha=1.5, deltaprior=abs(rnorm(1)),
         raw_a=rnorm(ns), raw_p=rnorm(ns)))
}
         
         
parameters = c("delta")

# The following command calls Stan with specific options.
# For a detailed description type "?rstan".
samples <- stan(model_code=model,   
                data=data, 
                init=myinits,  # If not specified, gives random inits
                pars=parameters,
                iter=myiterations, 
                chains=mychains, 
                thin=1,
                warmup=mywarmup,  # Stands for burn-in; Default = iter/2
                # seed=123  # Setting seed; Default is random seed
)
# Now the values for the monitored parameters are in the "samples" object, 
# ready for inspection.

delta.posterior <- extract(samples)$delta      

#============ BFs based on logspline fit ===========================
library(polspline) # this package can be installed from within R
fit.posterior <- logspline(delta.posterior, lbound=0)

# 95% confidence interval:
x0 <- qlogspline(0.025,fit.posterior)
x1 <- qlogspline(0.975,fit.posterior)

posterior     <- dlogspline(0, fit.posterior) # this gives the pdf at point delta = 0
prior         <- 2 * dnorm(0)  # height of order--restricted prior at delta = 0
BF01          <- posterior/prior
# BF01 = 0.22, BF10 = 4.49

#####################################################################
### Plot Prior and Posterior under order-restriction that delta > 0 
#####################################################################

par(cex.main = 1.5, mar = c(5, 6, 4, 5) + 0.1, mgp = c(3.5, 1, 0), cex.lab = 1.5,
    font.lab = 2, cex.axis = 1.3, bty = "n", las=1)
Nbreaks <- 80
y <- hist(delta.posterior, Nbreaks, plot=F)
plot(c(y$breaks, max(y$breaks)), c(0,y$density,0), type="S", lwd=2, lty=2,
     xlim=c(0,4), ylim=c(0,1.5), xlab=" ", ylab="Density", axes=F) 
lines(c(0,0), c(0,3), col="white", lwd=4)
axis(1, at = c(0,1,2,3,4), lab=c("0", "1", "2", "3", "4"))
axis(2)
mtext(expression(delta), side=1, line = 2.8, cex=2)
#now bring in log spline density estimation:
par(new=T)
plot(fit.posterior, ylim=c(0,1.5), xlim=c(0,4), lty=1, lwd=1, axes=F)
points(0, dlogspline(0, fit.posterior),pch=19, cex=2)
# plot the prior:
par(new=T)
plot (function( x ) 2*dnorm( x, 0, 1 ), ylim=c(0,1.5), xlim=c(0,4), lwd=2,
      lty=1, ylab=" ", xlab = " ") 
points(0, 2*dnorm(0), pch=19, cex=2)

text(0.8, 1, labels = "Posterior", cex = 1.5, pos=4)
text(1.7, 0.25, labels = "Prior", cex=1.5, pos=4)
