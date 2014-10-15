# clears workspace: 
rm(list=ls()) 

library(rstan)

#### Notes to Stan model #######################################################
## Implementation of this model can be difficult to understand for beginners. 
## Therefore I suggest either not trying to understand it and look on WinBUGS
## version or go deep into Stan manual.
################################################################################
model <- "
# ChaSaSoon Censored Data
data { 
  int<lower=0> nfails;
  int<lower=0> n;
  int<lower=0> z_observed;
} 
parameters { 
  real<lower=.25,upper=1> theta;  // Uniform Prior on Rate Theta
} 
model { 
  // Observed Data
  z_observed ~ binomial(n, theta); 
  
  // Unobserved Data
  increment_log_prob(nfails * log(binomial_cdf(25, n, theta) 
                                  - binomial_cdf(14, n, theta)));
}"
 
nfails <- 949  
n <- 50  # Number of questions  
z_observed <- 30  # Score on the successful trial

data <- list(nfails=nfails, n=n, z_observed=z_observed) # to be passed on to Stan

myinits <- list(
  list(theta=.5))

parameters <- c("theta")

# The following command calls Stan with specific options.
# For a detailed description type "?rstan".
samples <- stan(model_code=model,   
                data=data, 
                init=myinits,  # If not specified, gives random inits
                pars=parameters,
                iter=10000, 
                chains=1, 
                thin=1,
                # warmup = 100,  # Stands for burn-in; Default = iter/2
                # seed = 123  # Setting seed; Default is random seed
)
# Now the values for theta are in the "samples" object, ready for inspection.

# Collect all samples in "theta":
theta <- extract(samples)$theta 

# Plot the posterior for theta:
par(cex.main = 1.5, mar = c(5, 6, 4, 5) + 0.1, mgp = c(3.5, 1, 0), cex.lab = 1.5,
    font.lab = 2, cex.axis = 1.3, bty = "n", las=1)
plot(density(theta), ylim=c(0,200), xlim=c(0.2,0.5), lty=1, lwd=2, 
     axes=F, main=" ", xlab=" ", ylab="Posterior Density")
axis(1, at = c(0.20, 0.30, 0.40, 0.50), lab=c("0.20", "0.30", "0.40", "0.50"))
axis(2)
mtext(expression(paste(theta," for Cha Sa-soon")), side=1, line = 2.8, cex=2)

# plot 95% confidence interval
x0 <- quantile(theta,p=c(.025,.975))[[1]]
x1 <- quantile(theta,p=c(.025,.975))[[2]]
arrows(x0, 150, x1, 150, length = 0.05, angle = 90, code = 3, lwd=2)
text((x0+x1)/2, 170, labels="95%", cex = 1.5) 
text(x0-.03, 150, labels=as.character(round(x0,2)), cex = 1.5) 
text(x0+.04, 150, labels=as.character(round(x1,2)), cex = 1.5) 
