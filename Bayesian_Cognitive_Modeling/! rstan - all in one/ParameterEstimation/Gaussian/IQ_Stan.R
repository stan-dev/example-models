# clears workspace: 
rm(list=ls()) 

library(rstan)

#### Notes to Stan model #######################################################
## 1) We used implicit uniform distribution for sigma and vector mu
## 2) The second loop in this model (1:m) is not necessary here. You can try to
##    remove it, but don't forget to remove its index j. Remember vectorization!
################################################################################
model <- "
// Repeated Measures of IQ
data { 
  int<lower=1> n;
  int<lower=1> m;
  matrix[n, m] x;
}
parameters {
  vector<lower=0,upper=300>[n] mu;
  real<lower=0,upper=100> sigma;
} 
model {
  // Data Come From Gaussians With Different Means But Common Standard Deviation
  for (i in 1:n)
    for (j in 1:m)  
      x[i,j] ~ normal(mu[i], sigma);
}"

x <- matrix(c(90, 95, 100, 105, 110, 115, 150, 155, 160), 
      nrow=3, ncol=3, byrow=T) 
x

n <- nrow(x) # number of people
m <- ncol(x) # number of repeated measurements

data <- list(x=x, n=n, m=m) # to be passed on to Stan
myinits <- list(
  list(mu=rep(100, n), sigma=1))

# parameters to be monitored: 
parameters <- c("mu", "sigma")

# The following command calls Stan with specific options.
# For a detailed description type "?rstan".
samples <- stan(model_code=model,   
                data=data, 
                init=myinits,  # If not specified, gives random inits
                pars=parameters,
                iter=2000, 
                chains=1, 
                thin=1,
                # warmup = 100,  # Stands for burn-in; Default = iter/2
                # seed = 123  # Setting seed; Default is random seed
                )
# Now the values for the monitored parameters are in the "samples" object, 
# ready for inspection.
