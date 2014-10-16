# clears workspace: 
rm(list=ls()) 

library(rstan)

#### Notes to Stan model #######################################################
## 1) If parameter's prior distribution is not specified, Stan will assume that
##    you want it to be distributed uniformly with boundaries given by variable 
##    constraints. Here constrains <lower=0,upper=10> give uniform (0, 100)
## 2) In Stan, most of the sampling statements can be vectorized. In this example
##    you can see it in the statement for vector x. Instead of using for loop for  
##    each element of the vector, we can simple write it as above. This saves code, 
##    speeds up computation. For more information read Vectorization chapter in
##    the Stan manual (p.231 in version 2.4.0)
################################################################################
model <- "
// Inferring the Mean and Standard Deviation of a Gaussian
data { 
  int<lower=1> n;
  vector<lower=0>[n] x;
}
parameters {
  real mu;
  real<lower=0,upper=10> sigma; 
} 
model {
  // Priors
  mu ~ normal(0, sqrt(1000));

  // Data Come From A Gaussian
  x ~ normal(mu, sigma);
}"

x <- c(1.1, 1.9, 2.3, 1.8)
n <- length(x)

data <- list(x=x, n=n) # to be passed on to Stan
myinits <- list(
  list(mu=0, sigma=1))

# parameters to be monitored: 
parameters <- c("mu", "sigma")

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
# Now the values for the monitored parameters are in the "samples" object, 
# ready for inspection.

mu <- extract(samples)$mu
sigma <- extract(samples)$sigma 
