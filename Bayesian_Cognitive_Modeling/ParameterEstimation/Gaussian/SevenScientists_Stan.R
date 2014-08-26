# clears workspace: 
rm(list=ls()) 

library(rstan)

model <- "
// The Seven Scientists
data { 
  int<lower=1> n;
  vector[n] x;
}
parameters {
  real mu;
  vector<lower=0>[n] lambda;
} 
transformed parameters {
  vector[n] sigma;
  
  for (i in 1:n)
    sigma[i] <- inv_sqrt(lambda[i]);
}
model {
  // Priors
  mu ~ normal(0, sqrt(1000));
  lambda ~ gamma(.001, .001);
  
  // Data Come From Gaussians With Common Mean But Different Precisions
  x ~ normal(mu, sigma);
}"

x <- c(-27.020, 3.570, 8.191, 9.898, 9.603, 9.945, 10.056)
n <- length(x)

data <- list(x=x, n=n) # to be passed on to Stan
myinits <- list(
  list(mu=0, lambda=rep(1,n)))

# parameters to be monitored:  
parameters <- c("mu", "sigma")

# The following command calls Stan with specific options.
# For a detailed description type "?rstan".
samples <- stan(model_code=model,   
                data=data, 
                init=myinits,  # If not specified, gives random inits
                pars=parameters,
                iter=20000, 
                chains=1, 
                thin=1,
                # warmup = 100,  # Stands for burn-in; Default = iter/2
                # seed = 123  # Setting seed; Default is random seed
)
# Now the values for the monitored parameters are in the "samples" object, 
# ready for inspection.

