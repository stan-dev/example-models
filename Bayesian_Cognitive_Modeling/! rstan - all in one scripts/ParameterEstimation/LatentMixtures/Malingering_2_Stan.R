# clears workspace: 
rm(list=ls()) 

library(rstan)

model <- "
// Malingering, with Individual Differences
data { 
  int<lower=1> n;
  int<lower=1> p;
  int<lower=0,upper=n> k[p];
}
parameters {
  real<lower=0,upper=1> phi;
  real<lower=0,upper=1> mubon;
  real<lower=0> mudiff;
  real<lower=40,upper=800> lambdabon;
  real<lower=4,upper=100> lambdamal;
  matrix<lower=0,upper=1>[p,2] theta;
} 
transformed parameters {
  vector[2] lp_parts[p];
  vector<lower=0>[2] alpha;
  vector<lower=0>[2] beta;
  real<lower=0,upper=1> mumal;
    
  // Additivity on Logit Scale
  mumal <- inv_logit(logit(mubon) - mudiff);
  
  // Transformation to Group Mean and Precision
  alpha[1] <- mubon * lambdabon;
  beta[1] <- lambdabon * (1 - mubon);
  alpha[2] <- mumal * lambdamal;
  beta[2] <- lambdamal * (1 - mumal);
    
  // Data are Binomial with Rate Given by 
  // Each Person’s Group Assignment
  for (i in 1:p) {
    lp_parts[i,1] <- log1m(phi) + binomial_log(k[i], n, theta[i,1]);
    lp_parts[i,2] <- log(phi) + binomial_log(k[i], n, theta[i,2]);
  } 
}
model {
  // Priors
  mubon ~ beta(1, 1);  // can be removed
  mudiff ~ normal(0, 1 / sqrt(.5))T[0,];  // Constrained to be Positive
  // Relatively Uninformative Prior on Base Rate
  phi ~ beta(5, 5);
  
  for (i in 1:p)
    theta[i] ~ beta(alpha, beta);
    
  for (i in 1:p)
    increment_log_prob(log_sum_exp(lp_parts[i]));   
}
generated quantities {
  int<lower=0,upper=1> z[p];
    
  for (i in 1:p) {
    vector[2] prob;
    prob <- softmax(lp_parts[i]);
    // Each Person Belongs to One of Two Latent Groups
    z[i] <- bernoulli_rng(prob[2]);
  }
}"

k <- c(45, 45, 44, 45, 44, 45, 45, 45, 45, 45, 30,
       20, 6, 44, 44, 27, 25, 17, 14, 27, 35, 30)
p <- length(k) # number of people
n <- 45        # number of questions

data <- list(p=p, k=k, n=n) # To be passed on to Stan

myinits <- list(
    list(phi=.5, mubon=.5, mudiff=1, lambdabon=400, lambdamal=50,
         theta=matrix(rep(.5, p * 2), p, 2)))

# Parameters to be monitored
parameters <- c("mubon", "lambdabon", "mumal", "lambdamal", "mudiff", 
                "phi", "theta", "z")  

# The following command calls Stan with specific options.
# For a detailed description type "?rstan".
samples <- stan(model_code=model,   
                data=data, 
                init=myinits,  # If not specified, gives random inits
                pars=parameters,
                iter=3000, 
                chains=1, 
                thin=1,
                # warmup = 100,  # Stands for burn-in; Default = iter/2
                # seed = 123  # Setting seed; Default is random seed
)
# Now the values for the monitored parameters are in the "samples" object, 
# ready for inspection.

print(samples, digits=3)
