# clears workspace: 
rm(list=ls()) 

library(rstan)

#### Notes to Stan model #######################################################
## Stan has hypergeometric distribution implemented, so in one way the code
## is more intuitive. On the other hand, Stan can't sample discrete parameters,
## therefore we have to increment log probability manually (as we did in 
## Survey example).  
################################################################################
model <- "
// Planes
data { 
  int<lower=1> x;
  int<lower=1> n;
  int<lower=0,upper=n> k;
  int<lower=x> tmax;
}
transformed data {
  int<lower=x> tmin;
    
  tmin <- x + n - k;
}
transformed parameters {
  vector[tmax] lp_parts;

  for (t in 1:tmax) 
    if (t < tmin)
      lp_parts[t] <- log(1.0 / tmax) + negative_infinity();  // Zero probability
    else
      lp_parts[t] <- log(1.0 / tmax) + hypergeometric_log(k, n, x, t - x);
}
model {
  increment_log_prob(log_sum_exp(lp_parts));
}
generated quantities {
  int<lower=tmin,upper=tmax> t;
  simplex[tmax] tp;

  tp <- softmax(lp_parts);
  t <- categorical_rng(tp);
}"

x <- 10  # number of captures
k <- 4  # number of recaptures from n
n <- 5  # size of second sample
tmax <- 50  # maximum population size

data <- list(x=x, k=k, n=n, tmax=tmax) # to be passed on to Stan

parameters <- c("t")  # Parameters to be monitored

# The following command calls Stan with specific options.
# For a detailed description type "?rstan".
samples <- stan(model_code=model,   
                data=data, 
                init="random", 
                algorithm="Fixed_param",  # Since no parameters are sampled
                pars=parameters,
                iter=10000, 
                chains=1, 
                thin=1,
                # warmup = 100,  # Stands for burn-in; Default = iter/2
                # seed = 123  # Setting seed; Default is random seed
                )
# Now the values for the monitored parameters are in the "samples" object, 
# ready for inspection.

# Collect all samples in "t":
t <- extract(samples)$t

# Plot the posterior for theta:
windows(width=9,height=6) #Works only under Windows
par(cex.main = 1.5, mar = c(5, 6, 4, 5) + 0.1, mgp = c(3.5, 1, 0), cex.lab = 1.5,
    font.lab = 2, cex.axis = 1.3, bty = "n", las=1)
hist(t, xlim=c(x+n-k,tmax), lty=1, lwd=2, col="grey", prob=T, 
     breaks=(((x+n-k-1):tmax)+.5), axes=T, main=" ",
     xlab="Number of Planes", ylab="Posterior Mass")
