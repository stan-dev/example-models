# clears workspace: 
rm(list=ls()) 

library(rstan)

model <- "
// Difference Between Two Rates
data { 
  int<lower=1> n1; 
  int<lower=1> n2; 
  int<lower=0> k1;
  int<lower=0> k2;
} 
parameters {
  real<lower=0,upper=1> theta1;
  real<lower=0,upper=1> theta2;
} 
transformed parameters {
  real<lower=-1,upper=1> delta;
  delta <- theta1 - theta2;
}
model {
  // Prior Distribution for Rate Theta
  theta1 ~ beta(1, 1);
  theta2 ~ beta(1, 1);
  // Observed Counts
  k1 ~ binomial(n1, theta1);
  k2 ~ binomial(n2, theta2);
}"

k1 <- 5
k2 <- 7
n1 <- 10
n2 <- 10

data <- list(k1=k1, k2=k2, n1=n1, n2=n2)  # to be passed on to Stan

myinits <- list(
  list(theta1=0.1, theta2=0.9))

# parameters to be monitored:  
parameters <- c("delta", "theta1", "theta2")

# The following command calls Stan with specific options.
# For a detailed description type "?rstan".
samples <- stan(model_code=model,   
                data=data, 
                init=myinits,  # If not specified, gives random inits
                pars=parameters,
                iter=20000, 
                chains=1, 
                thin=1,
                # warmup=100,  # Stands for burn-in; Default = iter/2
                # seed=123  # Setting seed; Default is random seed
                )
# Now the values for the monitored parameters are in the "samples" object, 
# ready for inspection.

print(samples)

# Collect posterior samples:
delta <- extract(samples)$delta

# Now let's plot a histogram for delta. 
# First, some options to make the plot look better:
par(cex.main = 1.5, mar = c(5, 6, 4, 5) + 0.1, mgp = c(3.5, 1, 0), cex.lab = 1.5,
    font.lab = 2, cex.axis = 1.3, bty = "n", las=1)
Nbreaks <- 80
y       <- hist(delta, Nbreaks, plot=F)
plot(c(y$breaks, max(y$breaks)), c(0,y$density,0), type="S", lwd=2, lty=1,
     xlim=c(-1,1), ylim=c(0,10), xlab="Difference in Rates", 
     ylab="Posterior Density") 

# mean of delta:
mean(delta)
# median of delta:
median(delta)
# mode of delta, estimated from the "density" smoother:
density(delta)$x[which(density(delta)$y==max(density(delta)$y))]
# 95% credible interval for delta:
quantile(delta, c(.025,.975))
