# clears workspace: 
rm(list=ls()) 

library(rstan)

model <- "
// Inferring a Common Rate
data { 
  int<lower=1> n1; 
  int<lower=1> n2; 
  int<lower=0> k1;
  int<lower=0> k2;
} 
parameters {
  real<lower=0,upper=1> theta;
} 
model {
  // Prior on Single Rate Theta
  theta ~ beta(1, 1);
  // Observed Counts
  k1 ~ binomial(n1, theta);
  k2 ~ binomial(n2, theta);
}"

k1 <- 5
k2 <- 7
n1 <- 10
n2 <- 10

data <- list(k1=k1, k2=k2, n1=n1, n2=n2)  # to be passed on to Stan

myinits <- list(
  list(theta=0.5))

# parameters to be monitored:  
parameters <- c("theta")

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
print(samples)

theta <- extract(samples)$theta
 
# Now let's plot a histogram for theta. 
# First, some options to make the plot look better:
par(cex.main = 1.5, mar = c(5, 6, 4, 5) + 0.1, mgp = c(3.5, 1, 0), cex.lab = 1.5,
    font.lab = 2, cex.axis = 1.3, bty = "n", las=1)
Nbreaks <- 80
y       <- hist(theta, Nbreaks, plot=F)
plot(c(y$breaks, max(y$breaks)), c(0,y$density,0), type="S", lwd=2, lty=1,
     xlim=c(0,1), ylim=c(0,10), xlab="Rate", ylab="Posterior Density") 
