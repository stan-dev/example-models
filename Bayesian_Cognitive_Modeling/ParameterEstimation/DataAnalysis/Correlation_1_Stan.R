# clears workspace: 
rm(list=ls()) 

library(rstan)

#### Notes to Stan model #######################################################
## 1) Multivariate normal distribution in Stan uses covariance matrix instead of 
##    precision matrix.
## 2) Multivariate normal distribution can be (and is) also vectorized.
## 3) Warnings may occur during sampling, ignore them.
################################################################################
model <- "
// Pearson Correlation
data { 
  int<lower=0> n;
  vector[2] x[n];
}
parameters {
  vector[2] mu;
  vector<lower=0>[2] lambda;
  real<lower=-1,upper=1> r;
} 
transformed parameters {
  vector<lower=0>[2] sigma;
  cov_matrix[2] T;

  // Reparameterization
  sigma[1] <- inv_sqrt(lambda[1]);
  sigma[2] <- inv_sqrt(lambda[2]);
  T[1,1] <- square(sigma[1]);
  T[1,2] <- r * sigma[1] * sigma[2];
  T[2,1] <- r * sigma[1] * sigma[2];
  T[2,2] <- square(sigma[2]);
}
model {
  // Priors
  mu ~ normal(0, inv_sqrt(.001));
  lambda ~ gamma(.001, .001);
  
  // Data
  x ~ multi_normal(mu, T);
}"

# Choose a dataset:
dataset <- 1

# The datasets:
if (dataset == 1) { 
  x <- matrix(c( .8, 102, 
                1.0,  98, 
                 .5, 100,
                 .9, 105, 
                 .7, 103, 
                 .4, 110,
                1.2,  99, 
                1.4,  87,
                 .6, 113,
                1.1,  89,
                1.3,  93), nrow=11, ncol=2, byrow=T) 
}
if (dataset == 2) {
  x <- matrix(c( .8, 102, 
                1.0,  98,
                 .5, 100,
                 .9, 105,
                 .7, 103, 
                 .4, 110,
                1.2,  99,
                1.4,  87,
                 .6, 113, 
                1.1,  89,
                1.3,  93,
                 .8, 102, 
                1.0,  98, 
                 .5, 100, 
                 .9, 105, 
                 .7, 103, 
                 .4, 110, 
                1.2,  99, 
                1.4,  87, 
                 .6, 113, 
                1.1,  89, 
                1.3,  93), nrow=22,ncol=2,byrow=T) 
}

n <- nrow(x) # number of people/units measured

data <- list(x=x, n=n) # to be passed on to Stan
myinits <- list(
  list(r=0, mu=c(0, 0), lambda=c(1, 1)))

# parameters to be monitored: 
parameters <- c("r", "mu", "sigma")

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

r <- extract(samples)$r

#Frequentist point-estimate of r:
freq.r <- cor(x[,1],x[,2])

#make the two panel plot:
windows(width=9,height=6) #this command works only under Windows!
layout(matrix(c(1,2),1,2))
layout.show(2)
#some plotting options to make things look better:
par(cex.main=1.5, mar=c(5, 6, 4, 5) + 0.1, mgp=c(3.5, 1, 0), cex.lab=1.5,
    font.lab=2, cex.axis=1.3, bty = "n", las=1)
# data panel:    
plot(x[,1],x[,2], type="p", pch=19, cex=1)
# correlation panel:
plot(density(r, from=-1,to=1), main="", ylab="Posterior Density", 
     xlab="Correlation", lwd=2)
lines(c(freq.r, freq.r), c(0,100), lwd=2, lty=2)
