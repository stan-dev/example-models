# clears workspace: 
rm(list=ls()) 

library(rstan)

#### Notes to Stan model #######################################################
## 1) All notes from previous model Correlation_1 also apply to this model.
## 2) If you change sigmaerror to c(0.03, 10) as suggested in excercise 5.2.2
##    warning statements will be more frequent and posterior less smooth.
################################################################################
model <- "
// Pearson Correlation With Uncertainty in Measurement
data { 
  int<lower=0> n;
  vector[2] x[n];
  vector[2] sigmaerror;
}
parameters {
  vector[2] mu;
  vector<lower=0>[2] lambda;
  real<lower=-1,upper=1> r;
  vector[2] y[n];
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
  y ~ multi_normal(mu, T);
  for (i in 1:n)
    x[i] ~ normal(y[i], sigmaerror);
}"

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

n <- nrow(x) # number of people/units measured

# precision of measurement:
sigmaerror = c(.03, 1)
# sigmaerror = c(.03, 10)

data <- list(x=x, n=n, sigmaerror=sigmaerror) # to be passed on to Stan
myinits <- list(
  list(r=0, mu=c(0, 0), lambda=c(1, 1), y=matrix(c(rep(1, n), rep(100, n)), n, 2)))

# parameters to be monitored:  
parameters <- c("r", "mu", "sigma")

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

r    = extract(samples)$r

#Frequentist point-estimate of r:
freq.r <- cor(x[,1],x[,2])

#make the two panel plot:
windows(width=9,height=6) #this command works only under Windows!
layout(matrix(c(1,2),1,2))
layout.show(2)
#some plotting options to make things look better:
par(cex.main = 1.5, mar = c(5, 6, 4, 5) + 0.1, mgp = c(3.5, 1, 0), cex.lab = 1.5,
    font.lab = 2, cex.axis = 1.3, bty = "n", las=1)
# data panel:    
plot(x[,1],x[,2], type="p", pch=19, cex=1)
for (i in 1:n) {
  lines(c(x[i,1]-sigmaerror[1],x[i,1]+sigmaerror[1]), c(x[i,2],x[i,2]))
  lines(c(x[i,1],x[i,1]), c(x[i,2]-sigmaerror[2],x[i,2]+sigmaerror[2]))
}
# correlation panel:
plot(density(r, from=-1,to=1), main="", ylab="Posterior Density",
     xlab="Correlation", lwd=2)
lines(c(freq.r, freq.r), c(0,100), lwd=2, lty=2)
