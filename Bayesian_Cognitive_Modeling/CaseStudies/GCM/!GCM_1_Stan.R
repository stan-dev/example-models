# clears workspace: 
rm(list=ls()) 

library(rstan)

model <- "
// Generalized Context Model
data { 
  int nstim;
  int t;
  int a[nstim];
  int y[nstim];
  matrix[nstim,nstim] d1;
  matrix[nstim,nstim] d2;
}
transformed data {
  real b;
  b <- .5;
}
parameters {
  real<lower=0,upper=5> c;
  real<lower=0,upper=1> w;
} 
transformed parameters {
  vector<lower=0,upper=1>[nstim] r; 
  real tmp1[nstim,nstim,2];
  real tmp2[nstim,nstim,2];
  
  for (i in 1:nstim) {
    vector[nstim] numerator;
    vector[nstim] denominator;
    for (j in 1:nstim) {
      real s;
      // Similarities
      s <- exp(-c * (w * d1[i,j] + (1 - w) * d2[i,j])); 
      // Decision Probabilities
      tmp1[i,j,1] <- b * s;
      tmp1[i,j,2] <- 0;
      tmp2[i,j,1] <- 0;
      tmp2[i,j,2] <- (1 - b) * s;
      
      numerator[j] <- tmp1[i,j,a[j]];
      denominator[j] <- tmp1[i,j,a[j]] + tmp2[i,j,a[j]];
    }
    r[i] <- sum(numerator) / sum(denominator);
  }  
}
model {
  // Prior
  w ~ beta(1, 1);
  
  // Decision Data
  y ~ binomial(t, r);
}
generated quantities {
  vector[nstim] predy;
  
  for (i in 1:nstim)
    predy[i] <- binomial_rng(t, r[i]);
}"

load("KruschkeData.Rdata")  # Load Kruschke's data

x <- y
y <- apply(y, 1, sum)
t <- n * nsubj

data <- list(y=y, nstim=nstim, t=t, a=a, d1=d1, d2=d2) # To be passed on to Stan

myinits <- list(
  list(c=2.5, w=.5))

parameters <- c("c", "w", "predy")  # Parameters to be monitored

# The following command calls Stan with specific options.
# For a detailed description type "?stan".
samples <- stan(model_code=model,   
                data=data, 
                init=myinits,  # If not specified, gives random inits
                pars=parameters,
                iter=2000, 
                chains=1, 
                thin=1,
                warmup=1000,  # Stands for burn-in; Default = iter/2
                # seed=123  # Setting seed; Default is random seed
)
# Now the values for the monitored parameters are in the "samples" object, 
# ready for inspection.

c <- extract(samples)$c
w <- extract(samples)$w
predy <- extract(samples)$predy

#### Figure 17.3 ####
plot(c, w, xlim=c(0, 5), ylim=c(0,1), xlab="Generalization", pch=4, cex=.4,
     ylab="Attention Weight")

#### Figure 17.4 ####
breaks <- seq(0, t, by=2)

windows(10, 5)
par(mgp=c(2, 1, 0), mar=c(4, 4, 2, 2) + .1)
plot(NA, xlim=c(0.5, 8.5), ylim=c(0, 320), xlab="Stimulus", yaxt="n", xaxt="n",
     ylab="Category Decision")
axis(side=1, 1:8)
axis(side=2, c(0, t), c("B", "A"))

for (i in 1:nstim) {
  counts=hist(predy[, i], plot=FALSE, breaks=breaks)$counts
  breaks=hist(predy[, i], plot=FALSE, breaks=breaks)$breaks
  
  segments(i - counts * .003, breaks, i + counts * .003, breaks, col="gray",
           lwd=3.5, lend=2)
}
apply(x * 40, 2, lines, lty=3, col="gray")
lines(apply(x * 40, 1, mean), lwd=3)
