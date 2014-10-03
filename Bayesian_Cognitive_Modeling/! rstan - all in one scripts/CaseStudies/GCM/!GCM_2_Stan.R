# clears workspace: 
rm(list=ls()) 

library(rstan)

model <- "
# Generalized Context Model With Individual Differences
data { 
  int nstim;
  int nsubj;
  int n;
  int a[nstim];
  int y[nstim,nsubj];
  matrix[nstim,nstim] d1;
  matrix[nstim,nstim] d2;
}
transformed data {
  real b;
  
  b <- .5;
}
parameters {
  vector<lower=0,upper=5>[nsubj] c;
  vector<lower=0,upper=1>[nsubj] w;
} 
transformed parameters {
  matrix<lower=0,upper=1>[nstim,nsubj] r; 
  
  for (i in 1:nstim) {
    vector[nstim] numerator;
    vector[nstim] denominator;
    for (k in 1:nsubj) {
      for (j in 1:nstim) {
        real s;
        // Similarities
        s <- exp(-c[k] * (w[k] * d1[i,j] + (1 - w[k]) * d2[i,j])); 
        
        // Base Decision Probabilities
        numerator[j] <- (a[j] == 1) * b * s;
        denominator[j] <- (a[j] == 1) * b * s + (a[j] == 2) * (1 - b) * s;
      }
      // Decision Probabilities
      r[i,k] <- sum(numerator) / sum(denominator);
    } 
  }  
}
model {
  // Decision Data
  for (i in 1:nstim)
    y[i] ~ binomial(n, r[i]);
}
generated quantities {
  matrix[nstim,nsubj] predy;
  
  for (i in 1:nstim)
    for (k in 1:nsubj)
      predy[i,k] <- binomial_rng(n, r[i,k]);
}"

load("KruschkeData.Rdata")  # Load Kruschke's data

# To be passed on to Stan
data <- list(y=y, nstim=nstim, nsubj=nsubj, n=n, a=a, d1=d1, d2=d2) 

myinits <- list(
  list(c=rep(2.5, nsubj), w=rep(.5, nsubj)))

parameters <- c("c", "w", "predy")  # Parameters to be monitored

# The following command calls Stan with specific options.
# For a detailed description type "?stan".
samples <- stan(model_code=model,   
                data=data, 
                init=myinits,  # If not specified, gives random inits
                pars=parameters,
                iter=600, 
                chains=1, 
                thin=1,
                warmup=100,  # Stands for burn-in; Default = iter/2
                # seed=123  # Setting seed; Default is random seed
)
# Now the values for the monitored parameters are in the "samples" object, 
# ready for inspection.

#### Figure 17.5 ####
windows(12, 7)
layout(matrix(1:40, 5, 8, byrow=TRUE))
par(mar=c(1, 1, 2, 1) + .1, oma=c(3, 3, 0, 0), mgp=c(3, .3, 0))

for (i in 1:40) {
  plot(y[, i], type="l", ylim=c(0, 8), xaxt="n", yaxt="n", main=i, ylab="",
       xlab="")
  if (i == 33) {
    axis(side=1, at=1:8, tick=FALSE, line=0)
    axis(side=2, at=c(0, 8), labels=c("B", "A"), tick=FALSE)
  }
}
mtext("Stimulus", side=1, line=.8, at=.035, adj=0, cex=1, outer=TRUE)  
mtext("Category", side=2, line=.8, at=.035, adj=0, cex=1, outer=TRUE)  

#### Figure 17.7 ####
c <- extract(samples)$c
w <- extract(samples)$w

cMean <- apply(c, 2, mean)
wMean <- apply(w, 2, mean)
keep=sample(1:length(c[, 1]), size=20)

par(cex.lab=1.2)
plot("", xlim=c(0, 4), ylim=c(0,1), xlab="Generalization", xaxs="i", yaxs="i",
     ylab="Attention Weight")

for (i in 1:nsubj) {
  for (j in 1:length(keep)) {
    segments(cMean[i], wMean[i], c[keep[j], i], w[keep[j], i], col="gray")
  }
}
points(cMean, wMean, pch=16)

for (i in c(3, 31, 33))
  text(cMean[i], wMean[i], pos=4, labels = i, cex=1.3)

