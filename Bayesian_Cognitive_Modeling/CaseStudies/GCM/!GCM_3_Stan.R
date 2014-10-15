# clears workspace: 
rm(list=ls()) 

library(rstan)

model <- "
// Generalized Context Model With Contaminants and Two Attention Groups
data { 
  int nstim;
  int nsubj;
  int n;
  int a[nstim];
  int y[nsubj,nstim];
  matrix[nstim,nstim] d1;
  matrix[nstim,nstim] d2;
}
transformed data {
  real b;
  b <- .5;
}
parameters {
  vector<lower=0>[nsubj] c;
  vector<lower=0,upper=1>[nsubj] w;
  vector<lower=0>[2] cpredg;
  vector<lower=0,upper=1>[2] wpredg;  
  real<lower=0,upper=1> phic;
  real<lower=0,upper=1> phig;
  real<lower=0,upper=1> muctmp;
  real<lower=0,upper=1> muwtmp;
  real<lower=0,upper=1> delta;
  real<lower=0,upper=1> sigmactmp;
  real<lower=0,upper=1> sigmawtmp;
} 
transformed parameters {
  matrix<lower=0,upper=1>[nsubj,nstim] r; 
  real<lower=0,upper=5> muc;
  vector<lower=0,upper=1>[2] muw;
  real<lower=.01,upper=3> sigmac;
  real<lower=.01,upper=1> sigmaw;
  vector[2] lp_parts_c[nsubj];
  vector[2] lp_parts_g[nsubj];

  // Mean Generalization
  muc <- 5 * muctmp;
  // Mean Attention
  muw[1] <- muwtmp;
  muw[2] <- fmin(1, delta + muw[1]);
  // Standard Deviation Generalization
  sigmac <- fmax(.01, 3 * sigmactmp);
  // Standard Deviation Attention
  sigmaw <- fmax(.01, sigmawtmp);

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
      r[k,i] <- sum(numerator) / sum(denominator);
    } 
  }  
  for (k in 1:nsubj) { 
    lp_parts_g[k,1] <- log1m(phig) + normal_log(w[k], muw[1], sigmaw)
                       - log(normal_cdf(1, muw[1], sigmaw) 
                             - normal_cdf(0, muw[1], sigmaw));
    lp_parts_g[k,2] <- log(phig) + normal_log(w[k], muw[2], sigmaw)
                       - log(normal_cdf(1, muw[2], sigmaw) 
                             - normal_cdf(0, muw[2], sigmaw));
  }
  for (k in 1:nsubj) {
    lp_parts_c[k,1] <- log1m(phic) + binomial_log(y[k], n, r[k]);
    lp_parts_c[k,2] <- log(phic) + binomial_log(y[k], n, .5);
  }
}
model {
  // Subject Parameters
  for (k in 1:nsubj)
    c[k] ~ normal(muc, sigmac)T[0,];

  // Predicted Group Parameters
  for (g in 1:2) {
    wpredg[g] ~ normal(muw[g], sigmaw)T[0,1];
    cpredg[g] ~ normal(muc, sigmac)T[0,];
  }

  // Decision Data
  for (k in 1:nsubj) 
    increment_log_prob(log_sum_exp(lp_parts_g[k]));

  for (k in 1:nsubj)   
    increment_log_prob(log_sum_exp(lp_parts_c[k]));
}
generated quantities {
  matrix[3,nstim] predyg;
  int<lower=0,upper=3> z[nsubj];

  for (i in 1:nstim) {
    matrix[2,nstim] numeratorpredg;
    matrix[2,nstim] denominatorpredg;
    vector[3] rpredg;
    for (j in 1:nstim) { 
      for (g in 1:2) {
        real spredg;
        spredg <- exp(-cpredg[g] * (wpredg[g] * d1[i,j] 
                                    + (1 - wpredg[g]) * d2[i,j]));
        numeratorpredg[g,j]   <- (a[j] == 1) * b * spredg;
        denominatorpredg[g,j] <- (a[j] == 1) * b * spredg
                                  + (a[j] == 2) * (1 - b) * spredg;
      }      
    }
    for (g in 1:2)
      rpredg[g] <- sum(numeratorpredg[g]) / sum(denominatorpredg[g]); 

    rpredg[3] <- 0.5;
    
    // Groups
    for (g in 1:3)
      predyg[g,i] <- binomial_rng(n, rpredg[g]);
  }

  for (k in 1:nsubj) {
    vector[2] prob_c;
    vector[2] prob_g;
    int zc;
    int zg;

    prob_c <- softmax(lp_parts_c[k]);
    prob_g <- softmax(lp_parts_g[k]);
    zc <- bernoulli_rng(prob_c[2]);
    zg <- bernoulli_rng(prob_g[2]);
    z[k] <- (zc == 0) * (zg + 1) + 3 * (zc == 1);
  }
}"

load("KruschkeData.Rdata")  # Load Kruschke's data

y <- t(y)  # Transpose matrix (for simpler Stan implementation)        

# To be passed on to Stan
data <- list(y=y, nstim=nstim, nsubj=nsubj, n=n, a=a, d1=d1, d2=d2) 

myinits <- list(
  list(c=abs(rnorm(nsubj, 1, .5)), w=runif(nsubj, 0, 1), cpredg=c(2.5, 2.5),
       wpredg=c(.5, .5), phic=.5, phig=.5, muctmp=.5, muwtmp=.25, delta=.5,
       sigmactmp=.5, sigmawtmp=.5),
  list(c=abs(rnorm(nsubj, 1, .5)), w=runif(nsubj, 0, 1), cpredg=c(2.5, 2.5),
       wpredg=c(.5, .5), phic=.5, phig=.5, muctmp=.5, muwtmp=.25, delta=.5,
       sigmactmp=.5, sigmawtmp=.5))

# Parameters to be monitored
parameters <- c("delta", "phic", "phig", "c", "w", "muc", "muw", "sigmac",
                "sigmaw", "wpredg", "cpredg", "predyg", "z")  

# The following command calls Stan with specific options.
# For a detailed description type "?stan".
samples <- stan(model_code=model,   
                data=data, 
                init=myinits,  # If not specified, gives random inits
                pars=parameters,
                iter=700, 
                chains=2, 
                thin=1,
                warmup=200,  # Stands for burn-in; Default = iter/2
                seed=541  # Setting seed; Default is random seed
)
# Now the values for the monitored parameters are in the "samples" object, 
# ready for inspection.
# For better estimates increase number of iterations or run optimized  
# model in file GCM_3_optimized_Stan.R

z <- extract(samples)$z
predyg <- extract(samples)$predyg

############# Figure 17.9 ################
z1 <- c()
z2 <- c()
z3 <- c()
for (i in 1:40) {
    z1[i] <- sum(z[, i] == 1) / length(z[, 1])
    z2[i] <- sum(z[, i] == 2) / length(z[, 2])
    z3[i] <- sum(z[, i] == 3) / length(z[, 3])
}

ord1 <- order(z1, decreasing=TRUE)
ord2 <- order(z2, decreasing=TRUE)
ord3 <- order(z3, decreasing=TRUE)
ord1 <- ord1[z1[ord1] > .5]
ord2 <- ord2[z2[ord2] > .5]
ord3 <- ord3[z3[ord3] > .5]
ord <- c(ord3, ord2, ord1)

windows(7, 4)
plot(z1[ord], ylim=c(0, 1), type="b", pch=0, ylab="Membership Probability",
     xlab="Subject")
lines(z2[ord], type="b", pch=1)
lines(z3[ord], type="b", pch=2)
legend("center", c("Contaminant", "Attend Position", "Attend Height"),
       pch=c(2:0))

############# Figure 17.10 ################
textMain <- c("Attend Height", "Attend Position", "Contaminant")
textXlab <- c("", "", "Stimulus")
textYlab <- c("", "", "Category")
nsamples <- length(predyg[, 1, 1])
squaresize <- 8
whichGroup <- apply(cbind(z1, z2, z3), 1, which.max)
xaxes <- c("n", "n", "s")

windows(14,7)
layout(matrix(1:3, 1, 3))
par(cex.main=1.5, mar=c(5, 5, 4, 1) + 0.1, mgp=c(3.5, 1, 0), cex.lab=1.5,
    font.lab=2, cex.axis=1.3)

for (p in 3:1) {  # loop over the plots
  plot(0, type="n", main=textMain[p], xlab=textXlab[p], ylab=textYlab[p], 
       xlim=c(1, 8), ylim=c(0, 8), xaxt=xaxes[p], yaxt="n")
  if (p == 3)
    Axis(side=2, at=c(0, 8), labels=c("B", "A"), las=2, cex=2)
  
  for (n in 1:40) {  # plotting observed data
    if (whichGroup[n] == p)
      lines(y[n, ], type="l", col="gray", lwd=1.5)
  }
  
  # plotting mean of observed data
  tmpMean <- apply(y[whichGroup == p, ], 2, mean)
  lines(tmpMean, lwd=4)
  
  for (i in 1:8) {  # plotting posterior predictive distributions
    for (j in 0:8) {
      tmp <- sum(predyg[, p, i] == j) / nsamples
      if (tmp > 0)
        points(i, j, pch=0, cex=squaresize * sqrt(tmp), lwd=1.2)  
    }
  }
}
