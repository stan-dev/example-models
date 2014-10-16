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
  real<lower=0,upper=5> muc;  // Mean Generalization
  real<lower=0,upper=1> muwtmp;
  real<lower=0,upper=1> delta;
  real<lower=0,upper=3> sigmac;  // Standard Deviation Generalization
  real<lower=0,upper=1> sigmaw;  // Standard Deviation Attention
} 
transformed parameters {
  vector[nstim] log_r[nsubj]; 
  vector<lower=0,upper=1>[2] muw;
  vector[2] lp_parts_c[nsubj];
  vector[2] lp_parts_g[nsubj];

  // Mean Attention
  muw[1] <- muwtmp;
  muw[2] <- fmin(1, delta + muw[1]);

  for (i in 1:nstim) {
    vector[nstim] numerator;
    vector[nstim] denominator;
    for (k in 1:nsubj) {
      for (j in 1:nstim) {
        real log_s;
        // Similarities
        log_s <- -c[k] * (w[k] * d1[i,j] + (1 - w[k]) * d2[i,j]);
        if (a[j] == 1) {
          numerator[j] <- log(b) + log_s;
          denominator[j] <- numerator[j];
        } else {
          numerator[j] <- negative_infinity();
          denominator[j] <- log1m(b) + log_s;
        }
      }
      log_r[k,i] <- log_sum_exp(numerator) - log_sum_exp(denominator);
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
    vector[nstim] tmp;
    
     // Binomial distribution written as Log Probability Mass Function
    for (i in 1:nstim)
      tmp[i] <- lgamma(n + 1) - lgamma(y[k,i] + 1) - lgamma(n - y[k,i] + 1) 
                + y[k,i] * log_r[k,i] + (n - y[k,i]) * log1m_exp(log_r[k,i]);
  
    lp_parts_c[k,1] <- log1m(phic) + sum(tmp);
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
    vector[3] rpredg_log;
    for (g in 1:2) {
      vector[nstim] numeratorpredg;
      vector[nstim] denominatorpredg;
      for (j in 1:nstim) { 
        real spredg_log;
        spredg_log <- -cpredg[g] * (wpredg[g] * d1[i,j]
                                    + (1 - wpredg[g]) * d2[i,j]);
        if (a[j] == 1) {
          numeratorpredg[j] <- log(b) + spredg_log;
          denominatorpredg[j] <- numeratorpredg[j];
        } else {
          numeratorpredg[j] <- negative_infinity();
          denominatorpredg[j] <- log1m(b) + spredg_log;
        }
      }      
      rpredg_log[g] <- log_sum_exp(numeratorpredg) - log_sum_exp(denominatorpredg); 
    }
    rpredg_log[3] <- log(.5);
    
    // Groups
    for (g in 1:3) {
      real rpredg;
      rpredg <- exp(rpredg_log[g]);
      predyg[g,i] <- binomial_rng(n, rpredg);
    }
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

set.seed(1234)
myinits <- list(
  list(c=abs(rnorm(nsubj, 1, .5)), w=runif(nsubj, 0, 1), cpredg=c(2.5, 2.5),
       wpredg=c(.5, .5), phic=.5, phig=.5, muc=2.5, muwtmp=.25, delta=.5),
  list(c=abs(rnorm(nsubj, 1, .5)), w=runif(nsubj, 0, 1), cpredg=c(2.5, 2.5),
       wpredg=c(.5, .5), phic=.5, phig=.5, muc=2.5, muwtmp=.25, delta=.5))

# Parameters to be monitored
parameters <- c("delta", "phic", "phig", "c", "w", "muc", "muw", "sigmac",
                "sigmaw", "wpredg", "cpredg", "predyg", "z")  

# The following command calls Stan with specific options.
# For a detailed description type "?stan".
samples <- stan(model_code=model,   
                data=data, 
                init="random",  # If not specified, gives random inits
                pars=parameters,
                iter=600, 
                chains=2, 
                thin=1,
                warmup=100,  # Stands for burn-in; Default = iter/2
                # seed=541  # Setting seed; Default is random seed
)
# Now the values for the monitored parameters are in the "samples" object, 
# ready for inspection.
# For better estimates increase number of iterations. 

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
