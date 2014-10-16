# clears workspace: 
rm(list=ls()) 

# Set working directory!

library(rstan)

model <- "
// Stop
data {
  int ns;
  int nq;
  int nc;
  int y[ns,nq];
  int m[83,nc];
  int p[nq,2];
  vector[nc] v;
  vector[nc] x;
}
transformed data {
  int<lower=1,upper=nc> s[nc];
  int<lower=0,upper=3> t[ns,nq,2];
  
  // Cue Search Order Follows Validities
  s <- sort_indices_asc(v);
  
  // TTB Decision
  for (i in 1:ns) {  
    for (q in 1:nq) {	
      vector[nc] tmp1;
      real tmp2;
      int tmp3;

      // Add Cue Contributions To Mimic TTB Decision
      for (j in 1:nc)
        tmp1[j] <- (m[p[q,1],j] - m[p[q,2],j]) * 2 ^ (s[j] - 1);
        
      // Find if Cue Favors First, Second, or Neither Stimulus
      tmp2 <- sum(tmp1);
      tmp3 <- -1 * int_step(-tmp2) + int_step(tmp2);
      t[i,q,1] <- tmp3 + 2;
    }
  }
  // WADD Decision
  for (i in 1:ns) {
    for (q in 1:nq) {
      vector[nc] tmp4;
      real tmp5;
      int tmp6;
      
      for (j in 1:nc)
        tmp4[j] <- (m[p[q,1],j] - m[p[q,2],j]) * x[j];
      
      // Find if Cue Favors First, Second, or Neither Stimulus
      tmp5 <- sum(tmp4);
      tmp6 <- -1 * int_step(-tmp5) + int_step(tmp5);
      t[i,q,2] <- tmp6 + 2;
    }
  }
}
parameters {
  real<lower=.5,upper=1> gamma;
  real<lower=0,upper=1> phi;
}
transformed parameters {
  vector<lower=0,upper=1>[3] dec;
  vector[2] lp_parts[ns];

  // Follow Decision With Probability Gamma, or Guess
  dec[1] <- 1 - gamma;
  dec[2] <- .5;
  dec[3] <- gamma;
  
  // TTB and WADD Subjects in Latent Mixture
  for (i in 1:ns) {
    vector[nq] lp_tmp1;
    vector[nq] lp_tmp2;
    for (q in 1:nq) {
      lp_tmp1[q] <- bernoulli_log(y[i,q], dec[t[i,q,1]]);
      lp_tmp2[q] <- bernoulli_log(y[i,q], dec[t[i,q,2]]);
    }
    lp_parts[i,1] <- log1m(phi) + sum(lp_tmp1);
    lp_parts[i,2] <- log(phi) + sum(lp_tmp2);
  }
} 
model {
  // Prior
  phi ~ beta(1, 1);
  
  for (i in 1:ns)
    increment_log_prob(log_sum_exp(lp_parts[i]));
}
generated quantities {
  int<lower=0,upper=1> ypred[ns,nq];
  int<lower=0,upper=1> z[ns];

  for (i in 1:ns) {
    vector[2] prob;    
    prob <- softmax(lp_parts[i]);
    z[i] <- bernoulli_rng(prob[2]);
  }
  for (q in 1:nq)
    for (i in 1:ns)
      ypred[i,q] <- bernoulli_rng(dec[t[i,q,(z[i] + 1)]]);
}"

load("StopSearchData.RData")  # Load all data for the model

# To be passed on to Stan
data <- list(nc=nc, nq=nq, ns=ns, v=v, p=p, m=m, y=y, x=x) 

myinits <- list(
  list(gamma=.75, phi=.5))

parameters <- c("gamma", "ypred", "z", "phi")  # Parameters to be monitored

# For a detailed description type "?stan".
samples <- stan(model_code=model,   
                data=data, 
                init=myinits,  # If not specified, gives random inits
                pars=parameters,
                iter=3500, 
                chains=1, 
                thin=1,
                warmup=500,  # Stands for burn-in; Default = iter/2
                # seed=123  # Setting seed; Default is random seed
)
# Now the values for the monitored parameters are in the "samples" object, 
# ready for inspection.

#### Figure 18.4 ####
means <- get_posterior_mean(samples)

get_mean <- function(var, i, j="", dim=2) {
  if (dim == 2)
    means[paste(var, "[", i, ",", j, "]", sep=""), ]
  else
    means[paste(var, "[", i, "]", sep=""), ]
}

windows(9, 5)
par(mar=c(3, 3, 1, 1) + .1, xaxs="i", mgp=c(1.3, 0.2, 0))
plot("", xlim=c(0, 31), ylim=c(1, 20), xlab="Question", ylab="Subject", xaxt="n", yaxt="n")
axis(1, c(1, 25, 30), tck=0)
axis(2, c(1, 5, 10, 15, 20), tck=0, las=2)
abline(v=24.5, lty=2)

for (i in 1:ns) {
  for (j in 1:nq) {
    points(j, i, pch=3, cex=get_mean("ypred", i, j))
    points(j, i, pch=4, cex=y[i, j])
  }
}

#### Figure 18.5 ####
windows(8, 3)
par(mar=c(3, 2, 1, 1) + .1, mgp=c(1.3, 0.2, 0), oma=c(0,3,0,0))
barplot(get_mean("z", 1:20, dim=1), col="black", names.arg=1:20,
        yaxt="n", xlab="Subject", ylab="")
box()
title(ylab="Group", outer=T)
axis(2, at=c(0, 1), labels=c("TTB", "WADD"), tck=0, las=2)
