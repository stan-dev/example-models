# clears workspace: 
rm(list=ls()) 

# Set working directory!

library(rstan)

model <- "
// Take The Best
data {
  int ns;
  int nq;
  int nc;
  int y[ns,nq];
  int m[83,nc];
  int p[nq,2];
  vector[nc] v;
}
transformed data {
  int<lower=1,upper=nc> s[nc];
  int<lower=0,upper=3> t[nq];
  
  // Cue Search Order Follows Validities
  s <- sort_indices_asc(v);
  
  // TTB Model For Each Question
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
    t[q] <- tmp3 + 2;
  }
}
parameters {
  real<lower=.5,upper=1> gamma;
}
transformed parameters {
  vector<lower=0,upper=1>[3] ttb;

  // Choose TTB Decision With Probability Gamma, or Guess
  ttb[1] <- 1 - gamma;
  ttb[2] <- .5;
  ttb[3] <- gamma;
} 
model {
  // Data
  for (q in 1:nq)
    for (i in 1:ns)
      y[i,q] ~ bernoulli(ttb[t[q]]);
}
generated quantities {
  int<lower=0,upper=1> ypred[ns,nq];

  for (q in 1:nq)
    for (i in 1:ns)
      ypred[i,q] <- bernoulli_rng(ttb[t[q]]);
}"

load("StopSearchData.RData")  # Load all data for the model

data <- list(nc=nc, nq=nq, ns=ns, v=v, p=p, m=m, y=y) # To be passed on to Stan

myinits <- list(
  list(gamma=.75))

parameters <- c("gamma", "ypred")  # Parameters to be monitored

# For a detailed description type "?stan".
samples <- stan(model_code=model,   
                data=data, 
                init=myinits,  # If not specified, gives random inits
                pars=parameters,
                iter=5500, 
                chains=1, 
                thin=1,
                warmup=500,  # Stands for burn-in; Default = iter/2
                # seed=123  # Setting seed; Default is random seed
)
# Now the values for the monitored parameters are in the "samples" object, 
# ready for inspection.

#### Figure 18.2 ####
means <- get_posterior_mean(samples)

get_mean <- function(i, j) {
  means[paste("ypred[", i, ",", j, "]", sep=""), ]
}

windows(9,5)
par(mar=c(3, 3, 1, 1) + .1, xaxs="i", mgp=c(1.3, 0, 0))
plot("", xlim=c(0, 31), ylim=c(1, 20), xlab="Question", ylab="Subject", xaxt="n", yaxt="n")
axis(1, c(1, 25, 30), tck=0)
axis(2, c(1, 5, 10, 15, 20), tck=0, las=2)
abline(v=24.5, lty=2)


for (i in 1:ns) {
  for (j in 1:nq) {
    points(j, i, pch=3, cex=get_mean(i, j))
    points(j, i, pch=4, cex=y[i, j])
  }
}