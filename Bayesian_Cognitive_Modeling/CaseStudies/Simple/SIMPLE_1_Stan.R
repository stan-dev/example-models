# clears workspace: 
rm(list=ls()) 

# Set working directory!

library(rstan)

#### Notes to Stan model #######################################################
## 1) For resp[j] you have two possible lines to run. I suggest using the one
##    with inv_logit for better performance, results are the same.
################################################################################

model <- "
// SIMPLE Model
data { 
  int dsets;
  int y[40,dsets];
  int n[dsets];
  int listlength[dsets];
  int m[40,dsets];
}
parameters {
  vector<lower=0,upper=100>[dsets] c;
  vector<lower=0,upper=100>[dsets] s;
  vector<lower=0,upper=1>[dsets] t;
} 
transformed parameters {
  matrix[40,dsets] theta;

  // Similarities, Discriminabilities, and Response Probabilities
  for (x in 1:dsets) {
    for (i in 1:listlength[x]) {
      vector[listlength[x]] sim;
      vector[listlength[x]] resp;
      
      // Similarities
      for (j in 1:listlength[x])
        sim[j] <- exp(-c[x] * fabs(log(m[i,x]) - log(m[j,x])));
             
      for (j in 1:listlength[x]) {
        real disc;
        
        // Discriminabilities
        disc <- sim[j] / sum(sim);
        // Response Probabilities - use just one of following 2 lines
        resp[j] <- inv_logit(s[x] * (disc - t[x]));
//        resp[j] <- 1 / (1 + exp(-s[x] * (disc - t[x])));
      }
      // Free Recall Overall Response Probability
      theta[i,x] <- fmin(1, sum(resp));
    }
  }
}
model {
  // Prior
  t ~ beta(1, 1);  // can be removed
  
  // Observed Data
  for (x in 1:dsets)
    for (i in 1:listlength[x])
      y[i,x] ~ binomial(n[x], theta[i,x]);
}
generated quantities {
  real predpc[40,dsets];
  
  // Predicted Data
  for (x in 1:dsets) {
    for (i in 1:listlength[x]) {
      real predy;
      predy <- binomial_rng(n[x], theta[i,x]);
      predpc[i,x] <- predy / n[x];
    }
  }   
}"

y          <- matrix(scan("k_M.txt", sep=","), ncol=40, nrow=6, byrow=T) 
n          <- c(1440, 1280, 1520, 1520, 1200, 1280)
listlength <- c(10, 15, 20, 20, 30, 40)
pc         <- matrix(scan("pc_M.txt", sep=","), ncol=40, nrow=6, byrow=T) 

dsets <- 6 

# Loop over conditions
m <- y * 0
for (dset in 1:dsets) {
    if (dset==1) {
        nwords <- 10
        lag    <- 2
        offset <- 15
    } 
    if (dset==2) {
        nwords <- 15
        lag    <- 2
        offset <- 20
    } 
    if (dset==3) {
        nwords <- 20
        lag    <- 2
        offset <- 25
    } 
    if (dset==4) {
        nwords <- 20
        lag    <- 1
        offset <- 10
    } 
    if (dset==5) {
        nwords <- 30
        lag    <- 1
        offset <- 15
    } 
    if (dset==6) {
        nwords <- 40
        lag    <- 1
        offset <- 20
    } 
    # Temporal Offset For Free Recall
    m[dset, 1:nwords] <- offset + seq(from=(nwords - 1) * lag, by=-lag, to=0)
}

y <- t(y)
m <- t(m)

# To be passed on to Stan
data <- list(y=y, n=n, listlength=listlength, m=m, dsets=dsets) 

myinits <- list(  # beware, this model is very sensitive to initial values
  list(c=seq(14, 21, length=dsets), s=seq(7.5, 14, length=dsets), 
       t=rep(.6, .45, length=dsets)))

parameters <- c("c", "t", "s", "predpc")  # Parameters to be monitored

# The following command calls Stan with specific options.
# For a detailed description type "?stan".
samples <- stan(model_code=model,   
                data=data, 
                init=myinits,  # If not specified, gives random inits
                pars=parameters,
                iter=400, 
                chains=1, 
                thin=1,
                warmup=100,  # Stands for burn-in; Default = iter/2
                # seed=123  # Setting seed; Default is random seed
)
# NB. For better convergence, you might want to take a lunch break and run
# multiple chains, for more iterations.


#Figure 15.2
layout(matrix(c(
                7,1,2,3,
                7,4,5,6,
                8,8,8,8
                ),3,4,byrow=T), c(1,2,2,2), c(2,2,.5))
layout.show(8)
hm <- 20
ll <- listlength

for (dset in 1:dsets) {
    plot(-1,-1,xlim=c(0,40),ylim=c(0,1),xlab="",ylab="",las=1)
    for (i in 1:ll[dset]) { 
        data <- extract(samples)$predpc[, i, dset]
        points(i+runif(hm,0,1)*.1,data[ceiling(runif(hm,0,1)*length(data))],col="grey")
    }
    points(1:ll[dset],pc[dset,1:ll[dset]],xlim=c(0,40),ylim=c(0,1))
    lines(1:ll[dset],pc[dset,1:ll[dset]])
    
    box("plot")
}
par(mar=c(rep(0,4)))
plot.new()
text(.45,.5,"Probability Correct",cex=2.5,srt=90)
plot.new()
text(.5,.5,"Serial Position",cex=2.5,mar=c(rep(0,4)))
