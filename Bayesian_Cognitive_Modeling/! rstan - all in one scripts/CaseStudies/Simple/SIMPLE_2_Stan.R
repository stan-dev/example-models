# clears workspace: 
rm(list=ls()) 

library(rstan)

model <- "
// Hierarchical SIMPLE Model
data { 
  int dsets;
  int gsets;
  int y[40,dsets];
  int n[gsets];
  int listlength[gsets];
  matrix[50,gsets] m;
  int w[gsets];
}
parameters {
  real<lower=0,upper=100> c;
  real<lower=0,upper=100> s;
  vector<lower=-1,upper=1>[2] a;
} 
transformed parameters {
  matrix[50,gsets] theta;
  vector<lower=0,upper=1>[gsets] t;

  // Similarities, Discriminabilities, and Response Probabilities
  for (x in 1:gsets) {
    t[x] <- fmax(0, fmin(1, a[1] * w[x] + a[2]));

    for (i in 1:listlength[x]) {
      vector[listlength[x]] sim;
      vector[listlength[x]] resp;
      // Similarities
      for (j in 1:listlength[x])
        sim[j] <- exp(-c * fabs(log(m[i,x]) - log(m[j,x])));
             
      for (j in 1:listlength[x]) {
        real disc;
        // Discriminabilities
        disc <- sim[j] / sum(sim);
        // Response Probabilities - use just one of following 2 lines
        resp[j] <- inv_logit(s * (disc - t[x]));
//        resp[j] <- 1 / (1 + exp(-s * (disc - t[x])));
      }
      // Free Recall Overall Response Probability
      theta[i,x] <- fmin(1, sum(resp));
    }
  }
}
model {
  // Priors
  a[1] ~ uniform(-1, 0);
  a[2] ~ uniform(0, 1);
  
  // Observed Data
  for (x in 1:dsets)
    for (i in 1:listlength[x])
      y[i,x] ~ binomial(n[x], theta[i,x]);
}
generated quantities {
  real predpc[50,gsets];
  
  // Predicted Data
  for (x in 1:gsets) {
    for (i in 1:listlength[x]) {
      real predy;
      predy <- binomial_rng(n[x], theta[i,x]);
      predpc[i,x] <- predy / n[x];
    }
  }   
}"

# read the data
y          <- matrix(scan("k_M.txt", sep=","), ncol=40, nrow=6, byrow=T) 
n          <- c(1440, 1280, 1520, 1520, 1200, 1280)
listlength <- c(10, 15, 20, 20, 30, 40)
pc         <- matrix(scan("pc_M.txt", sep=","), ncol=40, nrow=6, byrow=T) 
labs       <- scan("labs_M.txt", sep="", what="character")
dsets <- 6 
gsets <- 9

# Set Dataset to Use
w <- array()
l <- array()
m <- matrix(rep(0, 50 * 9), ncol=50, nrow=9, byrow=T) 
for (dset in 1:gsets) {
  if (dset == 1) {
    nwords <- 10
    lag    <- 2
    offset <- 15
  } 
  if (dset == 2) {
    nwords <- 15
    lag    <- 2
    offset <- 20
  } 
  if (dset == 3) {
    nwords <- 20
    lag    <- 2
    offset <- 25
  } 
  if (dset == 4) {
    nwords <- 20
    lag    <- 1
    offset <- 10
  } 
  if (dset == 5) {
    nwords <- 30
    lag    <- 1
    offset <- 15
  } 
  if (dset == 6) {
    nwords <- 40
    lag    <- 1
    offset <- 20
  }
  #Generalization: 
  if (dset == 7) {
    nwords <- 10
    lag    <- 1
    offset <- 5
  } 
  if (dset == 8) {
    nwords <- 25
    lag    <- 1
    offset <- 12.5
  } 
  if (dset == 9) {
    nwords <- 50
    lag    <- 1
    offset <- 25
  } 
  # Temporal Offset For Free Recall
  m[dset, 1:nwords] <- offset + seq(from=(nwords - 1) * lag, by=-lag, to=0)
  w[dset]          <- nwords 
  l[dset]          <- lag
  listlength[dset] <- nwords
}

n[(dsets + 1):gsets] <- 1200
labs = c(labs,c("10-1", "25-1", "50-1"))

y <- t(y)
m <- t(m)

# To be passed on to Stan
data <- list(gsets=gsets, y=y, n=n, listlength=listlength, m=m, w=w, dsets=dsets) 

myinits <- list(  # beware, this model is very sensitive to initial values
  list(c=21, a=c(-.002, .64), s=10))

parameters <- c("c", "s", "t", "a", "predpc")  # Parameters to be monitored

# The following command calls Stan with specific options.
# For a detailed description type "?rstan".
samples <- stan(model_code=model,   
                data=data, 
                init=myinits,  # If not specified, gives random inits
                pars=parameters,
                iter=200, 
                chains=1, 
                thin=1,
                warmup=50,  # Stands for burn-in; Default = iter/2
                seed=12  # Setting seed; Default is random seed
)
# NB. For better convergence, you might want to take a lunch break and run
# multiple chains, for more iterations, and do some thinning.

###Figure 15.5

layout(matrix(c(
                10,1,2,3,
                10,4,5,6,
                10,7,8,9,
                11,11,11,11
                ),4,4,byrow=T), c(1,2,2,2), c(2,2,2,1))
layout.show(11)
hm <- 20
ll <- listlength

for (dset in 1:gsets) {
  plot(-1,-1,xlim=c(0,50),ylim=c(0,1),xlab="",ylab="",las=1)
  text(47,.96,labs[dset])
  for (i in 1:ll[dset]) {                 
    data <- extract(samples)$predpc[, i, dset] 
    points(i+runif(hm,0,1)*.1,data[ceiling(runif(hm,0,1)*length(data))],col="grey")
  }
  if (dset <= 6) {
    points(1:ll[dset],pc[dset,1:ll[dset]],xlim=c(0,50),ylim=c(0,1))
    lines(1:ll[dset],pc[dset,1:ll[dset]])
  }
  box("plot")
}
par(mar=c(rep(0,4)))
plot.new()
text(.45,.5,"Probability Correct",cex=2.5,srt=90)
plot.new()
text(.5,.5,"Serial Position",cex=2.5,mar=c(rep(0,4)))


###Figure 15.6
layout(matrix(1:3,1,3))
layout.show(3)
#Threshold
epss <- .05
sx  <- seq(9,11,epss)
sxe <- seq(9+epss/2,11-epss/2,epss)
S <- extract(samples)$s
count <- hist(S,breaks=sx,plot=F)
count <- count$counts
count <- count/sum(count)*epss
keep <- which(count>1e-12)

plot(sxe[keep],count[keep],type="l",ylim=c(0,0.015),xlim=c(9,11),
     xlab="Threshold Noise (s)", ylab="Posterior Density", cex.lab=1.2,axes=F)
axis(1)
axis(2,labels=F,lwd.ticks=0)
box("plot")

#Distinctiveness
epsc <- .17
cx <- seq(18,24,epsc) #mids are cxe in R
cxe <- seq(18+epsc/2,24-epsc/2,epsc)
C <- extract(samples)$c
count <- hist(C,breaks=cx,plot=F)
count <- count$counts
count <- count/sum(count)*epsc
keep <- which(count>1e-12)

plot(cxe[keep],count[keep],type="l",ylim=c(0,0.06),xlim=c(18,23),
     xlab="Distinctiveness (c)", ylab="Posterior Density", cex.lab=1.2, axes=F)
axis(1)
axis(2,labels=F,lwd.ticks=0)
box("plot")

###Threshold Parameter as a Function Of List Length
howmany <- 50
nsamples <- length(C)
keep <- ceiling(runif(howmany,min=0,max=1)*nsamples)
wdom <- seq(1,50,1)
plot(-1, -1, xlim=c(1,50), ylim=c(0,1),xlab="Item List Length (W)", 
     ylab="Treshold (t)")
for (i in 1:howmany){
  predt <- extract(samples)$a[, 1][[keep[i]]] * wdom +
           extract(samples)$a[, 2][[keep[i]]]
  points(wdom,predt,col="grey")
}
