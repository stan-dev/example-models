# clears workspace: 
rm(list=ls()) 

library(rstan)

model <- "
// Retention With No Individual Differences
data { 
  int ns;
  int nt;
  int k[ns - 1,nt - 1];
  int t[nt];
  int n;
}
parameters {
  real<lower=0,upper=1> alpha;
  real<lower=0,upper=1> beta;
} 
transformed parameters {
  matrix<lower=0,upper=1>[ns,nt] theta;
  
  // Retention Rate At Each Lag For Each Subject Decays Exponentially
  for (i in 1:ns) {
    for (j in 1:nt) {
      theta[i,j] <- fmin(1, exp(-alpha * t[j]) + beta);
    }
  }
}
model {
  // Priors
  alpha ~ beta(1, 1);  // can be removed
  beta ~ beta(1, 1);  // can be removed
  
  // Observed Data
  for (i in 1:(ns - 1))
    for (j in 1:(nt - 1))
      k[i,j] ~ binomial(n, theta[i,j]);
}
generated quantities {
  int<lower=0,upper=n> predk[ns,nt];
  
  // Predicted Data
  for (i in 1:ns)
    for (j in 1:nt)
      predk[i,j] <- binomial_rng(n, theta[i,j]);
}"

t     <- c(1, 2, 4, 7, 12, 21, 35, 59, 99, 200)
nt    <- length(t)
slist <- 1:4
ns    <- length(slist)

k1 <- matrix(c(18, 18, 16, 13, 9, 6, 4, 4, 4, NA,
               17, 13,  9,  6, 4, 4, 4, 4, 4, NA,
               14, 10,  6,  4, 4, 4, 4, 4, 4, NA,
               NA, NA, NA, NA,NA,NA,NA,NA,NA, NA), nrow=ns, ncol=nt, byrow=T)

k <- k1[1:(ns - 1), 1:(nt - 1)]   # Excluding NAs (for Stan solution)

n <- 18

data <- list(k=k, n=n, t=t, ns=ns, nt=nt) # To be passed on to Stan

myinits <- list(
  list(alpha=.5, beta=.1))

parameters <- c("alpha", "beta", "predk")  # Parameters to be monitored

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

##############Fancy plots################################################
##Figure 10.2

alpha <- extract(samples)$alpha
beta <- extract(samples)$beta
d.beta <- density(beta)

layout(matrix(c(1,2,3,0),2,2,byrow=T), width=c(2/3, 1/3), heights=c(2/3,1/3))
#layout.show()

par(mar=c(2,2,1,0))
plot(alpha,beta, xlab="", ylab="", xlim=c(0,1),ylim=c(0,1),axes=F)
box(lty=1)

par(mar=c(2,2,1,4))
plot(d.beta$y, d.beta$x, ylim=range(c(0,1)), xlim=rev(range(d.beta$y)), 
   type='l', axes=F, xlab="", ylab="")
axis(4, at=c(0,1))
mtext(expression(beta), side=4,line=1, cex=1.3)
box(lty=1)

par(mar=c(6,2,0,0))
plot(density(alpha),zero.line=F, main="", ylab="", xlab="", cex.lab=1.3,
   xlim=c(0,1), axes=F)
axis(1,at=c(0,1))
mtext(expression(alpha), side=1.2,line=1, cex=1.3)
box(lty=1)

##Figure 10.3

layout(matrix(c(1:4),2,2,byrow=T))
#layout.show()
sc <- 3.5
jj <- numeric()
xx <- numeric()

for (i in 1:ns) 
{
  plot(-1,100,xlim=c(0,10),ylim=c(0,18), main=(paste("Subject", i)),
     xlab=("Time Lags"), ylab=("Retention Count"),cex.lab=1.3, axes=F)
  axis(1, at=c(1,2,3,4,5,6,7,8,9,10), 
     lab=c("1","2","4","7","12","21","35","59","99","200"),cex.axis=0.7)
  axis(2, at=c(0,18),lab=c("0","18"),cex.axis=0.7)
  box(lty=1)
  for (j in 1:nt) 
  {
    count <- hist(extract(samples)$predk[,i,j],c(0:n),plot=F)
    count <- count$counts
    count <- count/sum(count)
    for (x in 1:n)
    {
      if (count[x]>0)
      {
        points(j,x,pch=22, col="black",cex=sc*sqrt(count[x]))
        if (!is.na(k1[i,j]) && k1[i,j]==x)
        {
          points(j,x,pch=22,bg="black",cex=sc*sqrt(count[x]))
          jj <- c(jj,j)
          xx <- c(xx,x)
        }
      }
    }
  }
  coords <- list(x=jj, y=xx)
  lines(coords,lwd=2)
  jj <- numeric()
  xx <- numeric()
}

