# clears workspace: 
rm(list=ls()) 

library(rstan)

model <- "
// Signal Detection Theory
data { 
  int<lower=1> k;
  int<lower=0> h[k];
  int<lower=0> f[k];
  int<lower=0> s[k];
  int<lower=0> n[k];
}
parameters {
  vector[k] d;
  vector[k] c;
} 
transformed parameters {
  real<lower=0,upper=1> thetah[k];
  real<lower=0,upper=1> thetaf[k];
  
  // Reparameterization Using Equal-Variance Gaussian SDT
  for(i in 1:k) {
    thetah[i] <- Phi(d[i] / 2 - c[i]);
    thetaf[i] <- Phi(-d[i] / 2 - c[i]);
  }
}
model {
  // These Priors over Discriminability and Bias Correspond 
  // to Uniform Priors over the Hit and False Alarm Rates
  d ~ normal(0, inv_sqrt(.5));
  c ~ normal(0, inv_sqrt(2));
  
  // Observed counts
  h ~ binomial(s, thetah);
  f ~ binomial(n, thetaf);
}"

dataset <- 1

if (dataset == 1) {  # Demo
  k <- 3 # number of cases
  data <- matrix(c(70, 50, 30, 50,
                    7,  5,  3,  5, 
                   10,  0,  0, 10), nrow=k, ncol=4, byrow=TRUE)
}

if (dataset == 2) {  # Lehrner et al. (1995) data 
  k <- 3 # number of cases
  data <- matrix(c(148, 29, 32, 151,
                   150, 40, 30, 140,
                   150, 51, 40, 139), nrow=k, ncol=4, byrow=TRUE)
}

h <- data[, 1]
f <- data[, 2]
MI <- data[, 3]
CR <- data[, 4]
s <- h + MI
n <- f + CR

data <- list(h=h, f=f, s=s, n=n, k=k) # To be passed on to Stan

myinits <- list(
  list(d=rep(0, k), c=rep(0, k)))

# Parameters to be monitored
parameters <- c("c", "d", "thetaf", "thetah")

# The following command calls Stan with specific options.
# For a detailed description type "?rstan".
samples <- stan(model_code=model,   
                data=data, 
                init=myinits,  # If not specified, gives random inits
                pars=parameters,
                iter=20000, 
                chains=1, 
                thin=1,
                # warmup=100,  # Stands for burn-in; Default = iter/2
                # seed=123  # Setting seed; Default is random seed
)
# Now the values for the monitored parameters are in the "samples" object, 
# ready for inspection.

d1 <- extract(samples)$d[,1]
d2 <- extract(samples)$d[,2]
d3 <- extract(samples)$d[,3]

c1 <- extract(samples)$c[,1]
c2 <- extract(samples)$c[,2]
c3 <- extract(samples)$c[,3]

h1 <- extract(samples)$thetah[,1]
h2 <- extract(samples)$thetah[,2]
h3 <- extract(samples)$thetah[,3]

f1 <- extract(samples)$thetaf[,1]
f2 <- extract(samples)$thetaf[,2]
f3 <- extract(samples)$thetaf[,3]

#make the four panel plot:
layout(matrix(c(1,2,3,4), 2, 2, byrow = TRUE))
#layout.show(4)
#some plotting options to make things look better:
par(cex.main = 1.5, mar = c(5, 6, 4, 5) + 0.1, mgp = c(3.5, 1, 0), cex.lab = 1.5,
    font.lab = 2, cex.axis = 1.3, bty = "n", las=1)
# Discriminability panel:    
plot(density(d1), lwd=2, col="red", main="", ylab="", xlab="", 
     xlim=c(-2,6), axes=F)
axis(1)
axis(2, labels=F, at=c(0,24))

lines(density(d2), lwd=2, col="green", lty=2)
lines(density(d3), lwd=2, col="blue", lty=2)

mtext("Probability Density", side=2, line = 2, cex=1.5, las=0)
mtext("Discriminability", side=1, line = 2.5, cex=1.5)

# Bias panel:    
plot(density(c1), lwd=2, col="red", main="", ylab="", xlab="", 
     xlim=c(-2,2), axes=F)
axis(1)
axis(2, labels=F, at=c(0,24))

lines(density(c2), lwd=2, col="green", lty=2)
lines(density(c3), lwd=2, col="blue", lty=2)

mtext("Probability Density", side=2, line = 2, cex=1.5, las=0)
mtext("Bias", side=1, line = 2.5, cex=1.5)

# Hit Rate panel:    
plot(density(h1), lwd=2, col="red", main="", ylab="", xlab="", 
     xlim=c(0,1), axes=F)
axis(1)
axis(2, labels=F, at=c(0,24))

lines(density(h2), lwd=2, col="green", lty=2)
lines(density(h3), lwd=2, col="blue", lty=2)

if (dataset == 1)
{
  lines(c(0, 0.1),c(7,7), lwd=2, lty=1, col="red")
  lines(c(0, 0.1),c(6,6), lwd=2, lty=2, col="green")
  lines(c(0, 0.1),c(5,5), lwd=2, lty=3, col="blue")
  
  text(0.15, 7, labels="first", offset=0, cex = 1.3, pos=4)
  text(0.15, 6, labels="second", offset=0, cex = 1.3, pos=4)
  text(0.15, 5, labels="third", offset=0, cex = 1.3,pos=4)
}

if (dataset == 2)
{
  lines(c(0, 0.1),c(7,7), lwd=2, lty=1, col="red")
  lines(c(0, 0.1),c(6,6), lwd=2, lty=2, col="green")
  lines(c(0, 0.1),c(5,5), lwd=2, lty=3, col="blue")
  
  text(0.15, 7, labels="Control", offset=0, cex = 1.3, pos=4)
  text(0.15, 6, labels="Group I", offset=0, cex = 1.3, pos=4)
  text(0.15, 5, labels="Group II", offset=0, cex = 1.3,pos=4)
}

mtext("Probability Density", side=2, line = 2, cex=1.5, las=0)
mtext("Hit Rate", side=1, line = 2.5, cex=1.5)

# False-Alarm Rate panel:    
plot(density(f1), lwd=2, col="red", main="", ylab="", xlab="", 
     xlim=c(0,1), axes=F)
axis(1)
axis(2, labels=F, at=c(0,24))

lines(density(f2), lwd=2, col="green", lty=2)
lines(density(f3), lwd=2, col="blue", lty=2)

mtext("Probability Density", side=2, line = 2, cex=1.5, las=0)
mtext("False-Alarm Rate", side=1, line = 2.5, cex=1.5)


