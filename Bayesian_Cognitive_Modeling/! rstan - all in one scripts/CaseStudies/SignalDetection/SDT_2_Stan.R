# clears workspace: 
rm(list=ls()) 

# Set working directory!

library(rstan)

model <- "
// Hierarchical Signal Detection Theory
data { 
  int<lower=1> k;
  int<lower=0> h[k];
  int<lower=0> f[k];
  int<lower=0> s;
  int<lower=0> n;
}
parameters {
  vector[k] d;
  vector[k] c;
  real muc;
  real mud;
  real<lower=0> lambdac;
  real<lower=0> lambdad;
} 

transformed parameters {
  real<lower=0,upper=1> thetah[k];
  real<lower=0,upper=1> thetaf[k];
  real<lower=0> sigmac;
  real<lower=0> sigmad;
  
  sigmac <- inv_sqrt(lambdac);
  sigmad <- inv_sqrt(lambdad);
  
  // Reparameterization Using Equal-Variance Gaussian SDT
  for(i in 1:k) {
    thetah[i] <- Phi(d[i] / 2 - c[i]);
    thetaf[i] <- Phi(-d[i] / 2 - c[i]);
  }
}
model {
  // Priors 
  muc ~ normal(0, inv_sqrt(.001));
  mud ~ normal(0, inv_sqrt(.001));
  lambdac ~ gamma(.001, .001);
  lambdad ~ gamma(.001, .001);
  
  // Discriminability and Bias
  c ~ normal(muc, sigmac);
  d ~ normal(mud, sigmad);
  // Observed counts
  h ~ binomial(s, thetah);
  f ~ binomial(n, thetaf);
}"

source("heit_rotello.RData") #loads the data

niter   <- 10000
nburnin <- 1000

for (dataset in 1:2) {  #analyze both conditions

  if (dataset == 1)
    data <- std_i # the induction data
  if (dataset == 2)
    data <- std_d # the deduction data
  
  h <- data[, 1]
  f <- data[, 2]
  MI <- data[, 3]
  CR <- data[, 4]
  s <- h + MI
  n <- f + CR
  s <- s[1]; n <- n[1] #Each subject gets same number of signal and noise trials 
  k <- nrow(data) 

  data <- list(h=h, f=f, s=s, n=n, k=k) # To be passed on to Stan

  myinits <- list(
    list(d=rep(0, k), c=rep(0, k), mud=0, muc=0, lambdad=1, lambdac=1)) 

  # Parameters to be monitored
  parameters <- c("mud", "muc", "sigmad", "sigmac")
  
  if (dataset == 1) {
    # The following command calls Stan with specific options.
    # For a detailed description type "?rstan".
    isamples <- stan(model_code=model,   
                     data=data, 
                     init=myinits,  # If not specified, gives random inits
                     pars=parameters,
                     iter=niter, 
                     chains=1, 
                     thin=1,
                     warmup=nburnin,  # Stands for burn-in; Default = iter/2
                     # seed=123  # Setting seed; Default is random seed
    )
  }
  if (dataset == 2) {
    # The following command calls Stan with specific options.
    # For a detailed description type "?rstan".
    dsamples <- stan(fit=isamples,   
                     data=data, 
                     init=myinits,  # If not specified, gives random inits
                     pars=parameters,
                     iter=niter, 
                     chains=1, 
                     thin=1,
                     warmup=nburnin,  # Stands for burn-in; Default = iter/2
                     # seed=123  # Setting seed; Default is random seed
    )
  }
}
# Now the values for the monitored parameters are in the "isamples" and 
# "dsamples "objects, ready for inspection.

#####Figure 11.5 & 11.6

keepi <- 1000
keep <- sample(niter, keepi)

imud <- extract(isamples)$mud
imuc <- extract(isamples)$muc
d.imuc <- density(imuc)

dmud <- extract(dsamples)$mud
dmuc <- extract(dsamples)$muc
d.dmuc <- density(dmuc)

layout(matrix(c(1,2,3,0),2,2,byrow=T), width=c(2/3, 1/3), heights=c(2/3,1/3))
#layout.show()

par(mar=c(2,2,1,0))
plot(imud[keep],imuc[keep], xlab="", ylab="", axes=F,xlim=c(-1,6), ylim=c(-3,3))
points(dmud[keep],dmuc[keep], col="grey")
box(lty=1)

par(mar=c(2,1,1,4))
plot(d.imuc$y, d.imuc$x, xlim=rev(c(0,2.5)),type='l', axes=F, xlab="", ylab="",
     ylim=c(-3,3))
lines(d.dmuc$y, d.dmuc$x, col="grey")
axis(4)
mtext(expression(paste(mu, "c")), side=4,line=2.3, cex=1.3)
box(lty=1)

par(mar=c(6,2,0,0))
plot(density(imud),zero.line=F ,main="", ylab="", xlab="", cex.lab=1.3, axes=F,
     xlim=c(-1,6),ylim=c(0,1))
lines(density(dmud), col="grey")
axis(1, at=c(-1, 0, 1, 2, 3, 4, 5, 6))
mtext(expression(paste(mu, "d")), side=1.2,line=2, cex=1.3)
box(lty=1)
