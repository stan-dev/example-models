# clears workspace: 
rm(list=ls()) 

# Set working directory!

library(rstan)

model <- "
// Hierarchical SDT With Parameter Expansion
data { 
  int<lower=1> k;
  int<lower=0> h[k];
  int<lower=0> f[k];
  int<lower=0> s;
  int<lower=0> n;
}
parameters {
  real muc;
  real mud;
  real<lower=0> lambdac;
  real<lower=0> lambdad;
  real<lower=0,upper=1> xic;
  real<lower=0,upper=1> xid;
  vector[k] deltac;
  vector[k] deltad;
}
transformed parameters {
  vector[k] d;
  vector[k] c;
  vector<lower=0,upper=1>[k] thetah;
  vector<lower=0,upper=1>[k] thetaf;
  real<lower=0> sigmacnew;
  real<lower=0> sigmadnew;
  real<lower=0> sigmac;
  real<lower=0> sigmad;

  sigmacnew <- inv_sqrt(lambdac);
  sigmadnew <- inv_sqrt(lambdad);
  sigmac <- fabs(xic) * sigmacnew
  sigmad <- fabs(xid) * sigmadnew
  
  // Discriminability and Bias
  c <- muc + xic * deltac;
  d <- mud + xid * deltad;
  
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
  xic ~ beta(1, 1);  // can be removed
  xid ~ beta(1, 1);  // can be removed
  lambdac ~ gamma(.1, .1);
  lambdad ~ gamma(.1, .1);
  
  deltac ~ normal(0, sigmacnew);
  deltad ~ normal(0, sigmadnew);
  
  // Observed counts
  h ~ binomial(s, thetah);
  f ~ binomial(n, thetaf);  
}"

source("heit_rotello.RData") #loads the data

niter   <- 10000
nburnin <- 1000

for (dataset in 1:2) {  # analyze both conditions

  if (dataset == 1)
    data <- std_i  # the induction data
  if (dataset == 2)
    data <- std_d  # the deduction data
  
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
    list(deltac=rep(0, k), deltad=rep(0, k), xic=.5, xid=.5,
         muc=0, mud=0, lambdac=1, lambdad=1))  
  
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
plot(imud[keep], imuc[keep], xlab = "", ylab="", axes = F,xlim = c(-1,6), 
     ylim = c(-3,3))
points(dmud[keep],dmuc[keep], col = "grey")
box(lty = 1)

par(mar=c(2,1,1,4))
plot(d.imuc$y, d.imuc$x, xlim = rev(c(0,2.5)),type = 'l', axes = F, xlab = "", 
     ylab = "",ylim = c(-3,3))
lines(d.dmuc$y, d.dmuc$x, col = "grey")
axis(4)
mtext(expression(paste(mu, "c")), side = 4,line = 2.3, cex = 1.3)
box(lty = 1)

par(mar=c(6,2,0,0))
plot(density(imud),zero.line = F ,main = "", ylab = "", xlab = "", 
     cex.lab = 1.3, axes = F, xlim = c(-1,6),ylim = c(0,1))
lines(density(dmud), col = "grey")
axis(1, at = c(-1, 0, 1, 2, 3, 4, 5, 6))
mtext(expression(paste(mu, "d")), side = 1.2,line = 2, cex = 1.3)
box(lty = 1)
