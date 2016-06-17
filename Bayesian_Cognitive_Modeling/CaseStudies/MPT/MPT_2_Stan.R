# clears workspace: 
rm(list=ls()) 

library(rstan)

model <- "
// Multinomial Processing Tree with Latent Traits
data { 
	int<lower=1> nsubjs; 
	int<lower=1> nparams; 
  int<lower=0,upper=20> k[nsubjs,4];
	cov_matrix[nparams] I;
}
transformed data {
	int df;
  vector[3] mudeltahat;
	
	df <- nparams + 1;
  mudeltahat[1] <- 0;
  mudeltahat[2] <- 0;
  mudeltahat[3] <- 0;
}
parameters {
	vector[nparams] deltahat[nsubjs]; 
	cov_matrix[nparams] Sigma;
	real muchat;
	real murhat;
	real muuhat;
	real<lower=0,upper=100> xichat;
	real<lower=0,upper=100> xirhat;
	real<lower=0,upper=100> xiuhat;
} 
transformed parameters {
  simplex[4] theta[nsubjs];
  vector<lower=0,upper=1>[nsubjs] c;
  vector<lower=0,upper=1>[nsubjs] r;
  vector<lower=0,upper=1>[nsubjs] u;
  matrix[nparams,nparams] rho;
		
  vector[nsubjs] deltachat;
  vector[nsubjs] deltarhat;
  vector[nsubjs] deltauhat;

	for (i in 1:nsubjs) {
    
		deltachat[i] <- deltahat[i,1];
		deltarhat[i] <- deltahat[i,2];
		deltauhat[i] <- deltahat[i,3];
		
		// Probitize Parameters c, r, and u 
		c[i] <- Phi(muchat + xichat * deltachat[i]);
		r[i] <- Phi(murhat + xirhat * deltarhat[i]);
		u[i] <- Phi(muuhat + xiuhat * deltauhat[i]);
		
		// MPT Category Probabilities for Word Pairs
		theta[i,1] <- c[i] * r[i];
		theta[i,2] <- (1 - c[i]) * (u[i]) ^ 2;
		theta[i,3] <- (1 - c[i]) * 2 * u[i] * (1 - u[i]);
		theta[i,4] <- c[i] * (1 - r[i]) + (1 - c[i]) * (1 - u[i]) ^ 2;
	}
  for (i1 in 1:nparams)
    for (i2 in 1:nparams)
      rho[i1,i2] <- Sigma[i1,i2] / sqrt(Sigma[i1,i1] * Sigma[i2,i2]);
}
model {
  // Priors
	muchat ~ normal(0, 1);
	murhat ~ normal(0, 1);
	muuhat ~ normal(0, 1);
	
	Sigma ~ wishart(df, I);
  // Individual Effects
  deltahat ~ multi_normal(mudeltahat, Sigma);
  // Data
	for (i in 1:nsubjs)
		k[i] ~ multinomial(theta[i]);
}
generated quantities {
  real<lower=0,upper=1> muc;
  real<lower=0,upper=1> mur;
  real<lower=0,upper=1> muu;
  real sigmac;
  real sigmar;
  real sigmau;

  // Post-Processing Means, Standard Deviations, Correlations
	muc <- Phi(muchat);
	mur <- Phi(murhat);
	muu <- Phi(muuhat);
  sigmac <- xichat * sqrt(Sigma[1,1]);
  sigmar <- xirhat * sqrt(Sigma[2,2]);
  sigmau <- xiuhat * sqrt(Sigma[3,3]);
}"

### Riefer et al (2002) data:
load("dataMPT.Rdata")

nparams <- 3			# Number of free parameters per participant: c_i, r_i, u_i 
I <- diag(3)			# Identity matrix for Wishart

myinits <- list(
  list(deltahat=matrix(rnorm(21 * 3), 21, 3), Sigma=diag(3),
       muchat=rnorm(1), murhat=rnorm(1), muuhat=rnorm(1),
       xichat=runif(1), xirhat=runif(1), xiuhat=runif(1)),
  list(deltahat=matrix(rnorm(21 * 3), 21, 3), Sigma=diag(3), 
       muchat=rnorm(1), murhat=rnorm(1), muuhat=rnorm(1),
       xichat=runif(1), xirhat=runif(1), xiuhat=runif(1)),
  list(deltahat=matrix(rnorm(21 * 3), 21, 3), Sigma=diag(3), 
       muchat=rnorm(1), murhat=rnorm(1), muuhat=rnorm(1),
       xichat=runif(1), xirhat=runif(1), xiuhat=runif(1)))

# Parameters to be monitored
parameters <- c("muc", "mur", "muu", "sigmac", "sigmar", "sigmau", "rho")  

# Run higher iterations for better estimate
myiterations <- 2100 
mywarmup <- 100

k <- response_1
nsubjs <- nrow(k)    	# Number of word pairs per participant	
data <- list(k=k, nparams=nparams, nsubjs=nsubjs, I=I) # To be passed on to Stan

# The following command calls Stan with specific options.
# For a detailed description type "?stan".
samples_1 <- stan(model_code=model,   
                  data=data, 
                  init=myinits,  # If not specified, gives random inits
                  pars=parameters,
                  iter=myiterations, 
                  chains=3, 
                  thin=1,
                  warmup=mywarmup,  # Stands for burn-in; Default = iter/2
                  control = list(adapt_delta = 0.999, stepsize = 0.001, max_treedepth = 20)
                  # seed=123  # Setting seed; Default is random seed
)
samples_1
traceplot(samples_1, pars = c("muc", "mur", "muu", "rho", "sigmac", "sigmar", "sigmau", "lp__"))

k <- response_2
nsubjs <- nrow(k) 	 	# Number of word pairs per participant	
data <- list(k=k, nparams=nparams, nsubjs=nsubjs, I=I) # To be passed on to Stan

samples_2 <- stan(fit=samples_1,   
                  data=data, 
                  init=myinits,  # If not specified, gives random inits
                  pars=parameters,
                  iter=myiterations, 
                  chains=3, 
                  thin=1,
                  warmup=mywarmup,  # Stands for burn-in; Default = iter/2
                  control = list(adapt_delta = 0.999, stepsize = 0.001, max_treedepth = 20)
)
samples_2
traceplot(samples_2, pars = c("muc", "mur", "muu", "rho", "sigmac", "sigmar", "sigmau", "lp__"))

k <- response_6
nsubjs <- nrow(k) 	 	# Number of word pairs per participant	
data <- list(k=k, nparams=nparams, nsubjs=nsubjs, I=I) # To be passed on to Stan

samples_6 <- stan(fit=samples_1,   
                  data=data, 
                  init=myinits,  # If not specified, gives random inits
                  pars=parameters,
                  iter=myiterations, 
                  chains=3, 
                  thin=1,
                  warmup=mywarmup,  # Stands for burn-in; Default = iter/2
                  control = list(adapt_delta = 0.999, stepsize = 0.001, max_treedepth = 20)
)
samples_6
traceplot(samples_6, pars = c("muc", "mur", "muu", "rho", "sigmac", "sigmar", "sigmau", "lp__"))
# Now the values for the monitored parameters are in the "samples" object, 
# ready for inspection.

muc1 <- extract(samples_1)$muc
mur1 <- extract(samples_1)$mur
muu1 <- extract(samples_1)$muu
muc2 <- extract(samples_2)$muc
mur2 <- extract(samples_2)$mur
muu2 <- extract(samples_2)$muu
muc6 <- extract(samples_6)$muc
mur6 <- extract(samples_6)$mur
muu6 <- extract(samples_6)$muu

rhocr1 <- extract(samples_1)$rho[, 1, 2]
rhocu1 <- extract(samples_1)$rho[, 1, 3]
rhoru1 <- extract(samples_1)$rho[, 2, 3]
rhocr2 <- extract(samples_2)$rho[, 1, 2]
rhocu2 <- extract(samples_2)$rho[, 1, 3]
rhoru2 <- extract(samples_2)$rho[, 2, 3]
rhocr6 <- extract(samples_6)$rho[, 1, 2]
rhocu6 <- extract(samples_6)$rho[, 1, 3]
rhoru6 <- extract(samples_6)$rho[, 2, 3]

#### Plots posteriors of the group--level c, r, and u parameters
windows(10, 5)
layout(matrix(1:3, 1, 3, byrow=TRUE))
par(cex=1.1, mar=c(2, 2, 1, 1), mgp=c(.8, .1, 0))

plot(density(muc6), xlim=c(0, 1), ylim=c(0, 7), lty="dotted",
     ylab="Probability Density", xlab=expression(mu[c]), main="", 
     yaxt="n", xaxt="n")
lines(density(muc2), lty="dashed")
lines(density(muc1))
axis(1, seq(0, 1, by=.2), tick=FALSE)

plot(density(mur6), xlim=c(0, 1), ylim=c(0, 15), lty="dotted", ylab="",
     xlab=expression(mu[r]), yaxt="n", xaxt="n", main="")
lines(density(mur2), lty="dashed")
lines(density(mur1))
axis(1, seq(0, 1, by=.2), tick=FALSE)
legend(0, 15.5, c("Trial 1", "Trial 2", "Trial 6"), lty = c(1, 2, 3), 
       col=c("black"), text.col = "black", bty = "n")

plot(density(muu6), xlim=c(0, 1), ylim=c(0, 9), lty="dotted", ylab="",
     xlab=expression(mu[u]), yaxt="n", xaxt="n", main="")
lines(density(muu2), lty="dashed")
lines(density(muu1))
axis(1, seq(0, 1, by=.2), tick=FALSE)

#### Plots posteriors for the correlations
windows(10, 5)
layout(matrix(1:3, 1, 3, byrow=TRUE))
par(cex=1.1, mar=c(2, 2, 1, 1), mgp=c(.8, .1, 0))

plot(density(rhocr6), xlim=c(-1, 1), ylim=c(0, 1.5), lty="dotted",
     ylab="Probability Density", xlab=expression(italic(rho[cr])), main="", 
     yaxt="n", xaxt="n")
lines(density(rhocr2), lty="dashed")
lines(density(rhocr1))
axis(1, seq(-1, 1, by=.5), tick=FALSE)

plot(density(rhocu6), xlim=c(-1, 1), ylim=c(0, 1.5), lty="dotted", ylab="",
     xlab=expression(italic(rho[cu])), yaxt="n", xaxt="n", main="")
lines(density(rhocu2), lty="dashed")
lines(density(rhocu1))
axis(1, seq(-1, 1, by=.5), tick=FALSE)
legend(-1, 1.55, c("Trial 1", "Trial 2", "Trial 6"), lty = c(1, 2, 3), 
       col=c("black"), text.col = "black", bty = "n")

plot(density(rhoru6), xlim=c(-1, 1), ylim=c(0, 1.5), lty="dotted", ylab="",
     xlab=expression(italic(rho[ru])), yaxt="n", xaxt="n", main="")
lines(density(rhoru2), lty="dashed")
lines(density(rhoru1))
axis(1, seq(-1, 1, by=.5), tick=FALSE)
