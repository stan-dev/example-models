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
		
	for (i in 1:nsubjs) {
    vector[nsubjs] deltachat;
    vector[nsubjs] deltarhat;
    vector[nsubjs] deltauhat;
    
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

trial_1 <- list(subjs=21, items=20, response=structure(
	.Data=c(2,4,4,10,2,1,3,14,2,2,5,11,6,0,4,10,1,0,4,15,1,0,2,17,1,2,4,13,4,1,
          6,9,5,1,4,10,1,0,9,10,5,0,3,12,0,1,6,13,1,5,7,7,1,1,4,14,2,2,3,13,2,
          1,5,12,2,0,6,12,1,0,5,14,2,1,8,9,3,0,2,15,1,2,3,14),
  .Dim = c(4, 21)))

trial_2 <- list(subjs=21, items=20, response=structure(
  .Data = c(7,5,3,5,5,2,3,10,6,2,7,5,9,4,2,5,2,2,7,9,1,3,3,13,5,0,5,10,7,3,4,
            6,7,3,6,4,4,1,10,5,9,1,2,8,3,1,6,10,3,5,9,3,2,0,6,12,8,0,3,9,3,2,
            7,8,7,1,5,7,2,1,6,11,5,3,5,7,5,0,6,9,6,2,2,10),
  .Dim = c(4, 21)))

trial_6 <- list(subjs=21, items=20, response=structure(
  .Data = c(14,3,1,2,12,3,1,4,18,0,1,1,15,3,0,2,7,1,10,2,3,6,11,0,8,4,3,5,17,
            1,1,1,13,4,3,0,11,6,1,2,16,1,2,1,10,1,3,6,7,13,0,0,8,4,3,5,16,1,
            1,2,5,4,7,4,15,0,5,0,6,3,6,5,17,2,0,1,17,1,0,2,8,3,6,3),
  .Dim = c(4, 21)))
	   
response_1 <- t(trial_1$response)
response_2 <- t(trial_2$response)
response_6 <- t(trial_6$response)

nparams <- 3			# Number of free parameters per participant: c_i, r_i, u_i 
I <- diag(3)			# Identity matrix for Wishart

k <- response_1
nsubjs <- nrow(k) 	 	# Number of word pairs per participant	

# To be passed on to Stan
data <- list(k=k, nparams=nparams, nsubjs=nsubjs, I=I)  

myinits <- list(
	 list(deltahat=matrix(0, 21, 3), Sigma=diag(3),
        muchat=0, murhat=0, muuhat=0, xichat=1, xirhat=1, xiuhat=1))
       
# Parameters to be monitored
parameters <- c("muc", "mur", "muu", "sigmac", "sigmar", "sigmau", "rho")  

# Run higher iterations for better estimate
myiterations <- 3300 
mywarmup <- 300

# The following command calls Stan with specific options.
# For a detailed description type "?stan".
samples_1 <- stan(model_code=model,   
                  data=data, 
                  init=myinits,  # If not specified, gives random inits
                  pars=parameters,
                  iter=myiterations, 
                  chains=1, 
                  thin=1,
                  warmup=mywarmup,  # Stands for burn-in; Default = iter/2
                  # seed=123  # Setting seed; Default is random seed
)

k <- response_2
nsubjs <- nrow(k) 	 	# Number of word pairs per participant	
# To be passed on to Stan
data <- list(k=k, nparams=nparams, nsubjs=nsubjs, I=I)  
samples_2 <- stan(fit=samples_1,   
                  data=data, 
                  init=myinits,  # If not specified, gives random inits
                  pars=parameters,
                  iter=myiterations, 
                  chains=1, 
                  thin=1,
                  warmup=mywarmup,  # Stands for burn-in; Default = iter/2
                  # seed=123  # Setting seed; Default is random seed
)

k <- response_6
nsubjs <- nrow(k) 	 	# Number of word pairs per participant	
# To be passed on to Stan
data <- list(k=k, nparams=nparams, nsubjs=nsubjs, I=I) 
samples_6 <- stan(fit=samples_1,   
                  data=data, 
                  init=myinits,  # If not specified, gives random inits
                  pars=parameters,
                  iter=myiterations, 
                  chains=1, 
                  thin=1,
                  warmup=mywarmup,  # Stands for burn-in; Default = iter/2
                  # seed=123  # Setting seed; Default is random seed
)
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

#### Plots posteriors of the group--level c, r, and u parameters
windows(14,7)
layout(matrix(1:3,1,3,byrow=T))
par(mar = c(4.25, 3, 1, 1))
par(cex.axis=0.9)
plot(density((muc1)),ylim = c(0,6),xlim=c(0,1), axes=F,xlab="c",ylab="",
     main="")
axis(1,at=seq(0,1,0.2))
axis(2,at = c(0,12),labels=F,lwd.ticks=0)
mtext("Density",2,1,cex=0.55,las=0)
par(cex=0.60)
legend(0.3, 10, c("Trial1","Trial 2","Trial 3"),lty = c(1,2,3),col=c("black"),
       text.col = "black")
par(cex=0.65)
lines(density((muc2)),lty=2)
lines(density((muc6)),lty=3)

plot(density((mur1),to=0.75),ylim = c(0,13),xlim=c(0,1), axes=F,xlab="r",
     ylab="",main="")
axis(1,at=seq(0,1,0.2))
axis(2,at = c(0,13),labels=F,lwd.ticks=0)
mtext("Density",2,1,cex=0.55,las=0)
lines(density((mur2)),lty=2)
lines(density((mur6)),lty=3)

plot(density((muu1)),ylim = c(0,9),xlim=c(0,1), axes=F,xlab="u",ylab="",
     main="")
axis(1,at=seq(0,1,0.2))
axis(2,at = c(0,10),labels=F,lwd.ticks=0)
mtext("Density",2,1,cex=0.55,las=0)
lines(density((muu2)),lty=2)
lines(density((muu6)),lty=3)