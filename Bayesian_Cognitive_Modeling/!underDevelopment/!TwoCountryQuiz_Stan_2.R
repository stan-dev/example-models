# clears workspace: 
rm(list=ls()) 

library(rstan)

model <- "

data { 
    int nx;
	int nz;
	int<lower=0, upper=1> k[nx, nz];
}
parameters {
	real<lower=0, upper=1> alpha;     # Match
	real<lower=0, upper=alpha> beta;  # Mismatch
	real<lower=0, upper=1> x_prob[nx];
	real<lower=0, upper=1> z_prob[nz];
} 
transformed parameters {
	vector[4] lp_parts[nx, nz];
	for (i in 1:nx) {
		for (j in 1:nz) {
			lp_parts[i, j, 1] <- log1m(x_prob[i]) + log1m(z_prob[j]) + 
									bernoulli_log(k[i, j], alpha);
			lp_parts[i, j, 2] <- log(x_prob[i]) + log(z_prob[j]) + 
									bernoulli_log(k[i, j], alpha);
			lp_parts[i, j, 3] <- log1m(x_prob[i]) + log(z_prob[j]) + 
									bernoulli_log(k[i, j], beta);
			lp_parts[i, j, 4] <- log(x_prob[i]) + log1m(z_prob[j]) + 
									bernoulli_log(k[i, j], beta);
    }
  }
}
model {
	for (i in 1:nx) {
		for (j in 1:nz) {
			increment_log_prob(log_sum_exp(lp_parts[i, j]));
		}
	}
}
generated quantities {
	int x[nx];
	int z[nz];
	real prob_x[nx];
	real prob_z[nz];
	
	for (i in 1:nx) {
		for (j in 1:nz) {  
			prob_x[i] <- (sum(exp(lp_parts[2, i])) + sum(exp(lp_parts[4, i]))) / 
						 (sum(exp(lp_parts[1, i])) + sum(exp(lp_parts[2, i])) +
						  sum(exp(lp_parts[3, i])) + sum(exp(lp_parts[4, i])));
		}
		x[i] <- bernoulli_rng(prob_x[i]);
	}
}"


dset <- 1

if (dset==1) {
	k <- c(1,0,0,1,1,0,0,1,
		   1,0,0,1,1,0,0,1,
		   0,1,1,0,0,1,0,0,
		   0,1,1,0,0,1,1,0,
		   1,0,0,1,1,0,0,1,
		   0,0,0,1,1,0,0,1,
		   0,1,0,0,0,1,1,0,
		   0,1,1,1,0,1,1,0)
	k <- matrix(k, nrow=8, byrow=T)
}


nx <- nrow(k)
nz <- ncol(k)

data <- list(nx=nx, nz=nz, k=k) # To be passed on to Stan

myinits <-list(
  list(z_prob=runif(nz), x_prob=runif(nx), alpha=0.75, beta=0.25))
  
parameters <- c("alpha", "beta", "x")  # Parameters to be monitored
# parameters <- c("z", "x", "alpha", "beta")  # Parameters to be monitored

# The following command calls Stan with specific options.
# For a detailed description type "?rstan".
samples <- stan(model_code=model,   
                data=data, 
               init=myinits,  # If not specified, gives random inits
                pars=parameters,
                iter=10000, 
                chains=1, 
                thin=1,
                # warmup = 100,  # Stands for burn-in; Don't set to 0 or 
				# ... low values, it can malfunction; Default = iter/2
                # seed = 123  # Setting seed; Default is random seed
)
# Now the values for the monitored parameters are in the "samples" object, 
# ready for inspection.

print(samples, digits=3)
