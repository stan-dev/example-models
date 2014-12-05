# clears workspace: 
rm(list=ls()) 

library(rstan)

###############################################################################
# Model 1 - is very simple, but is returning only parameter Phi 
# Model 2 - model is more complicated and is able to generate parameters z 
###############################################################################
# Your choice: ################################################################
choice <- 2    ################################################################
###############################################################################

if (choice == 1) {
model <- "
# Exam Scores
data { 
	int<lower=0> p;
	int<lower=0> k[p];
	int<lower=0> n;
}
transformed data {
	real psi;

	psi <- .5;
}
parameters {
	real<lower=.5, upper=1> phi;
	vector<lower=0, upper=1>[p] mix;
} 
model {
	phi ~ beta(1, 1)T[.5, 1];

	for (i in 1:p) {
		increment_log_prob(log_sum_exp(log(mix[i]) + binomial_log(k[i], n, phi),
									   log1m(mix[i]) + binomial_log(k[i], n, psi)));
    }
}"
} else if (choice == 2) {
model <- "
# Exam Scores
data { 
    int<lower=0> p;
    int<lower=0> k[p];
    int<lower=0> n;
}
transformed data {
	real psi;

	psi <- .5;
}
parameters {
    real<lower=.5, upper=1> phi;
    vector<lower=0, upper=1>[p] mix;
} 
transformed parameters {
    matrix[p, 2] lp_parts;

    for (i in 1:p) {
		lp_parts[i, 1] <- log(mix[i]) + binomial_log(k[i], n, phi);
		lp_parts[i, 2] <- log1m(mix[i]) + binomial_log(k[i], n, psi); 
    }
}
model {
    phi ~ beta(1, 1)T[.5, 1];

	for (i in 1:p) {
		increment_log_prob(log_sum_exp(lp_parts[i]));
	}
}
generated quantities {
	vector<lower=0, upper=1>[p] prob;
	int<lower=0, upper=1> z[p];
  
	for (i in 1:p) {
		prob[i] <- exp(lp_parts[i, 1]) / sum(exp(lp_parts[i]));
		z[i] <- bernoulli_rng(prob[i]);
	}
}"
}

k <- c(21, 17, 21, 18, 22, 31, 31, 34, 34, 35, 35, 36, 39, 36, 35)
p <- length(k)  # number of people
n <- 40  # number of questions

data <- list(p=p, k=k, n=n) # to be passed on to Stan

if (choice == 1) {
	myinits <- list(
		list(phi=.75, mix=rep(.5, p))) # Initial group assignment (for model 1)
	parameters <- c("phi")  # parameters to be monitored:	
} else if (choice == 2) {
	myinits <- list(
		list(phi=.75, mix=rep(.5, p))) # Initial group assignment (for model 2)
	parameters <- c("phi", "z")  # parameters to be monitored:	
}

# The following command calls Stan with specific options.
# For a detailed description type "?rstan".
samples <- stan(model_code=model,   
                data=data, 
                init=myinits,  # If not specified, gives random inits
                pars=parameters,
                iter=20000, 
                chains=1, 
                thin=1,
                # warmup = 100,  # Stands for burn-in; Don't set to 0 or low 
                # ... values, it can malfunction; Default = iter/2
                # seed = 123  # Setting seed; Default is random seed
)
# Now the values for the monitored parameters are in the "samples" object, 
# ready for inspection.

print(samples, digits_summary=3)
