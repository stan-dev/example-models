# clears workspace: 
rm(list=ls()) 

library(rstan)

model <- "
// BART Model of Risky Decision-Making
data { 
  int<lower=1> ntrials;
  real<lower=0,upper=1> p;
  int<lower=1> options[ntrials];
  int d[ntrials,30];
}
parameters {
  // Priors
  real<lower=0,upper=10> gplus;
  real<lower=0,upper=10> beta;
} 
transformed parameters {
  real<lower=0> omega;

  // Optimal Number of Pumps
  omega <- -gplus / log1m(p);  // log1m(p) equals log(1 - p), but faster
}
model {
  // Choice Data
  for (j in 1:ntrials) {
    for (k in 1:options[j]) {
      real theta;
      theta <- 1 - inv_logit(-beta * (k - omega));
      d[j,k] ~ bernoulli(theta);
    }
  }
}"

p       <- .15	# (Belief of) bursting probability
ntrials <- 90   # Number of trials for the BART

Data   <- matrix(data=as.numeric(as.matrix(read.table("GeorgeSober.txt"))[-1, ]),
                 ntrials, 8)
d      <- matrix(-99, ntrials, 30) # Data in binary format

cash    <-(Data[, 7] != 0) * 1	# Cash or burst?
npumps  <- Data[, 6]				    # Nr. of pumps
options <- cash + npumps        # Nr. of decision possibilities

for (j in 1:ntrials) {
	if (npumps[j] > 0) {d[j, 1:npumps[j]] <- rep(0, npumps[j])}
	if (cash[j] == 1) {d[j, (npumps[j] + 1)] <- 1}
}

# To be passed on to Stan:
data <- list(ntrials=ntrials, p=p, options=options, d=d) 

myinits <- list(
  list(gplus=1.2, beta=.5))

parameters <- c("gplus", "beta")  # Parameters to be monitored

# The following command calls Stan with specific options.
# For a detailed description type "?rstan".
samples <- stan(model_code=model,   
                data=data, 
                init=myinits,  # If not specified, gives random inits
                pars=parameters,
                iter=5000, 
                chains=1, 
                thin=1,
                warmup=2000,  # Stands for burn-in; Default = iter/2
                # seed=123  # Setting seed; Default is random seed
)
# Now the values for the monitored parameters are in the "samples" object, 
# ready for inspection.

gplus <- extract(samples)$gplus
beta  <- extract(samples)$beta

#################### PLOT RESULTS

par (cex.main = 2.5, cex.lab = 2, cex.axis = 1.5, mar = c(5, 5, 4, 0), las = 1)
layout (matrix (1:3, 1, 3))

hist(npumps, main = " ", xlab="Number of Pumps", ylab="Frequency", breaks=c(0:7), 
     col="lightblue", axes=F)
axis(2)
axis(1, at=seq(.5,6.5,by=1), labels = c("1", "2", "3", "4", "5", "6", "7"))
     
plot(density(gplus), xlab = expression (gamma^'+'),
	main = " ", bty = 'n', lwd=2, ylab="Posterior Density")
plot(density(beta), xlab = expression (beta),
	main = " ", bty = 'n', lwd=2, lab = c(5,3,5), ylab="Posterior Density")

