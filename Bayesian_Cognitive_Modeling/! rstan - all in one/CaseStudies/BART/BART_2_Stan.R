# clears workspace: 
rm(list=ls()) 

library(rstan)

model <- "
// BART Model of Risky Decision-Making
data {
  int<lower=1> nconds;
  int<lower=1> ntrials;
  real<lower=0,upper=1> p;
  int<lower=1> options[nconds,ntrials];
  int d[nconds,ntrials,30];
}
parameters {
  vector<lower=0>[nconds] gplus;
  vector<lower=0>[nconds] beta;
  
  // Priors
  real<lower=0,upper=10> mug;
  real<lower=0,upper=10> sigmag;
  real<lower=0,upper=10> mub;
  real<lower=0,upper=10> sigmab;
} 
transformed parameters {
  vector<lower=0>[nconds] omega;

  // Optimal Number of Pumps
  omega <- -gplus / log1m(p);
}
model {
  for (i in 1:nconds) {
    gplus[i] ~ normal(mug, sigmag)T[0,];
    beta[i] ~ normal(mub, sigmab)T[0,];
  }
  // Choice Data
  for (i in 1:nconds) {
    for (j in 1:ntrials) {
      for (k in 1:options[i,j]) {
        real theta;
        theta <- 1 - inv_logit(-beta[i] * (k - omega[i]));
        d[i,j,k] ~ bernoulli(theta);
      }
    }
  }
}"

p       <- .15  # (Belief of) bursting probability
ntrials <- 90  # Number of trials for the BART

################### READ IN THE DATA

Data <- list(
  matrix(data=as.numeric(as.matrix(read.table("GeorgeSober.txt"))[-1, ]),
         ntrials, 8),
  matrix(data=as.numeric(as.matrix(read.table("GeorgeTipsy.txt"))[-1, ]),
         ntrials, 8),
  matrix(data=as.numeric(as.matrix(read.table("GeorgeDrunk.txt"))[-1, ]),
         ntrials, 8)
)

nconds <- length(Data)                         
cash   <- npumps <- matrix (, nconds, ntrials)  # Cashes and nr. of pumps
d      <- array (-99, c(nconds, ntrials, 30))      # Data in binary format

for (i in 1:nconds) {
  cash[i, ]   <- (Data[[i]][, 7] != 0) * 1  # Cash or burst?
  npumps[i, ] <- Data[[i]][, 6]             # Nr. of pumps

  for (j in 1:ntrials) {
    if (npumps[i, j] > 0) {d[i, j, 1:npumps[i, j]] <- rep(0, npumps[i, j])}
    if (cash[i, j] == 1) {d[i, j, (npumps[i, j] + 1)] <- 1}
  }
}
options <- cash + npumps  # Nr. of decision possibilities

# To be passed on to Stan:
data <- list(nconds=nconds, ntrials=ntrials, p=p, options=options, d=d) 

myinits <- list(
  list(gplus=rep(1.2, 3), beta=rep(.5, 3), mug=1.2, sigmag=.1, mub=.8, sigmab=.8),
  list(gplus=rep(1.2, 3), beta=rep(.5, 3), mug=1.5, sigmag=.2, mub=1, sigmab=1.2))

# Parameters to be monitored:
parameters <- c("beta", "gplus", "mub", "mug", "sigmab", "sigmag")  

# The following command calls Stan with specific options.
# For a detailed description type "?rstan".
samples <- stan(model_code=model,   
                data=data, 
                init=myinits,  # If not specified, gives random inits
                pars=parameters,
                iter=1500, 
                chains=2, 
                thin=1,
                warmup=500,  # Stands for burn-in; Default = iter/2
                seed=1234  # Setting seed; Default is random seed
)
# Now the values for the monitored parameters are in the "samples" object, 
# ready for inspection.

gplus.sober <- extract(samples)$gplus[,1]
gplus.tipsy <- extract(samples)$gplus[,2]
gplus.drunk <- extract(samples)$gplus[,3]

beta.sober  <- extract(samples)$beta[,1]
beta.tipsy  <- extract(samples)$beta[,2]
beta.drunk  <- extract(samples)$beta[,3]

#################### PLOT SOME RESULTS
windows()
par (cex.main = 2.5, cex.lab = 2, cex.axis = 1.5, mar = c(5, 5, 4, 0), las = 1)
layout (matrix (1:9, 3, 3, byrow = T))

hist(npumps[1,], xlab=" ", main = "Sober: # pumps", breaks=c(0:max(npumps[1,])),
     xlim = range(c(0:9)), col="lightblue", axes=F)
axis(2)
axis(1, at=seq(.5,8.5,by=1), 
     labels = c("1", "2", "3", "4", "5", "6", "7", "8", "9"))

plot(density(gplus.sober), xlab = expression (gamma^'+'), xlim = c(0.6,1.8), 
     main = expression (paste ("Sober: Posterior ", gamma^'+')), bty = 'n')
plot (density(beta.sober), xlab = expression (beta), bty = 'n',
      main = expression (paste ("Sober: Posterior ", beta)), xlim = c(0.2,1.4))

hist(npumps[2,], xlab=" ", main = "Tipsy: # pumps", breaks=c(0:max(npumps[2,])),
     xlim = range(c(0:9)), col="lightblue", axes=F)
axis(2)
axis(1, at=seq(.5,8.5,by=1), 
     labels = c("1", "2", "3", "4", "5", "6", "7", "8", "9"))
plot(density (gplus.tipsy), xlab = expression (gamma^'+'), xlim = c(0.6,1.8),
     main = expression (paste ("Tipsy: Posterior ", gamma^'+')), bty = 'n')
plot(density (beta.tipsy), xlab = expression (beta), xlim = c(0.2,1.4),
     main = expression (paste ("Tipsy: Posterior ", beta)), bty = 'n')

hist(npumps[3,], xlab="Number of Pumps", main = "Drunk: # pumps",
     breaks=c(0:max(npumps[3,])), xlim = range(c(0:9)), col="lightblue", axes=F)
axis(2)
axis(1, at=seq(.5,8.5,by=1), 
     labels = c("1", "2", "3", "4", "5", "6", "7", "8", "9"))
plot(density(gplus.drunk), xlab = expression (gamma^'+'), xlim = c(0.6,1.8),
     main = expression(paste ("Drunk: Posterior ", gamma^'+')),  bty = 'n')
plot(density(beta.drunk), xlab = expression (beta), xlim = c(0.2,1.4),
     main = expression(paste ("Drunk: Posterior ", beta)),  bty = 'n')

