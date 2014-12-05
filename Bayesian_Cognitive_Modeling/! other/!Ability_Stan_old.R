# clears workspace: 
rm(list=ls(all=TRUE)) 

library(rstan)

model <- "
# Ability Correlation for ESP Replication
data { 
    int<lower=1> nsubjs;
    int<lower=1> ntrials;
    int<lower=0> k[nsubjs, 2];
}
parameters {
    vector[2] mu;
    vector<lower=0>[2] lambda;
    vector[2] thetap[nsubjs];
    real<lower=0, upper=1> r;
} 
transformed parameters {
    vector<lower=0>[2] sigma;
    vector[2] theta[nsubjs];
    cov_matrix[2] T;
    
    # Reparametrization
    sigma[1] <- 1 / (sqrt(lambda[1]));
    sigma[2] <- 1 / (sqrt(lambda[2]));

    T[1, 1] <- 1 / lambda[1];
    T[1, 2] <- r * sigma[1] * sigma[2];
    T[2, 1] <- r * sigma[1] * sigma[2];
    T[2, 2] <- 1 / lambda[2];

    for (i in 1:nsubjs) {
        for (j in 1:2) {
            theta[i, j] <- Phi(thetap[i, j]);
        }
    }
}
model {
    # Priors
    mu ~ normal(0, sqrt(1000));
    lambda ~ gamma(.001, .001);

    # Data
    for (i in 1:nsubjs) {
        thetap[i] ~ multi_normal(mu, T);
        for (j in 1:2) {
            k[i, j] ~ binomial(ntrials, theta[i, j]);
        }
    }
}"


# Proportion correct on erotic pictures, block 1 and block 2:
prc1.ero <- c(0.6000000, 0.5333333, 0.6000000, 0.6000000, 0.4666667, 
             0.6666667, 0.6666667, 0.4000000, 0.6000000, 0.6000000,
             0.4666667, 0.6666667, 0.4666667, 0.6000000, 0.3333333,
             0.4000000, 0.4000000, 0.2666667, 0.3333333, 0.5333333,
             0.6666667, 0.5333333, 0.6000000, 0.4000000, 0.4666667, 
             0.7333333, 0.6666667, 0.6000000, 0.6666667, 0.5333333,
             0.5333333, 0.6666667, 0.4666667, 0.3333333, 0.4000000,
             0.5333333, 0.4000000, 0.4000000, 0.3333333, 0.4666667,
             0.4000000, 0.4666667, 0.4666667, 0.5333333, 0.3333333,
             0.7333333, 0.2666667, 0.6000000, 0.5333333, 0.4666667,
             0.4000000, 0.5333333, 0.6666667, 0.4666667, 0.5333333,
             0.5333333, 0.4666667, 0.4000000, 0.4666667, 0.6666667,
             0.4666667, 0.3333333, 0.3333333, 0.3333333, 0.4000000,
             0.4000000, 0.6000000, 0.4666667, 0.3333333, 0.3333333,
             0.6666667, 0.5333333, 0.3333333, 0.6000000, 0.4666667,
             0.4666667, 0.4000000, 0.3333333, 0.4666667, 0.5333333,
             0.8000000, 0.4000000, 0.5333333, 0.5333333, 0.6666667,
             0.6666667, 0.6666667, 0.6000000, 0.6000000, 0.5333333,
             0.3333333, 0.4666667, 0.6666667, 0.5333333, 0.3333333,
             0.3333333, 0.2666667, 0.2666667, 0.4666667, 0.6666667)

prc2.ero <- c(0.3333333, 0.6000000, 0.5333333, 0.2666667, 0.6666667,
             0.5333333, 0.6666667, 0.4666667, 0.4666667, 0.6666667,
             0.4000000, 0.6666667, 0.2666667, 0.4000000, 0.4666667,
             0.3333333, 0.5333333, 0.6000000, 0.3333333, 0.4000000,
             0.4666667, 0.4666667, 0.6000000, 0.5333333, 0.5333333,
             0.6000000, 0.5333333, 0.6666667, 0.6000000, 0.2666667,
             0.4666667, 0.4000000, 0.6000000, 0.5333333, 0.4000000,
             0.4666667, 0.5333333, 0.3333333, 0.4000000, 0.4666667,
             0.8000000, 0.6000000, 0.2000000, 0.6000000, 0.4000000,
             0.4000000, 0.2666667, 0.2666667, 0.6000000, 0.4000000,
             0.4000000, 0.4000000, 0.4000000, 0.4000000, 0.6666667,
             0.7333333, 0.5333333, 0.5333333, 0.3333333, 0.6000000,
             0.5333333, 0.5333333, 0.4666667, 0.5333333, 0.4666667,
             0.5333333, 0.4000000, 0.4000000, 0.4666667, 0.6000000,
             0.6000000, 0.6000000, 0.4666667, 0.6000000, 0.6666667,
             0.5333333, 0.4666667, 0.6000000, 0.2000000, 0.5333333,
             0.4666667, 0.4000000, 0.5333333, 0.5333333, 0.5333333,
             0.5333333, 0.6000000, 0.6666667, 0.4000000, 0.4000000,
             0.5333333, 0.8000000, 0.6000000, 0.4000000, 0.2000000,
             0.6000000, 0.6666667, 0.4666667, 0.4666667, 0.4666667)             

k <- matrix(cbind(as.integer(round(prc1.ero * 60)), 
                  as.integer(round(prc2.ero * 60))), nrow=100)
nsubjs <- nrow(k) # number of participants
ntrials <- 60


data <- list(k=k, nsubjs=nsubjs, ntrials=ntrials) # To be passed on to Stan

# myinits <- list(
#     list(r=0, mu=c(0, 0), lambda=c(1, 1), 
#          thetap=matrix(rnorm(200), 100, 2)))

parameters <- c("r", "mu", "sigma")  # Parameters to be monitored

# The following command calls Stan with specific options.
# For a detailed description type "?rstan".
samples <- stan(model_code=model,   
                data=data, 
                init=myinits,  # If not specified, gives random inits
                pars=parameters,
                iter=1000, 
                chains=1, 
                thin=1,
                # warmup = 100,  # Stands for burn-in; Don't set to 0 or 
                # ... low values, it can malfunction; Default = iter/2
                # seed = 123  # Setting seed; Default is random seed
)
# Now the values for the monitored parameters are in the "samples" object, 
# ready for inspection.

print(samples, digits_summary=2)
