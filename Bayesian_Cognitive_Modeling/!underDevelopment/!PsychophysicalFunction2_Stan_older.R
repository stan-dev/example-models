# clears workspace: 
rm(list=ls(all=TRUE)) 

library(rstan)



model <- "
# Logistic Psychophysical Function
data { 
    int nsubjs;
    int nstim[nsubjs];
    int n[nsubjs, 28];
    int r[nsubjs, 28];
    int x[nsubjs, 28];
    vector[nsubjs] xmean; 
}
parameters {
    real mua;
    real mub;
    real mup;
    real<lower=0,upper=1000> sigmaa;
    real<lower=0,upper=1000> sigmab;
    real<lower=0,upper=3> sigmap;
    vector[nsubjs] alpha;
    vector[nsubjs] beta;
    vector[nsubjs] probitphi;
    matrix<lower=0,upper=1>[nsubjs, 28] pi;
} 
model {
    real theta; 
    vector[nsubjs] probitphilim;
    vector[nsubjs] phi;

    # Priors
    mua ~ normal(0, sqrt(1000));
    mub ~ normal(0, sqrt(1000));
    mup ~ normal(0, 1); 

    alpha ~ normal(mua, sigmaa);
    beta ~ normal(mub, sigmab);
    probitphi ~ normal(mup, sigmap);
    
    for (i in 1:nsubjs) {
       
        if (probitphi[i] < -5)
            probitphilim[i] <- -5;
        else if (probitphi[i] > 5)
            probitphilim[i] <- 5;
        else 
            probitphilim[i] <- probitphi[i]; 

        phi[i] <- Phi(probitphilim[i]);
        
        for (j in 1:nstim[i]) {   
            theta <- 1 / (1 + exp(-(alpha[i] + beta[i] * (x[i, j] - xmean[i]))));
            r[i, j] ~ binomial(n[i, j], theta);

            increment_log_prob(
                log_sum_exp(log(phi[i]) + binomial_log(r[i, j], n[i, j], theta),
                            log(1 - phi[i]) + binomial_log(r[i, j], n[i, j], pi[i, j])));
        
        }
    }
}"


x <- as.matrix(read.table("data_x.txt", sep="\t"))
x[is.na(x)] = -5  # transforming because Stan can't handle NAs 

n <- as.matrix(read.table("data_n.txt", sep="\t"))
n[is.na(n)] = -5

r <- as.matrix(read.table("data_r.txt", sep="\t"))
r[is.na(r)] = -5

rprop <- as.matrix(read.table("data_rprop.txt", sep="\t"))

xmean <- c(318.888, 311.0417, 284.4444, 301.5909, 
           296.2000, 305.7692, 294.6429, 280.3571)
nstim <- c(27, 24, 27, 22, 25, 26, 28, 28)
nsubjs <- 8

# to be passed on to Stan
data <- list(x=x, xmean=xmean, n=n, r=r, nsubjs=nsubjs, nstim=nstim) 

# myinits <- list(  # Doesn't work with this initial list
#     list(alpha=runif(nsubjs, -2, 2), beta=runif(nsubjs, 0, .5), 
#          mua=0, mub=0, sigmaa=1, sigmab=1),
#     list(alpha=runif(nsubjs, -2, 2), beta=runif(nsubjs, 0, .5), 
#          mua=0, mub=0, sigmaa=1, sigmab=1),
#     list(alpha=runif(nsubjs, -2, 2), beta=runif(nsubjs, 0, .5), 
#          mua=0, mub=0, sigmaa=1, sigmab=1))

parameters <- c("alpha", "beta")  # Parameters to be monitored

# The following command calls Stan with specific options.
# For a detailed description type "?rstan".
samples <- stan(model_code=model,
                data=data, 
                init=0,  # !!! has to be set on 0 for this example, 
#                 pars=parameters,
                iter=500, 
                chains=3, 
                thin=1,
                # warmup = 100,  # Stands for burn-in; Don't set to 0 or 
                # ... low values, it can malfunction; Default = iter/2
                # seed = 123  # Setting seed; Default is random seed
)
# Now the values for the monitored parameters are in the "samples" object, 
# ready for inspection.

print(samples, digits_summary=2)
