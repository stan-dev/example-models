library(gpclib);
library(maptools);  
library(spdep);
gpclibPermit()

library(CARBayes);

## Run CARBayes on exactly the Carlin data
#The function fits a Gaussian random effects models to spatial data,
#where the random effects are modelled by BYM
#conditional autoregressive (CAR) model (Besag et. al. 1991).
#The model represents the mean function for the set of Gaussian responses
#by a combination of covariates and two sets of random effects.
#For the latter, the first set are independent, while the second are
#spatially correlated and come from the IAR model.
#A set of offsets can also be included on the linear predictor scale.
#Inference is based on Markov Chain Monte Carlo (MCMC) simulation,
#using a combination of Gibbs sampling and Metropolis steps.

source("carlin_data.R");
x = 0.1 * data$x;
y = data$y;
E = log(data$E);

source("scotland_nb.R");
W.mat <- nb2mat(scotland_nbs, style="B")

# run CARbayes version of BYM
formula = y ~ x + offset(E);

gc()
carbayes_fit10k = S.CARbym(formula = formula, family="poisson", W=W.mat,
                        burnin=10000, n.sample=20000,
                        prior.sigma2=c(3.2761, 1/1.81),
                        prior.tau2=c(1.0,1.0));

gc()
carbayes_fit100k = S.CARbym(formula = formula, family="poisson", W=W.mat,
                        burnin=100000, n.sample=200000,
                        prior.sigma2=c(3.2761, 1/1.81),
                        prior.tau2=c(1.0,1.0));

gc()
carbayes_fit100k_thin = S.CARbym(formula = formula, family="poisson", W=W.mat,
                        burnin=100000, n.sample=200000, thin=20,
                        prior.sigma2=c(3.2761, 1/1.81),
                        prior.tau2=c(1.0,1.0));

gc()
carbayes_fit1M_thin = S.CARbym(formula = formula, family="poisson", W=W.mat,
                        burnin=1000000, n.sample=2000000, thin=20,
                        prior.sigma2=c(3.2761, 1/1.81),
                        prior.tau2=c(1.0,1.0));

gc()
carbayes_fit10k$summary.results
gc()
carbayes_fit100k$summary.results
gc()
carbayes_fit100k_thin$summary.results
gc()
carbayes_fit1M_thin$summary.results

