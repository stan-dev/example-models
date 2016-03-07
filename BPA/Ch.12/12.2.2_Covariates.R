## 12. Metapopulation modeling of abundance using
##     hierarchical Poisson regression: binomial mixture models
## 12.2. Generation and analysis of simulated data
## 12.2.2. Introducing covariates

library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
set.seed(123)

## Read data
## The data generation code is in bpa-code.txt, available at
## http://www.vogelwarte.ch/de/projekte/publikationen/bpa/complete-code-and-data-files-of-the-book.html
source("binmix_cov.data.R")

## Parameters monitored
params <- c("totalN", "alpha0", "alpha1", "beta0", "beta1")

## MCMC settings
ni <- 1000
nt <- 1
nb <- 500
nc <- 4

## Initial values
inits <- lapply(1:nc, function(i)
    list(alpha0 = runif(1, -1, 1),
         alpha1 = runif(1, -1, 1),
         beta0 = runif(1, -1, 1),
         beta1 = runif(1, -1, 1)))

## Call Stan from R
out <- stan("binmix_cov.stan",
            data = list(y = y, R = R, T = T, K = 100),
            init = inits, pars = params,
            chains = nc, iter = ni, warmup = nb, thin = nt,
            seed = 1,
            open_progress = FALSE)
print(out, digits = 3)
