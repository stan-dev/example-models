## 4. Introduction to random effects: Conventional Poisson GLMM for
## count data
## 4.3. Mixed models with random effects for variability among groups
## (site and year effects)
## 4.3.1. Generation and analysis of simulated data

library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
set.seed(123)

## Read data
## The data generation code is in bpa-code.txt, available at
## http://www.vogelwarte.ch/de/projekte/publikationen/bpa/complete-code-and-data-files-of-the-book.html
stan_data <- read_rdump("GLMM_Poisson2.data.R")

## Parameters monitored
params <- c("mu", "alpha", "beta", "sd_alpha", "sd_year")

# MCMC settings
ni <- 40000
nt <- 10
nb <- 30000
nc <- 4

## Initial values
inits <- lapply(1:nc, function(i) {
    list(mu = runif(1, 0, 2),
         alpha = runif(stan_data$nsite, -1, 1),
         beta = runif(3, -1, 1),
         sd_alpha = runif(1, 0, 0.1),
         sd_year = runif(1, 0, 0.1))})

## Call Stan from R
out <- stan("GLMM_Poisson2.stan", data = stan_data,
            init = inits, pars = params,
            chains = nc, iter = ni, warmup = nb, thin = nt,
            seed = 1,
            open_progress = FALSE)

## Note: this model converges very slowly, and there may be
## transitions after warmup that exceed the maximum treedepth.

## Summarize posteriors
print(out)
