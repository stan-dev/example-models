## 6. Estimation of the size of a closed population
## 6.2. Generation and analysis of simulated data with data
## augmentation
## 6.2.5. Combined effects: model Mth

library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
set.seed(1234)

## Read data
## The data generation code is in bpa-code.txt, available at
## http://www.vogelwarte.ch/de/projekte/publikationen/bpa/complete-code-and-data-files-of-the-book.html
source("Mth.data.R")

## Bundle data
stan_data <- list(y = y, M = M, T = T)

## Initial values
inits <- function() list(sigma = runif(1, 0.5, 0.9))

## Parameters monitored
params <- c("N", "mean_p", "mean_lp", "sigma", "omega")

## MCMC settings
ni <- 15000
nt <- 5
nb <- 10000
nc <- 4

## Call Stan from R
out <- stan("Mth.stan",
            data = stan_data, init = inits, pars = params,
            chains = nc, iter = ni, warmup = nb, thin = nt,
            seed = 1,
            open_progress = FALSE)

## Summarize posteriors
print(out, digits = 3)
