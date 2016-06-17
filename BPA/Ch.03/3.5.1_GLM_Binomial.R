## 3. Introduction to the generalized linear model (GLM): The simplest
## model for count data
## 3.5. Binomial GLM for modeling bounded counts or proportions
## 3.5.1. Generation and analysis of simulated data

library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
set.seed(123)

## Read data
## The data generation code is in bpa-code.txt, available at
## http://www.vogelwarte.ch/de/projekte/publikationen/bpa/complete-code-and-data-files-of-the-book.html
stan_data <- read_rdump("GLM_Binomial.data.R")

## Initial values
inits <- function() list(alpha = runif(1, -1, 1),
                         beta1 = runif(1, -1, 1),
                         beta2 = runif(1, -1, 1))

## Parameters monitored
params <- c("alpha", "beta1", "beta2", "p")

## MCMC settings
ni <- 2000
nt <- 1
nb <- 1000
nc <- 4

## Call Stan from R
out <- stan("GLM_Binomial.stan", data = stan_data,
            init = inits, pars = params,
            chains = nc, thin = nt, iter = ni, warmup = nb,
            seed = 1,
            open_progress = FALSE)

## Summarize posteriors
print(out)
