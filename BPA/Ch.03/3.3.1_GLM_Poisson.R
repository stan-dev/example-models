## 3. Introduction to the generalized linear model (GLM): The simplest
## model for count data
## 3.3. Poisson GLM in R and WinBUGS for modeling times series of
## counts
## 3.3.1. Generation and analysis of simulated data

library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
set.seed(123)

## Read data
## The data generation code is in bpa-code.txt, available at
## http://www.vogelwarte.ch/de/projekte/publikationen/bpa/complete-code-and-data-files-of-the-book.html
source("GLM_Poisson.data.R")

## Bundle data
mean.year <- mean(data$year)             # Mean of year covariate
sd.year <- sd(data$year)                 # SD of year covariate

stan_data <- list(C = C, n = n, year = year)

## Initial values
inits <- function() list(alpha = runif(1, -2, 2),
                         beta1 = runif(1, -3, 3))

## Parameters monitored
params <- c("alpha", "beta1", "beta2", "beta3", "lambda")

## MCMC settings
ni <- 2000
nt <- 1
nb <- 1000
nc <- 4

## Call Stan from R
out <- stan("GLM_Poisson.stan", data = stan_data,
            init = inits, pars = params,
            chains = nc, thin = nt, iter = ni, warmup = nb,
            seed = 1,
            open_progress = FALSE)

## Summarize posteriors
print(out)
