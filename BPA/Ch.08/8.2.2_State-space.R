## 8. Estimation of survival probabilities using mark-recovery data
## 8.2. The mark-recovery model as a state-space model
## 8.2.2. Analysis of a model with constant parameters

library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
set.seed(123)

## Read data
## The data generation code is in bpa-code.txt, available at
## http://www.vogelwarte.ch/de/projekte/publikationen/bpa/complete-code-and-data-files-of-the-book.html
stan_data <- read_rdump("mr_ss.data.R")

## Initial values
inits <- function() list(mean_s = runif(1, 0, 1),
                         mean_r = runif(1, 0, 1))

## Parameters monitored
params <- c("mean_s", "mean_r")

## MCMC settings
ni <- 2000
nt <- 1
nb <- 1000
nc <- 4

## Call Stan from R
mr_ss <- stan("mr_ss.stan",
              data = stan_data, init = inits, pars = params,
              chains = nc, iter = ni, warmup = nb, thin = nt,
              seed = 1,
              open_progress = FALSE)

## Summarize posteriors
print(mr_ss, digits = 3)
