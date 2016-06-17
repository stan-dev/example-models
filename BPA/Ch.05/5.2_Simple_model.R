## 5. State-space models
## 5.2. A simple model

library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
set.seed(123)

## Read data
## The data generation code is in bpa-code.txt, available at
## http://www.vogelwarte.ch/de/projekte/publikationen/bpa/complete-code-and-data-files-of-the-book.html
stan_data <- read_rdump("ssm.data.R")

## Parameters monitored
params <- c("lambda", "mean_lambda", "sigma2_obs", "sigma2_proc",
            "N_est")

## MCMC settings
ni <- 10000
nt <- 5
nb <- 5000
nc <- 4

## Initial values
inits <- lapply(1:nc, function(i) {
    list(sigma_proc = runif(1, 0, 5),
         mean_lambda = runif(1, 0.1, 2),
         sigma_obs = runif(1, 0, 10),
         N_est1 = runif(1, 20, 40))})

## Call Stan from R
ssm <- stan("ssm.stan",
            data = stan_data, init = inits, pars = params,
            chains = nc, iter = ni, warmup = nb, thin = nt,
            seed = 1,
            control = list(adapt_delta = 0.999),
            open_progress = FALSE)

## Note: there may be divergent transitions after warmup.

## Summarize posteriors
print(ssm)
