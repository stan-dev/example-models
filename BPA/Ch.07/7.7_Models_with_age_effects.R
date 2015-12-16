## 7. Estimation of survival probabilities using capture-recapture data
## 7.7. Models with age effects

library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
set.seed(123)

## Read data
## The data generation code is in bpa-code.txt, available at
## http://www.vogelwarte.ch/de/projekte/publikationen/bpa/complete-code-and-data-files-of-the-book.html
stan_data <- read_rdump("cjs_age.data.R")

## Parameters monitored
params <- c("beta", "mean_p")

## MCMC settings
ni <- 2000
nt <- 1
nb <- 1000
nc <- 4

## Initial values
inits <- function() list(beta = runif(2, 0, 1),
                         mean_p = runif(1, 0, 1))

## Call Stan from R
cjs_age  <- stan("cjs_age.stan",
                 data = stan_data, init = inits, pars = params,
                 chains = nc, iter = ni, warmup = nb, thin = nt,
                 seed = 1,
                 open_progress = FALSE)

## Summarize posteriors
print(cjs_age, digits = 3)
