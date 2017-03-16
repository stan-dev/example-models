## 8. Estimation of survival probabilities using mark-recovery data
## 8.4. Real data example: age-dependent survival in Swiss red kites

library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
set.seed(123)

## Read data
## The data generation code is in bpa-code.txt, available at
## http://www.vogelwarte.ch/de/projekte/publikationen/bpa/complete-code-and-data-files-of-the-book.html
stan_data <- read_rdump("mr_mnl_age3.data.R")

## Initial values
inits <- function() list(sjuv = runif(1, 0, 1),
                         ssub = runif(1, 0, 1),
                         sad = runif(1, 0, 1),
                         rjuv = runif(1, 0, 1),
                         rad = runif(1, 0, 1))

## Parameters monitored
params <- c("sjuv", "ssub", "sad", "rjuv", "rad")

## MCMC settings
ni <- 4000
nt <- 2
nb <- 2000
nc <- 4

## Call Stan from R
rk_ageA <- stan("mr_mnl_age3.stan",
                data = stan_data, init = inits, pars = params,
                chains = nc, iter = ni, warmup = nb, thin = nt,
                seed = 1,
                open_progress = FALSE)

## Summarize posteriors
print(rk_ageA, digits = 3)
