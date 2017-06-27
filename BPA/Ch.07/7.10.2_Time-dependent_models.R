## 7. Estimation of survival probabilities using capture-recapture data
## 7.10. Fitting the CJS to data in the m-array format: the multinomial likelihood
## 7.10.2. Time-dependent models

library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
set.seed(123)

## Read data
## The data generation code is in bpa-code.txt, available at
## http://www.vogelwarte.ch/de/projekte/publikationen/bpa/complete-code-and-data-files-of-the-book.html
stan_data <- read_rdump("cjs_mnl.data.R")

## Parameters monitored
params <- c("phi", "p", "fit", "fit_new")

## MCMC settings
ni <- 2000
nt <- 1
nb <- 1000
nc <- 4

## Initial values
inits <- function() list(phi= runif(dim(stan_data$marr)[2] - 1, 0, 1),
                         p = runif(dim(stan_data$marr)[2] - 1, 0, 1))

## Call Stan from R
cjs_mnl_age  <- stan("cjs_mnl.stan",
                     data = stan_data, init = inits, pars = params,
                     chains = nc, iter = ni, warmup = nb, thin = nt,
                     seed = 1,
                     open_progress = FALSE)

## Summarize posteriors
print(cjs_mnl_age, digits = 3)
