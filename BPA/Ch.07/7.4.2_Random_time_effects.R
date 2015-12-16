## 7. Estimation of survival probabilities using capture-recapture data
## 7.4. Models with time-variation
## 7.4.2. Random time effects

library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
set.seed(123)

## Read data
## The data generation code is in bpa-code.txt, available at
## http://www.vogelwarte.ch/de/projekte/publikationen/bpa/complete-code-and-data-files-of-the-book.html
stan_data <- read_rdump("cjs_temp_raneff.data.R")

## Initial values
inits <- function() list(mean_phi = runif(1, 0, 1),
                         sigma = runif(1, 0, 10),
                         mean_p = runif(1, 0, 1))

## Parameters monitored
params <- c("mean_phi", "mean_p", "sigma2")

## MCMC settings
ni <- 4000
nt <- 3
nb <- 1000
nc <- 4

## Call Stan from R
cjs_temp_raneff <- stan("cjs_temp_raneff.stan",
                        data = stan_data, init = inits, pars = params,
                        chains = nc, iter = ni, warmup = nb, thin = nt,
                        seed = 1,
                        open_progress = FALSE)

## Summarize posteriors
print(cjs_temp_raneff, digits = 3)
