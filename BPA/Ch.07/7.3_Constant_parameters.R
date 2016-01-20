## 7. Estimation of survival probabilities using capture-recapture data
## 7.3. Models with constant parameters

library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
set.seed(123)

## Read data
## The data generation code is in bpa-code.txt, available at
## http://www.vogelwarte.ch/de/projekte/publikationen/bpa/complete-code-and-data-files-of-the-book.html
stan_data <- read_rdump("cjs_c_c.data.R")

## Initial values
inits <- function() list(mean_phi = runif(1, 0, 1),
                         mean_p = runif(1, 0, 1))

## Parameters monitored
params <- c("mean_phi", "mean_p")

## MCMC settings
ni <- 2000
nt <- 1
nb <- 1000
nc <- 4

## Call Stan from R
cjs_c_c <- stan("cjs_c_c.stan",
                data = stan_data, init = inits, pars = params,
                chains = nc, iter = ni, warmup = nb, thin = nt,
                seed = 1,
                open_progress = FALSE)

## Summarize posteriors
print(cjs_c_c, digits = 3)
