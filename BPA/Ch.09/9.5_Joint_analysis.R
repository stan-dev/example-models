## 9. Multistate capture-recapture models
## 9.5. Joint analysis of capture-recapture and mark-recovery data
## 9.5.3. Analysis of the model

library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
set.seed(123)

## Read data
## The data generation code is in bpa-code.txt, available at
## http://www.vogelwarte.ch/de/projekte/publikationen/bpa/complete-code-and-data-files-of-the-book.html
stan_data <- read_rdump("lifedead.data.R")

## Initial values
inits <- function() list(mean_s = runif(1, 0, 1),
                         mean_f = runif(1, 0, 1),
                         mean_p = runif(1, 0, 1),
                         mean_r = runif(1, 0, 1))

## Parameters monitored
params <- c("mean_s", "mean_f", "mean_r", "mean_p")

## MCMC settings
ni <- 2000
nt <- 1
nb <- 1000
nc <- 4

## Call Stan from R
lifedead <- stan("lifedead.stan",
                 data = stan_data, init = inits, pars = params,
                 chains = nc, iter = ni, warmup = nb, thin = nt,
                 seed = 1,
                 open_progress = FALSE)
print(lifedead, digits = 3)
