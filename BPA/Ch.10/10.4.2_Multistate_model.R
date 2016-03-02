## 10.4. Models with constant survival and time-dependent entry
## 10.4.2 Analysis of the JS model as a multistate model

library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
set.seed(123)

## Read data
## The data generation code is in bpa-code.txt, available at
## http://www.vogelwarte.ch/de/projekte/publikationen/bpa/complete-code-and-data-files-of-the-book.html
source("js-grouse.data.R")
CH.ms <- cbind(rep(2, dim(y)[1]), y)
CH.ms[CH.ms == 0] <- 2
stan_data = list(y = CH.ms,
                 n_occasions = dim(CH.ms)[2],
                 M = dim(CH.ms)[1])

## Initial values
inits <- function() list(mean_phi = runif(1, 0, 1),
                         mean_p = runif(1, 0, 1))

## Parameters monitored
params <- c("psi", "mean_p", "mean_phi", "b",
            "Nsuper", "N", "B", "gamma")

## MCMC settings
ni <- 2000
nt <- 1
nb <- 1000
nc <- 4

## Call Stan from R
js_ms <- stan("js_ms.stan",
               data = stan_data, init = inits, pars = params,
               chains = nc, iter = ni, warmup = nb, thin = nt,
               seed = 1,
               open_progress = FALSE)
print(js_ms, digits = 3)
