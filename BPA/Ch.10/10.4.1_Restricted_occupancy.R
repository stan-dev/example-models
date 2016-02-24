## 10.4. Models with constant survival and time-dependent entry
## 10.4.1 Analysis of the JS model as a restricted occupancy model

library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
set.seed(123)

## Read data
## The data generation code is in bpa-code.txt, available at
## http://www.vogelwarte.ch/de/projekte/publikationen/bpa/complete-code-and-data-files-of-the-book.html
stan_data <- read_rdump("js-grouse.data.R")

## Parameters monitored
params <- c("psi", "mean_p", "mean_phi", "b",
            "Nsuper", "N", "B", "gamma")

## MCMC settings
ni <- 2000
nt <- 1
nb <- 1000
nc <- 4

## Initial values
inits <- lapply(1:nc, function(i)
    list(mean_phi = runif(1, 0, 1),
         mean_p = runif(1, 0, 1),
         gamma = runif(stan_data$n_occasions, 0, 1)))

## Call Stan from R
js_occ <- stan("js_rest_occ.stan",
               data = stan_data, init = inits, pars = params,
               chains = nc, iter = ni, warmup = nb, thin = nt,
               seed = 1,
               open_progress = FALSE)
print(js_occ, digits = 3)
