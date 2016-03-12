## 13. Metapopulation modeling of species distributions using
##     hierarchical logistic regression: Site-occupancy models
## 13.3. Generation and analysis of simulated data for
##       single-season occupancy
## 13.3.1. The simplest possible site-occupancy model

library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
set.seed(123)

## Read data
## The data generation code is in bpa-code.txt, available at
## http://www.vogelwarte.ch/de/projekte/publikationen/bpa/complete-code-and-data-files-of-the-book.html
stan_data <- read_rdump("site-occ.data.R")
## The true number of occupied sites is 160 in this data.

## Parameters monitored
params <- c("psi", "p", "occ_fs")

## MCMC settings
ni <- 2000
nt <- 1
nb <- 1000
nc <- 4

## Initial values
inits <- lapply(1:nc, function(i)
    list(psi = runif(1, 0, 1),
         p = runif(1, 0, 1)))

## Call Stan from R
out <- stan("site-occ.stan",
            data = stan_data,
            init = inits, pars = params,
            chains = nc, iter = ni, warmup = nb, thin = nt,
            seed = 1,
            control = list(adapt_delta = 0.8),
            open_progress = FALSE)
print(out, digits = 2)
