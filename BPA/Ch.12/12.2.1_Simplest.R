## 12. Metapopulation modeling of abundance using
##     hierarchical Poisson regression: binomial mixture models
## 12.2. Generation and analysis of simulated data
## 12.2.1. The simplest case with constant parameters

library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
set.seed(123)

## Read data
## The data generation code is in bpa-code.txt, available at
## http://www.vogelwarte.ch/de/projekte/publikationen/bpa/complete-code-and-data-files-of-the-book.html
source("binmix.data.R")

## Parameters monitored
params <- c("lambda", "p")

## MCMC settings
ni <- 1000
nt <- 1
nb <- 500
nc <- 4

## Initial values
inits <- lapply(1:nc, function(i)
    list(p = runif(1, 0, 1),
         lambda = runif(1, 0, 1)))

## Call Stan from R
out <- stan("binmix.stan",
            data = list(y = y, R = R, T = T, K = 100),
            init = inits, pars = params,
            chains = nc, iter = ni, warmup = nb, thin = nt,
            seed = 1,
            open_progress = FALSE)
print(out, digits = 3)
