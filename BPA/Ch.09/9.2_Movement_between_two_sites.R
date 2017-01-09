## 9. Multistate capture-recapture models
## 9.2. Estimation of movement between two sites
## 9.2.3. Analysis of the model

library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
set.seed(123)

## Read data
## The data generation code is in bpa-code.txt, available at
## http://www.vogelwarte.ch/de/projekte/publikationen/bpa/complete-code-and-data-files-of-the-book.html
stan_data <- read_rdump("ms.data.R")

## Initial values
inits <- function() list(mean_phi = runif(2, 0, 1),
                         mean_psi = runif(2, 0, 1),
                         mean_p = runif(2, 0, 1))

## Parameters monitored
params <- c("mean_phi", "mean_psi", "mean_p")

## MCMC settings
ni <- 2000
nt <- 1
nb <- 1000
nc <- 4

## Call Stan from R
ms <- stan("ms.stan",
           data = stan_data, init = inits, pars = params,
           chains = nc, iter = ni, warmup = nb, thin = nt,
           seed = 1,
           open_progress = FALSE)
print(ms, digits = 3)

## Alternative model 1
ms1 <- stan("ms_alternative1.stan",
            data = stan_data, init = inits, pars = params,
            chains = nc, iter = ni, warmup = nb, thin = nt,
            seed = 1,
            open_progress = FALSE)
print(ms1, digits = 3)

## Alternative model 2
params <- c("phiA", "phiB", "psiAB", "psiBA", "pA", "pB")
inits <- function() list(phiA = runif(1, 0, 1),
                         phiBA = runif(1, 0, 1),
                         psiAB = runif(1, 0, 1),
                         psiBA = runif(1, 0, 1),
                         pA = runif(1, 0, 1),
                         pB = runif(1, 0, 1))
ms2 <- stan("ms_alternative2.stan",
            data = stan_data, init = inits, pars = params,
            chains = nc, iter = ni, warmup = nb, thin = nt,
            seed = 1,
            open_progress = FALSE)
print(ms2, digits = 3)
