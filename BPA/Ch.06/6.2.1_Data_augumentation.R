## 6. Estimation of the size of a closed population
## 6.2. Generation and analysis of simulated data with data
## augmentation
## 6.2.1. Introduction to data augmentation for the simplest case:
## model M0

library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
set.seed(1234)

## Generate simulated data
## data.fn() is defined in bpa-code.txt, available at
## http://www.vogelwarte.ch/de/projekte/publikationen/bpa/complete-code-and-data-files-of-the-book.html
data <- data.fn()

## Augment data set by 150 potential individuals
nz <- 150
yaug <- rbind(data$yobs, array(0, dim = c(nz, data$T)))

## Bundle data
stan_data <- list(y = yaug, M = nrow(yaug), T = ncol(yaug))

## Initial values
inits <- function() {
    list(p = runif(1, 0, 1), omega = 0.5)}

## Parameters monitored
params <- c("N", "p", "omega")

## MCMC settings
ni <- 2000
nt <- 1
nb <- 1000
nc <- 4

## Call Stan from R
out <- stan("M0.stan",
            data = stan_data, init = inits, pars = params,
            chains = nc, iter = ni, warmup = nb, thin = nt,
            seed = 2,
            open_progress = FALSE)

## Summarize posteriors
print(out, digits = 3)
