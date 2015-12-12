## 6. Estimation of the size of a closed population
## 6.2. Generation and analysis of simulated data with data
## augmentation
## 6.2.4. Individual (random) effects: the heterogeneity model Mh

library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
set.seed(1234)

## Generate simulated data
## data.fn() is defined in bpa-code.txt, available at
## http://www.vogelwarte.ch/de/projekte/publikationen/bpa/complete-code-and-data-files-of-the-book.html
data <- data.fn()

## Aggregate capture-histories and augment data set
y <- sort(apply(data$yobs, 1, sum), decreasing = TRUE)
nz <- 300
yaug <- c(y, rep(0, nz))

## Bundle data
stan_data <- list(y = yaug, M = length(yaug), T = ncol(data$yobs))

## Initial values
inits <- function() list(mean_p = 0.5, sigma = runif(1, 0.5, 0.9))

## Parameters monitored
params <- c("N", "mean_p", "sigma", "omega")

## MCMC settings
ni <- 15000
nt <- 5
nb <- 10000
nc <- 4

## Call Stan from R
out <- stan("Mh.stan",
            data = stan_data, init = inits, pars = params,
            chains = nc, iter = ni, warmup = nb, thin = nt,
            seed = 1,
            open_progress = FALSE)
## Note: There may be divergent transitions after warmup.

## Summarize posteriors
print(out, digits = 3)
