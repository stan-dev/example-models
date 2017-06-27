## 6. Estimation of the size of a closed population
## 6.2. Generation and analysis of simulated data with data
## augmentation
## 6.3. Analysis of a real data set: model Mtbh for species richness estimation

library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
set.seed(1234)

## Read in data
## "p610.txt" is available at
## http://www.vogelwarte.ch/de/projekte/publikationen/bpa/complete-code-and-data-files-of-the-book.html
p610 <- read.table("p610.txt", header = TRUE)
y <- p610[, 5:9]                        # Grab counts
y[y > 1] <- 1                           # Counts to det-nondetections

## Bundle data
stan_data <- list(y = as.matrix(y), M = nrow(y), T = ncol(y))

## Parameters monitored
params <- c("N", "mean_p", "gamma", "sigma", "omega")

## MCMC settings
ni <- 10000
nt <- 9
nb <- 1000
nc <- 4

## Initial values
inits <- lapply(1:nc, function(i) {
    list(mean_p = runif(ncol(y), 0.1, 0.5),
         sigma = runif(1, 0.5, 1.0),
         omega = runif(1, 0.2, 0.5),
         gamma = runif(1, -0.5, 0.5),
         eps = runif(nrow(y), -1, 1))})

## Call Stan from R
out <- stan("Mtbh.stan",
            data = stan_data, init = inits, pars = params,
            chains = nc, iter = ni, warmup = nb, thin = nt,
            seed = 1,
            control = list(adapt_delta = 0.8),
            open_progress = FALSE)

## Summarize posteriors
print(out, digits = 3)

## Model M0

## Define parameters to be monitored
params <- c("N", "p", "omega")

## Call Stan from R
out0 <- stan("M0.stan",
             data = stan_data, pars = params,
             chains = nc, iter = ni, warmup = nb, thin = nt,
             seed = 1,
             open_progress = FALSE)

## Summarize posteriors
print(out0, digits = 3)
