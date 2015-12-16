## 7. Estimation of survival probabilities using capture-recapture data
## 7.6. Models with time and group effects
## 7.6.1. Fixed group and time effects

library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
set.seed(123)

## Read data
## The data generation code is in bpa-code.txt, available at
## http://www.vogelwarte.ch/de/projekte/publikationen/bpa/complete-code-and-data-files-of-the-book.html
stan_data <- read_rdump("cjs_add.data.R")

## Parameters monitored
params <- c("phi_g1", "phi_g2", "p_g", "beta")

## MCMC settings
ni <- 3000
nt <- 2
nb <- 1000
nc <- 4

## Initial values
inits <- lapply(1:nc, function(i) {
    list(gamma = rnorm(stan_data$n_occasions - 1),
         beta = c(0, rnorm(1)),
         p_g = runif(length(unique(stan_data$group)), 0, 1))})

## Call Stan from R
cjs_add  <- stan("cjs_add.stan",
                 data = stan_data, init = inits, pars = params,
                 chains = nc, iter = ni, warmup = nb, thin = nt,
                 seed = 1,
                 open_progress = FALSE)

## Summarize posteriors
print(cjs_add, digits = 3)
