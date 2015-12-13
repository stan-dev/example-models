## 7. Estimation of survival probabilities using capture-recapture data
## 7.6. Models with time and group effects
## 7.6.2. Fixed group and random time effects

library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
set.seed(123)

## Read data
## The data generation code is in bpa-code.txt, available at
## http://www.vogelwarte.ch/de/projekte/publikationen/bpa/complete-code-and-data-files-of-the-book.html
source("cjs_add.data.R")

## Bundle data
stan_data <- list(y = y, nind = nind, n_occasions = n_occasions,
                  g = g, group = group, R = R, df = df)

## Parameters monitored
params <- c("eta_phi", "p_g", "Sigma", "mean_phi")

## MCMC settings
ni <- 2000
nt <- 1
nb <- 1000
nc <- 4

## Initial values
inits <- function() list(p_g = runif(length(unique(group)), 0, 1),
                         Omega = matrix(c(1, 0, 0, 1), ncol = 2))

## Call Stan from R
cjs_temp_corr <- stan("cjs_temp_corr.stan",
                      data = stan_data, init = inits, pars = params,
                      chains = nc, iter = ni, warmup = nb, thin = nt,
                      seed = 1,
                      open_progress = FALSE)

## Summarize posteriors
print(cjs_temp_corr, digits = 3)
