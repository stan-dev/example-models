## 7. Estimation of survival probabilities using capture-recapture data
## 7.5. Models with individual variation
## 7.5.1. Fixed group effects

library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
set.seed(123)

## Read data
## The data generation code is in bpa-code.txt, available at
## http://www.vogelwarte.ch/de/projekte/publikationen/bpa/complete-code-and-data-files-of-the-book.html
stan_data <- read_rdump("cjs_group.data.R")

## Parameters monitored
params <- c("phi_g", "p_g")

## MCMC settings
ni <- 4000
nt <- 3
nb <- 1000
nc <- 4

## Initial values
inits <- lapply(1:nc, function(i) {
    list(phi_g = runif(length(unique(stan_data$group)), 0, 1),
         p_g = runif(length(unique(stan_data$group)), 0, 1))})

## Call Stan from R
cjs_group <- stan("cjs_group.stan",
                  data = stan_data, init = inits, pars = params,
                  chains = nc, iter = ni, warmup = nb, thin = nt,
                  seed = 1,
                  open_progress = FALSE)

## Summarize posteriors
print(cjs_group, digits = 3)
