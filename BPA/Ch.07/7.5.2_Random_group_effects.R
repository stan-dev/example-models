## 7. Estimation of survival probabilities using capture-recapture data
## 7.5. Models with individual variation
## 7.5.2. Random group effects

library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
set.seed(123)

## Read data
## The data generation code is in bpa-code.txt, available at
## http://www.vogelwarte.ch/de/projekte/publikationen/bpa/complete-code-and-data-files-of-the-book.html
stan_data <- read_rdump("cjs_group.data.R")

## Parameters monitored
params <- c("mean_phi", "mean_p", "phi_g", "sigma")

## MCMC settings
ni <- 4000
nt <- 1
nb <- 3000
nc <- 4

## Initial values
inits <- lapply(1:nc, function(i) {
    list(phi_g = runif(length(unique(stan_data$group)), 0, 1),
         p_g = runif(length(unique(stan_data$group)), 0, 1))})

## Call Stan from R
cjs_group_raneff <- stan("cjs_group_raneff.stan",
                         data = stan_data, init = inits, pars = params,
                         chains = nc, iter = ni, warmup = nb, thin = nt,
                         control = list(adapt_delta = 0.99),
                         seed = 1,
                         open_progress = FALSE)
## Note: there may be divergent transitions after warmup.

## Summarize posteriors
print(cjs_group_raneff, digits = 3)
