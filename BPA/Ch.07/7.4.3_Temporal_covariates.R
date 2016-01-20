## 7. Estimation of survival probabilities using capture-recapture data
## 7.4. Models with time-variation
## 7.4.3. Temporal covariates

library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
set.seed(123)

## Read data
## The data generation code is in bpa-code.txt, available at
## http://www.vogelwarte.ch/de/projekte/publikationen/bpa/complete-code-and-data-files-of-the-book.html
stan_data <- read_rdump("cjs_cov_raneff.data.R")

## Initial values
inits <- function(){list(mu = rnorm(1),
                         sigma = runif(1, 0, 5),
                         beta = runif(1, -5, 5),
                         mean_p = runif(1, 0, 1))}  

## Parameters monitored
params <- c("mean_phi", "mean_p", "phi_est", "sigma2", "beta")

# MCMC settings
ni <- 4000
nt <- 3
nb <- 1000
nc <- 4

## Call Stan from R
cjs_cov <- stan("cjs_cov_raneff.stan",
                data = stan_data, init = inits, pars = params,
                chains = nc, iter = ni, warmup = nb, thin = nt,
                seed = 1,
                open_progress = FALSE)

## Summarize posteriors
print(cjs_cov, digits = 3)
