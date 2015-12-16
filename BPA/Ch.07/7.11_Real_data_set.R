## 7. Estimation of survival probabilities using capture-recapture data
## 7.11. Analysis of a real data set: survival of female Leisler’s bats

library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
set.seed(123)

## Read data
## The data generation code is in bpa-code.txt, available at
## http://www.vogelwarte.ch/de/projekte/publikationen/bpa/complete-code-and-data-files-of-the-book.html
stan_data <- read_rdump("cjs_mnl_ran.data.R")

## Parameters monitored
params <- c("phi", "mean_p", "mean_phi", "sigma2", "sigma2_real",
            "fit", "fit_new")

## MCMC settings
ni <- 3000
nt <- 2
nb <- 1000
nc <- 4

## Initial values
inits <- function() list(mean_phi = runif(1, 0, 1),
                         sigma = runif(1, 0, 5),
                         mean_p = runif(1, 0, 1))

## Call Stan from R
leis_result <- stan("cjs_mnl_ran.stan",
                    data = stan_data, init = inits, pars = params,
                    chains = nc, iter = ni, warmup = nb, thin = nt,
                    seed = 1,
                    open_progress = FALSE)

## Summarize posteriors
print(leis_result, digits = 3)
