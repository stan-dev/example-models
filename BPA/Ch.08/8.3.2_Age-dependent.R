## 8. Estimation of survival probabilities using mark-recovery data
## 8.3. The mark-recovery model fitted with the multinomial likelihood 
## 8.3.2. Age-dependent parameters

library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
set.seed(123)

## Read data
## The data generation code is in bpa-code.txt, available at
## http://www.vogelwarte.ch/de/projekte/publikationen/bpa/complete-code-and-data-files-of-the-book.html
stan_data <- read_rdump("mr_mnl_age.data.R")

## Initial values
inits <- function() list(mean_sj = runif(1, 0, 1),
                         mean_sa = runif(1, 0, 1),
                         mean_rj = runif(1, 0, 1),
                         mean_ra = runif(1, 0, 1))

## Parameters monitored
params <- c("mean_sj", "mean_rj", "mean_sa", "mean_ra")

## MCMC settings
ni <- 2000
nt <- 1
nb <- 1000
nc <- 4

## Call Stan from R
mr_age <- stan("mr_mnl_age.stan",
               data = stan_data, init = inits, pars = params,
               chains = nc, iter = ni, warmup = nb, thin = nt,
               seed = 1,
               open_progress = FALSE)

## Summarize posteriors
print(mr_age, digits = 3)
