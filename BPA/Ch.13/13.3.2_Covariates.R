## 13. Metapopulation modeling of species distributions using
##     hierarchical logistic regression: Site-occupancy models
## 13.3. Generation and analysis of simulated data for
##       single-season occupancy
## 13.3.2. Site-occupancy models with covariates

library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
set.seed(123)

## Read data
## The data generation code is in bpa-code.txt, available at
## http://www.vogelwarte.ch/de/projekte/publikationen/bpa/complete-code-and-data-files-of-the-book.html
stan_data <- read_rdump("site-occ_cov.data.R")
## The true number of occupied sites is 73 in this data.

## Parameters monitored
params <- c("alpha_occ", "beta_occ", "alpha_p", "beta_p", "occ_fs")

## MCMC settings
ni <- 2000
nt <- 1
nb <- 1000
nc <- 4

## Initial values
inits <- lapply(1:nc, function(i)
    list(alpha_occ = runif(1, -3, 3),
         beta_occ = runif(1, -3, 3),
         alpha_p = runif(1, -3, 3),
         beta_p = runif(1, -3, 3)))

## Call Stan from R
out <- stan("site-occ_cov.stan",
            data = stan_data,
            init = inits, pars = params,
            chains = nc, iter = ni, warmup = nb, thin = nt,
            seed = 1,
            open_progress = FALSE)
print(out, digits = 2)
