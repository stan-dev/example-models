## 9. Multistate capture-recapture models
## 9.4. Estimation of age-specific probability of first breeding
## 9.4.3. Analysis of the model

library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
set.seed(123)

## Read data
## The data generation code is in bpa-code.txt, available at
## http://www.vogelwarte.ch/de/projekte/publikationen/bpa/complete-code-and-data-files-of-the-book.html
stan_data <- read_rdump("agerecruitment.data.R")

## Initial values
inits <- function() list(mean_phi1 = runif(1, 0, 1),
                         mean_phi2 = runif(1, 0, 1),
                         mean_phiad = runif(1, 0, 1),
                         mean_alpha1 = runif(1, 0, 1),
                         mean_alpha2 = runif(1, 0, 1),
                         mean_pNB = runif(1, 0, 1),
                         mean_pB = runif(1, 0, 1))

## Parameters monitored
params <- c("mean_phi1", "mean_phi2", "mean_phiad",
            "mean_alpha1", "mean_alpha2",
            "mean_pNB", "mean_pB")

## MCMC settings
ni <- 2000
nt <- 1
nb <- 1000
nc <- 4

## Call Stan from R
agefirst <- stan("agerecruitment.stan",
                 data = stan_data, init = inits, pars = params,
                 chains = nc, iter = ni, warmup = nb, thin = nt,
                 seed = 1,
                 open_progress = FALSE)
print(agefirst, digits = 3)
