## 9. Multistate capture-recapture models
## 9.6. Estimation of movement among three sites
## 9.6.3. Analysis of the model

library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
set.seed(123)

## Read data
## The data generation code is in bpa-code.txt, available at
## http://www.vogelwarte.ch/de/projekte/publikationen/bpa/complete-code-and-data-files-of-the-book.html
stan_data <- read_rdump("ms3_multinomlogit.data.R")

## Initial values
inits <- function() list(phiA = runif(1, 0, 1),
                         phiB = runif(1, 0, 1),
                         phiC = runif(1, 0, 1),
                         lpsiA = rnorm(2),
                         lpsiB = rnorm(2),
                         lpsiC = rnorm(2),
                         pA = runif(1, 0, 1),
                         pB = runif(1, 0, 1),
                         pC = runif(1, 0, 1))

## Parameters monitored
params <- c("phiA", "phiB", "phiC",
            "psiA", "psiB", "psiC",
            "pA", "pB", "pC")

## MCMC settings
ni <- 2000
nt <- 1
nb <- 1000
nc <- 4

## Call Stan from R
ms3 <- stan("ms3_multinomlogit.stan",
            data = stan_data, init = inits, pars = params,
            chains = nc, iter = ni, warmup = nb, thin = nt,
            seed = 1,
            open_progress = FALSE)
print(ms3, digits = 3)
