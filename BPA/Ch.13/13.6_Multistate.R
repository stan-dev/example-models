## 13. Metapopulation modeling of species distributions using
##     hierarchical logistic regression: Site-occupancy models
## 13.6. Multistate occupancy models

library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
set.seed(123)

## Read data
## The data file "owls.txt" is available at
## http://www.vogelwarte.ch/de/projekte/publikationen/bpa/complete-code-and-data-files-of-the-book.html
owls <- read.table("owls.txt", header = TRUE)
y <- as.matrix(owls[, 2:6])
y <- y + 1
y[is.na(y)] <- 0
stan_data <- list(y = y, R = dim(y)[1], T = dim(y)[2])

## Parameters monitored
params <- c("p2", "p3", "r", "psi", "n_occ")

## MCMC settings
ni <- 4000
nt <- 3
nb <- 1000
nc <- 4

## Model 1
## Initial values
inits <- lapply(1:nc, function(i)
    list(p2 = runif(1, 0, 1),
         psi = runif(1, 0, 1),
         r = runif(1, 0, 1),
         beta = runif(3, 0, 2)))

## Call Stan from R
out1 <- stan("owls_ms1.stan",
             data = stan_data,
             init = inits, pars = params,
             chains = nc, iter = ni, warmup = nb, thin = nt,
             seed = 1,
             control = list(adapt_delta = 0.8),
             open_progress = FALSE)
print(out1, digits = 2)

## Model 2
## Initial values
inits <- lapply(1:nc, function(i)
    list(p2 = runif(stan_data$T, 0, 1),
         psi = runif(1, 0, 1),
         r = runif(1, 0, 1)))

out2 <- stan("owls_ms2.stan",
             data = stan_data,
             init = inits, pars = params,
             chains = nc, iter = ni, warmup = nb, thin = nt,
             seed = 1,
             control = list(adapt_delta = 0.8),
             open_progress = FALSE)
print(out2, digits = 2)
