## 6. Estimation of the size of a closed population
## 6.4. Capture-recapture models with individual covariates: model
## Mt+X
## 6.4.1. Individual covariate model for species richness estimation

library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
set.seed(1234)

## Generate simulated data
## data.fn() is defined in bpa-code.txt, available at
## http://www.vogelwarte.ch/de/projekte/publikationen/bpa/complete-code-and-data-files-of-the-book.html
p610 <- read.table("p610.txt", header = TRUE)
y <- p610[,5:9]                         # Grab counts
y[y > 1] <- 1                           # Convert to det-nondetections
ever.observed <- apply(y, 1, max)
wt <- p610$bm[ever.observed == 1]       # Body mass
yy <- as.matrix(y[ever.observed == 1,]) # Detection histories
dimnames(yy) <- NULL

## Augment both data sets
nz <- 150
yaug <- rbind(yy, array(0, dim = c(nz, ncol(yy))))
logwt3 <- c(log(wt^(1/3)), rep(NA, nz))

## Bundle data
bsize <- logwt3[1:nrow(yy)]
stan_data <- list(y = yaug,
                  bsize = bsize - mean(bsize),
                  M = nrow(yaug),
                  T = ncol(yaug),
                  C = nrow(yy),
                  prior_sd_upper = 3)

## Initial values
inits <- function() list(beta = runif(1, 0, 1),
                         mu_size = rnorm(1, 0, 1))

## Parameters monitored
params <- c("N", "mean_p", "beta", "omega", "mu_size", "sd_size")

## MCMC settings
ni <- 15000
nt <- 10
nb <- 5000
nc <- 4

## Call Stan from R
outX <- stan("MtX.stan",
             data = stan_data, init = inits, pars = params,
             chains = nc, iter = ni, warmup = nb, thin = nt,
             seed = 1,
             control = list(max_treedepth = 15),
             open_progress = FALSE)

## Summarize posteriors
print(outX, digits = 3)
