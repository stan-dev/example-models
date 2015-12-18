## 4.3.2. Analysis of real data set

library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
set.seed(123)

## Read data
## "tits.txt" is available at
## http://www.vogelwarte.ch/de/projekte/publikationen/bpa/complete-code-and-data-files-of-the-book.html
## The data conversion is derived from original bpa-code.txt

tits <- read.table("tits.txt", header = TRUE)
C <- as.matrix(tits[5:13])

obs <- as.matrix(tits[14:22])
a <- as.numeric(levels(factor(obs)))     # All the levels, numeric
newobs <- obs                            # Gets ObsID from 1:271
for (j in 1:length(a)) {
  newobs[which(obs == a[j])] <- j
}
newobs[is.na(newobs)] <- 272

first <- as.matrix(tits[23:31])
first[is.na(first)] <- 0

## Separate missing data
df <- expand.grid(site = 1:nrow(C), year = 1:ncol(C))
df$count <- c(C)
obsdata <- subset(df, !is.na(count))
misdata <- subset(df, is.na(count))

## Bundle data
stan_data <- list(nobs = nrow(obsdata),
                  nmis = nrow(misdata),
                  nyear = ncol(C),
                  nsite = nrow(C),
                  obs = obsdata$count,
                  obsyear = obsdata$year,
                  obssite = obsdata$site,
                  misyear = misdata$year,
                  missite = misdata$site)

##
##  (a) Null or intercept-only model
##

## Initial values
inits <- function() list(alpha = runif(1, -10, 10))

## Parameters monitored
params <- c("alpha", "mis")

## MCMC settings
ni <- 2000
nt <- 1
nb <- 1000
nc <- 4

## Call Stan from R
out0 <- stan("GLM0.stan", data = stan_data,
             init = inits, pars = params,
             chains = nc, iter = ni, warmup = nb, thin = nt,
             seed = 1,
             open_progress = FALSE)

## Summarize posteriors
print(out0, dig = 2)

##
##  (b) Fixed site effects
##

## Bundle data
stan_data <- list(nobs = nrow(obsdata),
                  nmis = nrow(misdata),
                  nyear = ncol(C),
                  nsite = nrow(C),
                  obs = obsdata$count,
                  obsyear = obsdata$year,
                  obssite = obsdata$site,
                  misyear = misdata$year,
                  missite = misdata$site)

## Initial values
inits <- function() list(alpha = runif(235, -1, 1))

## Parameters monitored
params <- c("alpha", "mis")

## MCMC settings
ni <- 2000
nt <- 1
nb <- 1000
nc <- 4

## Call Stan from R
out1 <- stan("GLM1.stan", data = stan_data,
             init = inits, pars = params,
             chains = nc, iter = ni, warmup = nb, thin = nt,
             seed = 1,
             open_progress = FALSE)

## Summarize posteriors
print(out1, digits = 2)

##
##  (c) Fixed site and fixed year effects
##

## Bundle data
stan_data <- list(nobs = nrow(obsdata),
                  nmis = nrow(misdata),
                  nyear = ncol(C),
                  nsite = nrow(C),
                  obs = obsdata$count,
                  obsyear = obsdata$year,
                  obssite = obsdata$site,
                  misyear = misdata$year,
                  missite = misdata$site)

## Initial values
inits <- function() list(alpha = runif(235, -1, 1),
                         eps = c(runif(8, -1, 1)))


## Parameters monitored
params <- c("alpha", "eps", "mis")

## MCMC settings
ni <- 2000
nt <- 1
nb <- 1000
nc <- 4

## Call Stan from R
out2 <- stan("GlM2.stan", data = stan_data,
             init = inits, pars = params,
             chains = nc, iter = ni, warmup = nb, thin = nt,
             seed = 1,
             open_progress = FALSE)

## Summarize posteriors
print(out2, digits = 2)

##
##  (d) Random site effects (no year effects)
##

## Bundle data
stan_data <- list(nobs = nrow(obsdata),
                  nmis = nrow(misdata),
                  nyear = ncol(C),
                  nsite = nrow(C),
                  obs = obsdata$count,
                  obsyear = obsdata$year,
                  obssite = obsdata$site,
                  misyear = misdata$year,
                  missite = misdata$site)

## Initial values
inits <- function() list(mu.alpha = runif(1, 2, 3))

## Parameters monitored
params <- c("alpha", "mu_alpha", "sd_alpha", "mis")

## MCMC settings
ni <- 2000
nt <- 1
nb <- 1000
nc <- 4

## Call Stan from R
out3 <- stan("GLMM1.stan", data = stan_data,
             init = inits, pars = params,
             chains = nc, iter = ni, warmup = nb, thin = nt,
             seed = 1,
             open_progress = FALSE)

## Summarize posteriors
print(out3, digits = 3)

##
##  (e) Random site and random year effects
##

## Bundle data
stan_data <- list(nobs = nrow(obsdata),
                  nmis = nrow(misdata),
                  nyear = ncol(C),
                  nsite = nrow(C),
                  obs = obsdata$count,
                  obsyear = obsdata$year,
                  obssite = obsdata$site,
                  misyear = misdata$year,
                  missite = misdata$site)

## Initial values (not required for all)
inits <- function() list(mu = runif(1, 0, 4),
                         alpha = runif(235, -2, 2),
                         eps = runif(9, -1, 1))

## Parameters monitored
params <- c("mu", "alpha", "eps", "sd_alpha", "sd_eps", "mis")

## MCMC settings
ni <- 6000
nt <- 5
nb <- 1000
nc <- 4

out4 <- stan("GLMM2.stan", data = stan_data,
             init = inits, pars = params,
             chains = nc, iter = ni, warmup = nb, thin = nt,
             seed = 1,
             open_progress = FALSE)

## Summarize posteriors
print(out4, digits = 2)

##
## (f) Random site and random year effects and first-year fixed
## observer effect
##

## Bundle data
stan_data <- list(nobs = nrow(obsdata),
                  nmis = nrow(misdata),
                  nyear = ncol(C),
                  nsite = nrow(C),
                  obs = obsdata$count,
                  obsyear = obsdata$year,
                  obssite = obsdata$site,
                  misyear = misdata$year,
                  missite = misdata$site,
                  first = t(first))

## Initial values
inits <- function() list(mu = runif(1, 0, 4),
                         beta2 = runif(1, -1, 1),
                         alpha = runif(235, -2, 2),
                         eps = runif(9, -1, 1))

## Parameters monitored
params <- c("mu", "beta2", "alpha", "eps", "sd_alpha", "sd_eps",
            "mis")

## MCMC settings
ni <- 6000
nt <- 5
nb <- 1000
nc <- 4

## Call Stan from R
out5 <- stan("GLMM3.stan", data = stan_data,
             init = inits, pars = params,
             chains = nc, thin = nt, iter = ni, warmup = nb,
             seed = 1,
             open_progress = FALSE)

## Summarize posteriors
print(out5, digits = 2)

##
##  (g) Random site and random year effects, first-year fixed observer
##  effect and overall linear time trend
##

## Bundle data
stan_data <- list(nobs = nrow(obsdata),
                  nmis = nrow(misdata),
                  nyear = ncol(C),
                  nsite = nrow(C),
                  obs = obsdata$count,
                  obsyear = obsdata$year,
                  obssite = obsdata$site,
                  misyear = misdata$year,
                  missite = misdata$site,
                  first = t(first),
                  year = ((1:9) - 5) / 4)

## Initial values
inits <- function() list(mu = runif(1, 0, 4),
                         beta1 = runif(1, -1, 1),
                         beta2 = runif(1, -1, 1),
                         alpha = runif(235, -2, 2),
                         eps = runif(9, -1, 1))

## Parameters monitored
params <- c("mu", "beta1", "beta2", "alpha", "eps",
            "sd_alpha", "sd_eps", "mis")

# MCMC settings
ni <- 9000
nt <- 8
nb <- 1000
nc <- 4

## Call Stan from R
out6 <- stan("GLMM4.stan", data = stan_data,
             init = inits, pars = params,
             chains = nc, iter = ni, warmup = nb, thin = nt,
             seed = 1,
             open_progress = FALSE)

## Summarize posteriors
print(out6, digits = 2)


##
## (h) The full model 
##

## Bundle data
stan_data <- list(nobs = nrow(obsdata),
                  nmis = nrow(misdata),
                  nyear = ncol(C),
                  nsite = nrow(C),
                  obs = obsdata$count,
                  obsyear = obsdata$year,
                  obssite = obsdata$site,
                  misyear = misdata$year,
                  missite = misdata$site,
                  nnewobs = 272,
                  newobs = t(newobs),
                  first = t(first),
                  year = ((1:9)-5) / 4)

## Initial values
inits <- function() list(mu = runif(1, 0, 4),
                         beta1 = runif(1, -1, 1),
                         beta2 = runif(1, -1, 1),
                         alpha = runif(235, -1, 1),
                         gamma = runif(272, -1, 1),
                         eps = runif(9, -1, 1))

## Parameters monitored
params <- c("mu", "beta1", "beta2", "alpha", "gamma", "eps",
            "sd_alpha", "sd_gamma", "sd_eps", "mis")

## MCMC settings
ni <- 6000
nt <- 4
nb <- 2000
nc <- 4

## Call Stan from R
out7 <- stan("GLMM5.stan", data = stan_data,
             init = inits, pars = params,
             chains = nc, iter = ni, warmup = nb, thin = nt,
             seed = 1,
             control = list(adapt_delta = 0.9),
             open_progress = FALSE)

## Summarize posteriors
print(out7, digit = 2)
