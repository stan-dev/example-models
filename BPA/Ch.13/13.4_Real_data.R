## 13. Metapopulation modeling of species distributions using
##     hierarchical logistic regression: Site-occupancy models
## 13.4. Analysis of real data set: Single-season occupancy model

library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
set.seed(123)

## Read data
## The data file "bluebug.txt" is available at
## http://www.vogelwarte.ch/de/projekte/publikationen/bpa/complete-code-and-data-files-of-the-book.html
data <- read.table("bluebug.txt", header = TRUE)

# Collect the data into suitable structures
y <- as.matrix(data[, 4:9])
y[y > 1] <- 1
edge <- data$forest_edge
dates <- as.matrix(data[, 10:15])
hours <- as.matrix(data[, 16:21])

# Standardize covariates
mean.date <- mean(dates, na.rm = TRUE)
sd.date <- sd(dates[!is.na(dates)])
DATES <- (dates-mean.date) / sd.date
DATES[is.na(DATES)] <- 0

mean.hour <- mean(hours, na.rm = TRUE)
sd.hour <- sd(hours[!is.na(hours)])
HOURS <- (hours-mean.hour) / sd.hour
HOURS[is.na(HOURS)] <- 0

last <- sapply(1:dim(y)[1],
               function(i) max(grep(FALSE, is.na(y[i, ]))))
y[is.na(y)] <- 0

stan_data <- list(y = y, R = nrow(y), T = ncol(y), edge = edge,
                  DATES = DATES, HOURS = HOURS, last = last)

## Parameters monitored
params <- c("alpha_psi", "beta_psi", "mean_p", "occ_fs",
            "alpha_p", "beta1_p", "beta2_p", "beta3_p", "beta4_p")

## MCMC settings
ni <- 6000
nt <- 5
nb <- 1000
nc <- 4

## Initial values
inits <- lapply(1:nc, function(i)
    list(alpha_psi = runif(1, -3, 3),
         alpha_p = runif(1, -3, 3)))

## Call Stan from R
out <- stan("bluebug.stan",
            data = stan_data,
            init = inits, pars = params,
            chains = nc, iter = ni, warmup = nb, thin = nt,
            seed = 1,
            control = list(adapt_delta = 0.8),
            open_progress = FALSE)
print(out, digits = 2)

## Posteriors of alpha_psi and beta_psi will be somewhat different
## from those in the book. This may be because convergences of these
## parameters in WinBUGS are not good, as described in the text.
## JAGS will produce more similar results to Stan.

hist(extract(out, pars = "occ_fs")$occ_fs, nclass = 30, col = "gray")
