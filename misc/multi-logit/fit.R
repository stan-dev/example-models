library(rstan);
fit <- stan("multi_logit.stan", data=c("K","D","N","x","y"),
            chains=4, iter=500);
