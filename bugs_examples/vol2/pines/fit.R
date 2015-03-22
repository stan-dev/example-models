library(rstan);
source('pines.data.R');
# fit_independent <- stan('pines-independent.stan', data=c("N", "y", "x", "z"));
fit <- stan('pines.stan', data=c("N", "y", "x", "z"), iter=100000)
