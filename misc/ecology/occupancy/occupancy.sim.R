R <- 250;
T <- 10;
p <- 0.55;
psi <- 0.4;

z <- rbinom(R,1,psi);

y <- matrix(NA,R,T);
for (r in 1:R)
  y[r,] <- rbinom(T,1,z[r] * p);

# to fit in Stan, use this:
# library('rstan')
# fit <- stan('occupancy.stan', data=c("R","T","y"));

