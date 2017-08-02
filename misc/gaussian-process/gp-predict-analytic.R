library(rstan)

stan_dat <- read_rdump('gp-predict.data.R')

fit_predict_analytic <- stan(file="gp-predict-analytic.stan",   
                             data=stan_dat,
                             iter=200, chains=3)

y2 <- rstan::extract(fit_predict_analytic, pars = c('y2'))$y2
print(fit_predict_analytic, pars = c('rho','alpha','sigma'))

# plot fits vs. simulated value from which fits drawn
df <- data.frame(x2=stan_dat$x2, y2_samp=colMeans(y2))
plot <- qplot(x2,y2_samp, data=df)


