library("rstan")
library("ggplot2")

set.seed(123)

x <- (-50:50)/25
N <- length(x)

comp_sim <- stan_model(file="gp-sim.stan")
fit_sim <- sampling(comp_sim, data = list(x = x, N = N), iter = 100, chains = 3, seed = 123)

fit_sim_ss <- extract(fit_sim)

# DUMP DATA: gp-sim.stan
dump(c("N","x"),file="gp-sim.data.R")

# DUMP DATA: gp-fit.stan
y <- fit_sim_ss$y[1,] # any sample will do
dump(c("N","x","y"),file="gp-fit.data.R")  # need sample


# DUMP DATA FOR gp-predict.stan
N1 <- N
x1 <- x
y1 <- fit_sim_ss$y[100,]
x2 <- (-80:80)/25
N2 <- length(x2)

dump(c("N1","x1","y1","N2","x2"),file="gp-predict.data.R")

# DUMP DATA: gp-logit-predict.stan
y1 <- fit_sim_ss$y[1,] # any sample will do
z1 <- rep(0,length(y1))
z1 <- rbinom(N1, 1,plogis(y1))
dump(c("N1","x1","z1","N2","x2"),file="gp-logit-predict.data.R")

y <- rpois(N, lambda = exp(1.5 + fit_sim_ss$y[1,] - mean(fit_sim_ss$y[1,])))
dump(c('x','N','y'),file = 'gp-fit-pois.data.R')

z <- rbinom(N, 1, plogis(1 + fit_sim_ss$y[1,] - mean(fit_sim_ss$y[1,])))
dump(c('x','N','z'),file = 'gp-fit-logit.data.R')

x1 <- x
x2 <- x
x_mat <- expand.grid(x1, x2)
x <- x_mat[sample(dim(x_mat)[1], N),]
D <- 2

gp_mult <- stan(file = 'gp-sim-multi.stan', data = list(x = x, D = D, N = N), chains = 3, iter = 200, seed = 123)
gp_mult_ss <- rstan::extract(gp_mult)
y <- gp_mult_ss$y[200,]
dump(c('x','N','D','y'), file = 'gp-fit-ARD.data.R')

x <- x1
D <- 3
gp_mult_out <- stan(file = 'gp-sim-multi-output.stan', data = list(x = x, D = D, N = N),
                    chains = 3, iter = 200, seed = 123, control = list(max_treedepth = 15))
gp_mult_out_ss <- rstan::extract(gp_mult_out)
y <- gp_mult_out_ss$y[200,,]
Omega <- gp_mult_out_ss$Omega[200,,]

dump(c('x','N','D','y','Omega'), file = 'gp-fit-multi-output.data.R')
