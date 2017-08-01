library("rstan");
library("ggplot2");

x <- (-50:50)/25;
N <- length(x);

comp_sim <- stan_model(file="gp-sim.stan");
fit_sim <- sampling(comp_sim, data = list(x = x, N = N), iter = 300, chains = 3, seed = 123)

fit_sim_ss <- extract(fit_sim);

df <- data.frame(x=x,y_sim=colMeans(fit_sim_ss$y));
plot <- qplot(x,y_sim, data=df, xlim=c(-2,2),geom = 'line');

# DUMP DATA: gp-sim.stan
dump(c("N","x"),file="gp-sim.data.R");

# DUMP DATA: gp-fit.stan
y <- fit_sim_ss$y[1,]; # any sample will do
dump(c("N","x","y"),file="gp-fit.data.R");  # need sample

# DUMP DATA FOR gp-predict.stan
N1 <- N;
x1 <- x;
y1 <- fit_sim_ss$y[100,];
x2 <- (-80:80)/25;
N2 <- length(x2);

dump(c("N1","x1","y1","N2","x2"),file="gp-predict.data.R");

# DUMP DATA: gp-logit-predict.stan
y1 <- fit_sim_ss$y[1,]; # any sample will do
z1 <- rep(0,length(y1));
set.seed(123)
for (n in 1:N1)
  z1[n] <- rbinom(1,1,plogis(y1[n]));
dump(c("N1","x1","z1","N2","x2"),file="gp-logit-predict.data.R");



