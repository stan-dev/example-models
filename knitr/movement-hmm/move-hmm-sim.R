rcat <- function(theta) sample(1:length(theta), size=1, prob=theta);

N <- 10000;
theta <- matrix(c(0.90, 0.10, 0.10, 0.90),
                2, 2, byrow=TRUE);
z <- rep(NA, N);
z[1] <- 1;
for (n in 2:N) z[n] <- rcat(theta[z[n - 1], ]);

library(CircStats);
# rvm:  von mises on (0, 2pi)
# for (-pi, pi), adjust mu + pi, turn - pi
mu_turn <- c(-0.2, 0) + pi;
kappa_turn <- c(0.5, 4);
turn <- rep(NA, N);
for (n in 1:N) turn[n] <- rvm(1, mu_turn[z[n]], kappa_turn[z[n]]);
turn <- turn - pi;

alpha_dist <- c(1, 3);
sigma_dist <- c(0.25, 0.5);
dist <- rep(NA, N);
for (n in 1:N) dist[n] <- rlnorm(1, alpha_dist[z[n]], sigma_dist[z[n]]);
# was: rweibull(...)

df <- data.frame(turn, dist);
library(ggplot2);
plot <-
  ggplot(df, aes(x = turn, y = dist)) +
  geom_point(alpha = 0.1) +  
#  coord_polar() +  # values OK, labels are a mess with this
  scale_y_log10() +
  scale_x_continuous(breaks=c(-pi/2, 0, pi/2),
                     labels=c("-pi/2", "0", "pi/2"));
plot;


library(rstan);
K <- 2;
m <- stan_model("move-hmm.stan");
fit <- sampling(m, data=c("N", "turn", "dist", "K"),
                chains = 1, iter = 200,
                control=list(stepsize=0.05, adapt_delta=0.9), refresh=2);
