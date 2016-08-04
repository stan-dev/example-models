theta <- 0.3;

K <- 2;
mu <- c(1, 3);
sigma <- c(0.5, 1);
N <- 1000;
y <- c(rlnorm(theta * N, mu[1], sigma[1]),
       rlnorm((1 - theta) * N, mu[2], sigma[2]));

fit <- stan("mix-normal.stan",
            data = c("y", "N"),
            chains = 1,
            control = list(stepsize=0.01, adapt_delta=0.99), fit = fit);
