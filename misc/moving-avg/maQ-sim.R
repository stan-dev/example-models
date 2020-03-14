mu <- -1.25;
sigma <- 0.75;
theta <- c(0.7,0.3,0.1);
T <- 1000;
Q <- 3;
y <- rep(0,T);
epsilon <- rep(0,T);

predictor <- mu;
y[1] <- rnorm(1,predictor,sigma);
epsilon[1] <- y[1] - predictor;

predictor <- mu + theta[1] * epsilon[1];
y[2] <- rnorm(1,predictor,sigma);
epsilon[2] <- y[2] - predictor;

predictor <- mu + theta[1] * epsilon[1] + theta[2] * epsilon[2];
y[3] <- rnorm(1,predictor,sigma);
epsilon[3] <- y[3] - predictor;

for (t in 4:T) {
  predictor <- mu + theta[1] * epsilon[t - 1] + theta[2] * epsilon[t - 2] + theta[3] * epsilon[t - 3];
  y[t] <- rnorm(1, predictor, sigma);
  epsilon[t] <- y[t] - predictor;
}