## ---- sim-1pl ----
inv_logit <- function(u) {
  return(1 / (1 + exp(-u)));
}

I <- 20;
J <- 100;
theta <- rnorm(J, 0, 1);
b <- rnorm(I, -1, 2);

y <- matrix(NA, I, J);
for (i in 1:I)
  y[i,] <- rbinom(J, 1, inv_logit(theta - b[i]));

## ---- sim-1pl-hist ----
library(ggplot2);
hist_theta_sim <-  
  ggplot(data=data.frame(theta), aes(theta)) + 
  geom_histogram(binwidth=0.5, colour="black", fill="white");
hist_theta_sim;

hist_b_sim <-
  ggplot(data=data.frame(b), aes(b)) +
  geom_histogram(binwidth=1, colour="black", fill="white");
hist_b_sim;
