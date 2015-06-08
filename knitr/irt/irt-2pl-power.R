## ---- power-2pl-compile ----
library(rstan);
model <- stan_model("irt_2pl_power.stan");

sample_to_df <- function(test_name) {
  fit <- sampling(model, algorithm="Fixed_param",
                  data=c("I", "a", "b"), chains=1, iter=20000,
                  refresh=10000, seed=1234);
  sims <- extract(fit)$z_sim;
  questions <- c();
  mean <- c();
  sd <- c();
  theta_sim <- c();
  five <- c();
  fifty <- c();
  ninety_five <- c();
  for (j in 1:dim(sims)[2]) {
    questions[j] <- test_name;
    theta_sim[j] <- (j - 50) / 10;
    mean[j] <- mean(sims[,j]);
    sd[j] <- sd(sims[,j]);
    five[j] <- quantile(sims[,j], 0.05);
    fifty[j] <- quantile(sims[,j], 0.50);
    ninety_five[j] <- quantile(sims[,j], 0.95);
  }
  df <- data.frame(questions, mean, sd, theta_sim, five, fifty, ninety_five);
  return(df);
}

## ---- power-2pl-test-1 ----

## TEST 1a: low discrim
b <- ((0:20) - 10) / 2;
I <- length(b);
a <- rep(0.25, I);
df_1a <- sample_to_df("low discrim (a = 0.25)");

## TEST 1b: medium discrim
b <- ((0:20) - 10) / 2;
a <- rep(1, I);
df_1b <- sample_to_df("medium discrim (a = 1)");

## TEST 1c: high discrim
b <- ((0:20) - 10) / 2;
a <- rep(4, length(b));
df_1c <- sample_to_df("high discrim (a = 4)");

df_1abc <- rbind(df_1a,df_1b,df_1c);

## ---- power-2pl-plot-1

library(ggplot2);
plot_1abc <-
  ggplot(df_1abc, aes(x=theta_sim, y=mean)) +
  facet_grid(. ~ questions) +
  geom_ribbon(aes(ymin=(mean - 2 * sd), ymax=(mean + 2 * sd)),
              colour="darkgray", fill="lightgray") +
  geom_line(size=0.5) +
  xlab("ability") +
  ylab("expected correct +/- 2 std dev") +
  ggtitle("21 Questions of Mixed Difficulty (b in [-5, 5])");

plot(plot_1abc);

## ---- power-2pl-test-2 ----

## TEST 2a: 20 all easy
b <- rep(-3, I);
a <- rep(1, I);
df_2a <- sample_to_df("easy (b = -3)");

## TEST 2b: 20 all medium
b <- rep(0, I);
a <- rep(1, I);
df_2b <- sample_to_df("moderate (b = 0)");

## TEST 2c: 20 all hard
b <- rep(3, I);
a <- rep(1, I);
df_2c <- sample_to_df("difficult (b = 3)");

df_2abc <- rbind(df_2a, df_2b, df_2c);

plot_2abc <-
  ggplot(df_2abc, aes(x=theta_sim, y=mean)) +
  facet_grid(. ~ questions) +
  geom_ribbon(aes(ymin=(mean - 2 * sd), ymax=(mean + 2 * sd)),
              colour="darkgray", fill="lightgray") +
  geom_line(size=0.5) +
  xlab("ability") +
  ylab("expected correct +/- 2 std dev") +
  ggtitle("21 Questions of Same Difficulty (discrim a = 1)")
plot(plot_2abc);

