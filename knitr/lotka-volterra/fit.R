library(ggplot2)
library(rstan)
library(reshape)

lynx_hare_df <- read.csv("hudson-bay-lynx-hare.csv", comment.char="#")
lynx_hare_df$Lynx <- lynx_hare_df$Lynx # * 1000
lynx_hare_df$Hare <- lynx_hare_df$Hare # * 1000
population_plot1 <-
  ggplot(data = lynx_hare_df, aes(x = Lynx, y = Hare, color = Year)) +
  geom_path() +
  geom_point()
population_plot1

lynx_hare_melted_df <- melt(as.matrix(lynx_hare_df[, 2:3]));
colnames(lynx_hare_melted_df) <- c("year", "species", "population")
lynx_hare_melted_df$year <-
  lynx_hare_melted_df$year + rep(1899,
  length(lynx_hare_melted_df$year))
population_plot2 <-
  ggplot(data = lynx_hare_melted_df, aes(x = year, y = population, color = species)) +
  geom_line() +
  geom_point()
population_plot2

N <- length(lynx_hare_df$Year) - 1
ts <- 1:N
y0 <- c(lynx_hare_df$Lynx[1], lynx_hare_df$Hare[1])
y <- as.matrix(lynx_hare_df[2:(N + 1), 2:3])
lynx_hare_data <- list(N, ts, y0, y)

init_fun <- function(chain_id) {
  return(list(theta = c(0.5, 0.8, 0.025, 0.025),
              z0 = y0,
              sigma = c(0.5, 0.5)));
}

# approx solution: a1 = .5486, a2 = .8375, b1 = .0283, b2 = .0264, and

model <- stan_model("lotka-volterra-4ln.stan")
fit <- sampling(model, data = lynx_hare_data,
                chains = 1, iter = 1000,
                control = list(stepsize = 5, adapt_delta = 0.9),
		seed=123)

# fit <- optimizing(model, data = lynx_hare_data, init = init_fun)
