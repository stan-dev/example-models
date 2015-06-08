## ---- fit-vs-sim-function ----
fit_vs_sim <- function(var_name, interval_width = 0.9) {
  ss <- extract(fit);
  theta <- get(var_name);
  theta_hat <- c();
  theta_lower <- c();
  theta_upper <- c();
  for (i in 1:length(theta)) {
    ss_theta_i <- ss[[var_name]][,i];
    theta_hat[i] <- mean(ss_theta_i);
    theta_lower[i] <- quantile(ss_theta_i, probs=0.5 - interval_width/2);
    theta_upper[i] <- quantile(ss_theta_i, probs=0.5 + interval_width/2);
  }
  df <- data.frame(theta, theta_hat, theta_lower, theta_upper);
  ggp <-  
    ggplot(df, aes(x=theta, y=theta_hat)) +
    geom_point() +
    geom_errorbar(aes(ymin=theta_lower, ymax=theta_upper), width=.05) +
    geom_abline(intercept=0, slope=1, colour="red") +
    xlab(paste(var_name,"(simulated)")) +
    ylab(paste(var_name," (posterior mean and ", 100 * interval_width,
               "% interval)", sep="")) +
    ggtitle(paste("Posterior Estimate vs. Simulation: ", var_name));
  return(ggp);
}
