
dp_dt <- function(q, k, star) {
  r = sqrt((q - star) %*% (q - star))
  - k * (q - star) / c(r^3)
}

solve_trajectory <- function(q0, p0, dt, k, m, n_obs, ts, star = c(0, 0)) {
  q_s = array(NA, c(n_obs, 2))
  
  q = q0
  p = p0
  t = 0
  n = 1
  
  eps <- dt / 10
  
  while (t <= ts[n_obs]) {

    if (abs(t - ts[n]) < eps) {
      q_s[n, ] = q
      n = n + 1
    }

    p = p + 0.5 * dt * dp_dt(q, k, star)
    q = q + p / m * dt
    p = p + 0.5 * dt * dp_dt(q, k, star)

    t = t + dt
  }

  q_s
}


ppc_plot <- function(fit, chains, pred = "qx_pred", pars = "k",
                     data_pred) {
  qx_pred <- rstan::extract(fit, pars = c(pred), permuted = FALSE)
  
  for (i in 1:chains) { 
    qx_pred_chain <- qx_pred[, i, ]
    pred_new <- as.data.frame(qx_pred_chain) %>% gather(factor_key = TRUE) %>%  
      group_by(key) %>%
      summarize(lb = quantile(value, probs = 0.05), 
                median = quantile(value, probs = 0.5),
                ub = quantile(value, probs = 0.95)) %>% 
      bind_cols(data_pred)
    
    if (i == 1) pred <- pred_new
    if (i != 1) pred <- bind_rows(pred, pred_new) 
  }
  
  estimated_k <- rep(NA, chains)
  for (i in 1:chains) estimated_k[i] <-
    paste0(pars, " = ", summary(r_fit, pars = pars)$c_summary[1, 1, i]) 
  pred$chain <-  rep(estimated_k, each = 40)

  p1 <- ggplot(pred, aes(x = t, y = q)) +
    geom_point() +
    geom_line(aes(x = t, y = median)) + 
    geom_ribbon(aes(ymin = lb, ymax = ub), alpha = 0.25) + 
    facet_wrap(~chain) + theme_bw() 

  p1
}


ppc_plot2D <- function(fit, pred = c("qx_pred", "qy_pred"), pars = "k", data_pred,
                       plot_star = FALSE) {
  qx_pred <- fit$draws(variables = pred[1])
  qy_pred <- fit$draws(variables = pred[2])

  chains <- dim(qx_pred)[2]
  if (plot_star) {
    # extract the median estimate for the star's position (for each chain)
    star <- array(NA, c(2, chains))
    for (c in 1:chains) star[, c] <-
      pull(summary(as_draws_df(fit$draws(variables = "star")[, c, ])),
           "median")
  }

  for (i in 1:chains) {
    qx_pred_chain <- as_draws_df(qx_pred[, i, ])
    qy_pred_chain <- as_draws_df(qy_pred[, i, ])
    
    pred_new <- rbind(summary(qx_pred_chain)[, c(3, 6, 7)],
                      summary(qy_pred_chain)[, c(3, 6, 7)])

    if (i == 1) pred <- pred_new
    if (i != 1) pred <- bind_rows(pred, pred_new)
  }

  sequence_x <- c()
  sequence_y <- c()
  for (c in 1:chains) {
    start <- (c - 1) * 80 + 1
    sequence_x <- c(sequence_x, start:(start + 39))

    start <- (c - 1) * 80 + 41
    sequence_y <- c(sequence_y, start:(start + 39))
  }

  plot_data <- data.frame(qx_pred = pred$median[sequence_x], 
                          qy_pred = pred$median[sequence_y])
  plot_data$chains <- rep(paste0("chain ", 1:chains), each = 40)
  plot_data <- cbind(plot_data, data_pred)

  p1 <- ggplot(plot_data, aes(x = qx, y = qy)) +
    geom_point(shape = 3, alpha = 0.5) + theme_bw() +
    geom_path(aes(x = qx_pred, y = qy_pred)) +
    facet_wrap(~chains)
  
  if (plot_star) {
    data_add <- as.data.frame(t(star))
    names(data_add) <- c("star_x", "stary_y")
    data_add$chains <- paste0("chain ", 1:8)
    names(data_add) <- c("star_x", "star_y", "chains")

    p1 <- p1 + geom_point(data = data_add,
                          aes(x = star_x, y = star_y),
                          shape = 8)
  }
  
  p1
}

