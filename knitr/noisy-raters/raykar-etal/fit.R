

init_fun <- function() {
  return(list(alpha=rep(0.75,R), beta=rep(0.75,R), w=rep(0,D), w0=0));
}


library('rstan')
fit <- stan("raykar-marginal.stan",
            data = c("N","R","D","x","y"),
            init = function(n) list(alpha = rep(0.8, R), beta = rep(0.8, R), w = rep(0, D), w0 = 0))
