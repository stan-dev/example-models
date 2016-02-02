hits_2006_df <- read.csv("baseball-hits-2006.csv", sep=",", comment.char="#");
y <- hits_2006_df$y;
K <- hits_2006_df$K;
N <- length(y);
K_new <- rep(0, N);
y_new <- rep(0.0, N);

fit_hits_2006_pool <- stan("pool.stan",
                       data = c("N", "K", "y", "K_new", "y_new"));

fit_hits_2006_no_pool <- stan("no-pool.stan",
                         data = c("N", "K", "y", "K_new", "y_new"));

fit_hits_2006_hier <- stan("hier.stan",
                      data = c("N", "K", "y", "K_new", "y_new"));

fit_hits_2006_hier_logit <- stan("hier-logit.stan",
                            data = c("N", "K", "y", "K_new", "y_new"));

pars_to_print <- c("p_max", "p_mean", "p_sd",
                   "theta",
                   "some_ability_gt_350");
pars_to_print2 <- c(pars_to_print,
                    "is_best");
                   
print(fit_hits_2006_pool, c(pars_to_print, "phi"));
print(fit_hits_2006_no_pool, pars_to_print2);
print(fit_hits_2006_hier, c(pars_to_print2, "phi", "kappa"));
print(fit_hits_2006_hier_logit, c(pars_to_print2, "mu", "sigma"));
