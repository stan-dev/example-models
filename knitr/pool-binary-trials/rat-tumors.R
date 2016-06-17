rat_tumors_df <- read.csv("rat-tumors.csv", sep=",", comment.char="#");
y <- rat_tumors_df$y;
K <- rat_tumors_df$N;
N <- length(y);
K_new <- rep(0, N);
y_new <- rep(0.0, N);

fit_rats_pool <- stan("pool.stan",
                       data = c("N", "K", "y", "K_new", "y_new"));

fit_rats_no_pool <- stan("no-pool.stan",
                         data = c("N", "K", "y", "K_new", "y_new"));

fit_rats_hier <- stan("hier.stan",
                      data = c("N", "K", "y", "K_new", "y_new"));

fit_rats_hier_logit <- stan("hier-logit.stan",
                            data = c("N", "K", "y", "K_new", "y_new"));

pars_to_print <- c("p_max", "p_mean", "p_sd",
                   "theta[1]", "theta[33]", "theta[54]", "theta[70]",
                   "some_ability_gt_350");
pars_to_print2 <- c(pars_to_print,
                    "is_best[1]", "is_best[33]", "is_best[55]", "is_best[70]");
                   
print(fit_rats_pool, c(pars_to_print, "phi"));
print(fit_rats_no_pool, pars_to_print2);
print(fit_rats_hier, c(pars_to_print2, "phi", "kappa"));
print(fit_rats_hier_logit, c(pars_to_print2, "mu", "sigma"));
