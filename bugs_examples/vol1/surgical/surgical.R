# BUGS surgical example (Vol 1, Example 4)
# http://www.mrc-bsu.cam.ac.uk/wp-content/uploads/WinBUGS_Vol1.pdf
library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores(logical = FALSE))
chains = 4
iter = 1000
set.seed(4711)

sourceToList = function(file){
  source(file, local = TRUE)
  d = mget(ls())
  d$file = NULL
  d
}
# Data are the same for all models
data = sourceToList("surgical.data.R")
init = rep(list(sourceToList("surgical.init.R")), chains)

# Original init
surgical = stan(file = "surgical.stan", data = data, chains = chains,
             init = init, iter = iter)
plot(surgical, ask = FALSE)

surgical_summary = summary(surgical,probs = c(.025,0.975),
                        pars = "p")$summary[,c(1,3,4,5,6)]

# Stanified, default random init
surgical_stanified = stan(file = "surgical_stanified.stan", data = data, chains = chains,
                iter = iter)
plot(surgical_stanified, ask = FALSE)
surgical_stanified_summary =
  summary(surgical_stanified,probs = c(.025,0.975),
          pars = "p")$summary[,c(1,3,4,5,6)]

# Stanified has better n_eff, no big differences otherwise
(surgical_summary)
(surgical_stanified_summary)

# Plot Figure 6 in manual
library(ggplot2)
theme_set(theme_bw())
sss = data.frame(surgical_stanified_summary[,c(1,3,4)])
sss$hospital = as.factor(paste0(LETTERS[1:length(data$n)],"(",data$r,"/",data$n,")"))
colnames(sss)[2:3] = c("min","max")
sss$hospital = reorder(sss$hospital, sss$min)
pop_mean = get_posterior_mean(surgical_stanified, pars="pop_mean")[,"mean-all chains"]
ggplot(data = sss) +
  geom_errorbarh(aes(y = hospital, x = mean, xmin = min, xmax = max), size = 1, height = .2) +
  geom_point(aes(y = hospital, x = mean), shape = 21, fill = "white", size = 4) +
  xlab("Proportion of deaths") + geom_vline(xintercept = pop_mean, col = "gray") +
  ggtitle("BUGS Vol 1, Example 4: Hospital mortality")

# Figure 7 in manual
ssr = data.frame(summary(surgical_stanified,probs = c(.025,0.975),
          pars = "ranks")$summary[,c(1,4,5)])
ssr$hospital = as.factor(paste0(LETTERS[1:length(data$n)],"(",data$r,"/",data$n,")"))
colnames(ssr)[2:3] = c("min","max")
ssr$hospital = reorder(ssr$hospital, ssr$mean)
ggplot(data = ssr) +
  geom_errorbarh(aes(y = hospital, x = mean, xmin = min, xmax = max), size = 1, height = .2) +
  geom_point(aes(y = hospital, x = mean), shape = 21, fill = "white", size = 4) +
  xlab("Rank of Hospital") +
  ggtitle("BUGS Vol 1, Example 4: Hospital mortality")

