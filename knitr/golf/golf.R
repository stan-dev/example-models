library("arm")
library("rstan")
library("rstanarm")
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

## Set up the data
golf <- read.table("golf.txt", header=TRUE, skip=2)
x <- golf$x
y <- golf$y
n <- golf$n
J <- length(y)
r <- (1.68/2)/12
R <- (4.25/2)/12
se <- sqrt((y/n)*(1-y/n)/n)

## Fit the model
golf_data <- list(x=x, y=y, n=n, J=J)
fit_trig <- stan("golf1.stan", data=golf_data)
print(fit_trig)

## Post-processing
sims1 <- extract(fit1)
sigma_hat <- median(sims1$sigma)

## Naive logistic regression
fit_logistic <- stan("golf_logistic.stan", data=golf_data)
print(fit_logistic)

## Post-processing
sims_logistic <- extract(fit_logistic)
a_hat <- median(sims_logistic$a)
b_hat <- median(sims_logistic$b)

## Plots
pdf("golf0.pdf", height=5, width=7)
par(mar=c(3,3,2,1), mgp=c(1.7,.5,0), tck=-.02)
plot(x, y/n, xlim=c(0, 1.1*max(x)), ylim=c(0, 1.02), xaxs="i", yaxs="i", pch=20, bty="l", xlab="Distance from hole (feet)", ylab="Probability of success", main="Data on putts in pro golf")
segments(x, y/n + se, x, y/n-se, lwd=.5)
text(x + .4, y/n + se + .02, paste(y, "/", n,sep=""), cex=.6, col="gray40") 
dev.off()

pdf("golf1.pdf", height=5, width=7)
par(mar=c(3,3,2,1), mgp=c(1.7,.5,0), tck=-.02)
plot(x, y/n, xlim=c(0, 1.1*max(x)), ylim=c(0, 1.02), xaxs="i", yaxs="i", pch=20, bty="l", xlab="Distance from hole (feet)", ylab="Probability of success", main="Fitted logistic regression")
segments(x, y/n + se, x, y/n-se, lwd=.5)
curve(invlogit(a_hat + b_hat*x), from=0, to=1.1*max(x), add=TRUE)
text(10.6, .57, paste("Logistic regression,\n    a = ", fround(a_hat, 2), ", b = ", fround(b_hat, 2), sep=""))
dev.off()

pdf("golf2.pdf", height=5, width=7)
par(mar=c(3,3,2,1), mgp=c(1.7,.5,0), tck=-.02)
plot(x, y/n, xlim=c(0, 1.1*max(x)), ylim=c(0, 1.02), xaxs="i", yaxs="i", pch=20, bty="l", xlab="Distance from hole (feet)", ylab="Probability of success", main="Custom nonlinear model fit in Stan")
segments(x, y/n + se, x, y/n-se, lwd=.5)
x_grid <- seq(R-r, 1.1*max(x), .01)
p_grid <- 2*pnorm(asin((R-r)/x_grid) / sigma_hat) - 1
lines(c(0, R-r, x_grid), c(1, 1, p_grid))
text(18.5, .26, paste("Geometry-based model,\n sigma = ", fround(sigma_hat*180/pi, 1), " degrees", sep=""))
dev.off()


pdf("golf2a.pdf", height=5, width=7)
par(mar=c(3,3,2,1), mgp=c(1.7,.5,0), tck=-.02)
plot(x, y/n, xlim=c(0, 1.1*max(x)), ylim=c(0, 1.02), xaxs="i", yaxs="i", pch=20, bty="l", xlab="Distance from hole (feet)", ylab="Probability of success", main="Custom nonlinear model fit in Stan")
segments(x, y/n + se, x, y/n-se, lwd=.5)
x_grid <- seq(R-r, 1.1*max(x), .01)
n_sims <- length(sims1$sigma)
for (s in sample(n_sims, 20)){
  p_grid <- 2*pnorm(asin((R-r)/x_grid) / sims1$sigma[s]) - 1
  lines(c(0, R-r, x_grid), c(1, 1, p_grid), lwd=0.5)
}
text(18.5, .26, "Geometry-based model,\n post draws of sigma")
dev.off()

pdf("golf3.pdf", height=5, width=7)
par(mar=c(3,3,2,1), mgp=c(1.7,.5,0), tck=-.02)
plot(x, y/n, xlim=c(0, 1.1*max(x)), ylim=c(0, 1.02), xaxs="i", yaxs="i", pch=20, bty="l", xlab="Distance from hole (feet)", ylab="Probability of success", main="Two models fit to the golf putting data")
segments(x, y/n + se, x, y/n-se, lwd=.5)
curve(invlogit(a_hat + b_hat*x), from=0, to=1.1*max(x), add=TRUE)
x_grid <- seq(R-r, 1.1*max(x), .01)
p_grid <- 2*pnorm(asin((R-r)/x_grid) / sigma_hat) - 1
lines(c(0, R-r, x_grid), c(1, 1, p_grid))
text(10.3, .58, "Logistic regression")
text(18.5, .24, "Geometry-based model")
dev.off()
