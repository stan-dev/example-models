library(rstan)
library(ggplot2)
library(lme4)
library(arm)  # for the disaplay() function

# Data are at http://www.stat.columbia.edu/~gelman/arm/examples/radon

srrs2 <- read.table ("ARM/Ch.12/srrs2.dat", header=T, sep=",")
mn <- srrs2$state=="MN"
radon <- srrs2$activity[mn]
log.radon <- log (ifelse (radon==0, .1, radon))
floor <- srrs2$floor[mn]       # 0 for basement, 1 for first floor
n <- length(radon)
y <- log.radon
x <- floor

# get county index variable
county.name <- as.vector(srrs2$county[mn])
uniq <- unique(county.name)
J <- length(uniq)
county <- rep (NA, J)
for (i in 1:J){
  county[county.name==uniq[i]] <- i
}

 # no predictors
ybarbar = mean(y)

sample.size <- as.vector (table (county))
sample.size.jittered <- sample.size*exp (runif (J, -.1, .1))
cty.mns = tapply(y,county,mean)
cty.vars = tapply(y,county,var)
cty.sds = mean(sqrt(cty.vars[!is.na(cty.vars)]))/sqrt(sample.size)
cty.sds.sep = sqrt(tapply(y,county,var)/sample.size)

## Varying-intercept model w/ no predictors
dataList.1 <- list(N = length(y), y = y, county = county, J = J)
radon_intercept.sf1 <- stan(file = 'ARM/Ch.12/radon_intercept.stan',
                            data = dataList.1,
                            iter = 500, chains = 4,
                            control = list(stepsize = 0.05))
print(radon_intercept.sf1)

## Including x as a predictor
dataList.2 <- list(N = length(y), y = y, x = x, county = county, J = J)
radon_no_pool.sf1 <- stan(file = 'ARM/Ch.12/radon_no_pool.stan',
                          data = dataList.2,
                          iter = 500, chains = 4,
                          control = list(stepsize = 0.05))
print(radon_no_pool.sf1)

M1 <- extract(radon_no_pool.sf1)
M1.alpha <- colMeans(M1$a)
M1.beta <- mean(M1$beta)

# equivalent model in lmer
M1_lmer <- lmer (y ~ x + (1 | county))
display(M1_lmer)
coef(M1_lmer)
fixef(M1_lmer)
ranef(M1_lmer)
se.fixef(M1_lmer)
se.ranef(M1_lmer)

# 95% CI for the slope
M1.beta + c(-2,2) * sd(M1$beta)

# same as above, but with lmer
fixef(M1_lmer)["x"] + c(-2,2) * se.fixef(M1_lmer)["x"]

# 95% CI for the intercept in county 26
M1.alpha[26] + c(-2, 2) * sd(M1$a[, 26])

# same as above, but with lmer
coef(M1_lmer)$county[26, 1] + c(-2, 2) * se.ranef(M1_lmer)$county[26]

# to plot Figure 12.4
a.hat.M1 <- coef(M1_lmer)$county[, 1]                # 1st column is the intercept
b.hat.M1 <- coef(M1_lmer)$county[, 2]                # 2nd element is the slope

x.jitter <- x + runif(n,-.05,.05)
display8 <- c (36, 1, 35, 21, 14, 71, 61, 70)  # counties to be displayed
y.range <- range (y[!is.na(match(county,display8))])

radon.data <- data.frame(y, x.jitter, county)
radon8.data <- subset(radon.data, county %in% display8)
radon8.data$county.name <- radon8.data$county
radon8.data$county.name <- factor(radon8.data$county.name,
                                  levels = c("36","1","35","21","14","71","61","70"),
                                  labels = c("LAC QUI PARLE", "AITKIN", "KOOCHICHING",
                                             "DOUGLAS", "CLAY", "STEARNS", "RAMSEY", "ST LOUIS"))
radon8.data$pooled.int <- fixef(M1_lmer)[1]
radon8.data$pooled.slope <- fixef(M1_lmer)[2]
radon8.data$unpooled.int <- M1$a[radon8.data$county]
radon8.data$unpooled.slope <- M1.beta
radon8.data$coef.int <- a.hat.M1[radon8.data$county]
radon8.data$coef.slope <- b.hat.M1[radon8.data$county]

p1 <- ggplot(radon8.data, aes(x.jitter, y)) +
    geom_jitter(position = position_jitter(width = .05, height = 0)) +
    scale_x_continuous(breaks=c(0,1), labels=c("0", "1")) +
    geom_abline(aes(intercept = pooled.int, slope = pooled.slope), linetype = "dashed") +
    geom_abline(aes(intercept = unpooled.int, slope = unpooled.slope), size = 0.25) +
    geom_abline(aes(intercept = coef.int, slope = coef.slope), size = 0.25) +
    facet_wrap(~ county.name, ncol = 4)
print(p1)

## Multilevel model ests vs. sample size (plot on the right on figure 12.3)
a.se.M1 <- se.coef(M1_lmer)$county
sample.size <- as.vector (table (county))
sample.size.jittered <- sample.size*exp (runif (J, -.1, .1))
max1 <- a.hat.M1+a.se.M1
min1 <- a.hat.M1-a.se.M1
frame3 <- data.frame(x1 = sample.size.jittered,
                     y1 = colMeans(M1$a),
                     max1 = max1[, 1],
                     min1 = min1[, 1])

p2 <- ggplot(frame3, aes(x=x1,y=y1,ymin= min1,ymax=max1)) +
      geom_linerange() +
      geom_point() +
      scale_y_continuous("estimated intercept alpha (no pooling)") +
      scale_x_log10("Sample Size in County j") +
      theme_bw() +
      geom_hline(aes(yintercept=fixef(M1_lmer)[1]))
print(p2)
