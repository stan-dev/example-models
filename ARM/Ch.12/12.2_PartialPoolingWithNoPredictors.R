library(rstan)
library(lme4)
library(ggplot2)

srrs2 <- read.table ("ARM/Ch.12/srrs2.dat", header = TRUE, sep = ",")
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
  county[county.name == uniq[i]] <- i
}

 # no predictors
ybarbar <- mean(y)

sample.size <- as.vector (table (county))
sample.size.jittered <- sample.size * exp(runif(J, -.1, .1))
cty.mns <- tapply(y, county, mean)
cty.vars <- tapply(y, county, var)
cty.sds <- mean(sqrt(cty.vars[!is.na(cty.vars)])) / sqrt(sample.size)
cty.sds.sep <- sqrt(tapply(y, county, var) / sample.size)

 # lme4 varying-intercept model, no predictors
M0 <- lmer(y ~ 1 + (1 | county))
summary(M0)
M0_coef <- coef(M0)$county[, 1]

dataList.1 <- list(N = length(y), J = J, y = y, county = county)
radon_intercept.sf1 <- stan(file = 'ARM/Ch.12/radon_intercept.stan', data = dataList.1,
                            iter = 200, chains = 4, control = list(stepsize = 0.05))
print(radon_intercept.sf1)
post <- extract(radon_intercept.sf1)
mean.a <- rep (NA, J)
sd.a <- rep (NA, J)
for (n in 1:J) {
  mean.a[n] <- mean(post$a[ ,n])
  sd.a[n] <- sd(post$a[ ,n])
}

# comparing lme4 and stan posterior mean estimates
sqrt(mean(M0_coef - mean.a)^2)

## Figure 12.1 (a)
frame1 = data.frame(x1=sample.size.jittered,y1=cty.mns,x2=sample.size.jittered[36],y2=cty.mns[36])
limits <- aes(ymax=cty.mns + cty.sds,ymin=cty.mns - cty.sds)
p1 <- ggplot(frame1,aes(x=x1,y=y1)) +
      geom_point(aes(x=x2,y=y2),shape=1,size=30) +
      scale_y_continuous("Avg. Log Radon in County j") +
      scale_x_log10("Sample Size in County j") +
      theme_bw() +
      geom_pointrange(limits) +
      labs(title="No Pooling")
print(p1)

## Figure 12.1 (b)
frame2 = data.frame(x1=sample.size.jittered,y1=mean.a,x2=sample.size.jittered[36],y2=mean.a[36])
limits <- aes(ymax=mean.a+sd.a, ymin=mean.a-sd.a)
p2 <- ggplot(frame2,aes(x=x1,y=y1)) +
      geom_point(aes(x=x2,y=y2),shape=1,size=30) +
      scale_y_continuous("Avg. Log Radon in County j") +
      scale_x_log10("Sample Size in County j") +
      theme_bw() +
      geom_pointrange(limits) +
      labs(title="Multilevel Model")
print(p2)
