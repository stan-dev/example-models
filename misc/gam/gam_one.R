##
# Generalized Additive Model : y ~ s(X)
# using cubic splines with negative binomial and log-link
# The spline design matrix is centered and augmented with an intercept
# Author: Yannick G Spill
# References:
# - Generalized Additive Models: an introduction with R, Simon N. Wood, CRC press (2006)
# - Eilers and Marx, Stat Sci, 1996
# This work is licensed under a Creative Commons Attribution-ShareAlike 4.0 International License
##

library(data.table)
library(parallel)
library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
library(ggplot2)
library(shinystan)
library(mgcv)

#generate data according to a function
generate_data = function(fun=sin, xmin=0, xmax=1, npoints=100, sd=1) {
  x=seq(from=xmin, to=xmax, length.out=npoints)
  f=sapply(x,fun)
  #y=rnorm(npoints, mean=f, sd=sd)
  y=rnbinom(npoints, mu=f, size=sd)
  return(data.table(x=x,f=f,y=y))
}


stan_matrix_to_datatable = function(opt, x) {
  vals=data.table(opt)
  melt(data.table(vals), id.vars="x")
}


#plot initial data
data=generate_data(fun=function(x){3+2*sin(5*x)+3*x}, sd=100, xmin=0,xmax=3, npoints=1000)
ggplot(data)+geom_point(aes(x,y))+geom_line(aes(x,f))#+scale_y_log10()

##make dump
#N=data[,.N]
#K=10
#y=data[,y]
#x=data[,x]
#dump(c("N","K","y","x"), file="gam_one.data.R")

#fit it with stan
smf="gam_one_centered_design.stan"
sm = stan_model(file = smf)
op = optimizing(sm, data = standata, as_vector=F, hessian=F, iter=10000)
#sf = stan(file=smf, data = standata, iter = 1000)
#launch_shinystan(sf)

#statistics
c(op$par$lambda,op$par$edf,op$value)

#extract weighted basis functions
#mat=stan_matrix_to_datatable(op$par$designmat, data[,x])
mat=stan_matrix_to_datatable(op$par$weighted, data[,x])
data[,pred:=op$par$pred]

#compare to gam
fit=gam(data=data, formula = y~s(x,bs="ps", m=c(2,2), k=10), family=nb())
data[,gam:=fit$fitted.values]

#plot result
ggplot(data)+geom_point(aes(x,y),alpha=0.2)+geom_line(aes(x,f),colour="black")+
  geom_line(data=mat, aes(x,value,colour=variable))+
  geom_line(aes(x,gam),colour="blue")+
  geom_line(aes(x,pred),colour="red")#+scale_y_log10()


