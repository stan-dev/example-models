library(data.table)
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
  vals[,x:=x]
  melt(data.table(vals), id.vars="x")
}




#plot initial data
data=generate_data(fun=function(x){3+2*sin(5*x)+3*x}, sd=100, xmin=0,xmax=3, npoints=1000)
data[,pos:=x]
data[,x:=x*1e5/3.+35400000]
ggplot(data)+geom_point(aes(x,y))+geom_line(aes(x,f))#+scale_y_log10()

#fit it with stan
smf="gam_one_centered_design.stan"
sm = stan_model(file = smf)
op = optimizing(sm, data = list(N=data[,.N], K=10, y=data[,y], x=data[,x]),
                as_vector=F, hessian=F, iter=10000)
#sf = stan(file=smf, data = list(N=data[,.N], K=10, y=data[,y], x=data[,x]))
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


