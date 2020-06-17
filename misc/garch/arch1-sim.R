T=100
r=c(rep(0,T))
for (n in 1:T) {
  r[n]=rnorm(1, mean = 0, sd = 0.01*n)
}