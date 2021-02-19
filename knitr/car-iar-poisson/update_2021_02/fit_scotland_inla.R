library(INLA)
data("Scotland")
data_scotland = data.frame(y = Scotland$Counts, E = Scotland$E, X = Scotland$X/10, Region=Scotland$Region)

adj = sparseMatrix(i=node1,j=node2,x=1,dims=c(N,N))
adj =( adj + t(adj))/2 + Diagonal(N)

formula  = y ~ 1 + X + f(Region, model="besag",graph=adj,
                              hyper = list(prec= list( param=c(1,1))))

result =inla(formula, data = data_scotland, family="poisson", E = E)

phis = data.frame( mean=   result$summary.random$Region$mean, sd = result$summary.random$Region$sd, 
                   lower_quant = result$summary.random$Region$`0.025quant`, 
                   upper_quant = result$summary.random$Region$`0.975quant`,
                   median = result$summary.random$Region$`0.5quant`)
summary(result)
