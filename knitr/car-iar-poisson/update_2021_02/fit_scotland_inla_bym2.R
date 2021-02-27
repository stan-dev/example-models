library(INLA)
data("Scotland")
data_scotland = data.frame(y = Scotland$Counts, E = Scotland$E, X = Scotland$X/10, Region=Scotland$Region)

graph.scot = system.file("demodata/scotland.graph", package="INLA")

g = inla.read.graph(graph.scot)

id.singletons = c(6,8,11)
G = inla.graph2matrix(g)
G[id.singletons,] = 0
G[ , id.singletons ] = 0

#generates the disconnected graph
g_disc = inla.read.graph(G)


#sepcify the latent structure usig a fromula object
formula.bym2 = Counts ~ 1 + f(Region, model="bym2",
                           scale.model=TRUE,
                           adjust.for.con.comp=TRUE,
                           graph=g_disc,
                           hyper=list(
                               phi=list(
                                   prior = "pc",
                                   param=c (0.5 ,0.5)
                               ),
                               prec=list( prior="pc.prec",
                                         param=c(0.2 / 0.31 ,0.05) ,
                                         initial =5))
                           )
                                      
result= inla(formula.bym2, data=Scotland, E=E, family="poisson")

#phis = data.frame( mean= result$summary.random$Region$mean, sd = result$summary.random$Region$sd, 
                   #lower_quant = result$summary.random$Region$`0.025quant`, 
                   #upper_quant = result$summary.random$Region$`0.975quant`,
                   #median = result$summary.random$Region$`0.5quant`)
summary(result)
