# clears workspace: 
rm(list=ls()) 

library(rstan)

model <- "
// Cheating Latent Mixture Model
data { 
  int<lower=1> n;
  int<lower=1> p;
  int<lower=1,upper=n> k[p];
  int<lower=0,upper=1> truth[p];
}
parameters {
  real<lower=0,upper=1> phi;
  real<lower=0,upper=1> mubon;
  real<lower=0> mudiff;
  real<lower=5,upper=50> lambdabon;
  real<lower=5,upper=50> lambdache;
  matrix<lower=0,upper=1>[2,p] theta;
} 
transformed parameters {
  vector[2] lp_parts[p];
  vector<lower=0>[2] alpha;
  vector<lower=0>[2] beta;
  real<lower=0,upper=1> muche;
    
  // Additivity on Logit Scale
  muche <- inv_logit(logit(mubon) + mudiff);
  
  // Transformation to Group Mean and Precision
  alpha[1] <- mubon * lambdabon;
  beta[1] <- lambdabon * (1 - mubon);
  alpha[2] <- muche * lambdache;
  beta[2]  <- lambdache * (1 - muche);
    
  // Data are Binomial with Rate Given by 
  // Each Person’s Group Assignment
  for (i in 1:p) {
    lp_parts[i,1] <- log1m(phi) + binomial_log(k[i], n, theta[1,i]);
    lp_parts[i,2] <- log(phi) + binomial_log(k[i], n, theta[2,i]);
  } 
}
model {
  // Priors
  mubon ~ beta(1, 1);  // can be removed
  mudiff ~ normal(0, 1 / sqrt(.5))T[0,];  // Constrained to be Positive
  // Relatively Uninformative Prior on Base Rate
  phi ~ beta(5, 5);
    
  theta[1] ~ beta(alpha[1], beta[1]);
  theta[2] ~ beta(alpha[2], beta[2]);  
  
  for (i in 1:p)  
    increment_log_prob(log_sum_exp(lp_parts[i]));    
}
generated quantities {
  int<lower=0,upper=1> z[p];
  real pc;
  vector[p] pct;
  
  for (i in 1:p) {
    vector[2] prob;
    prob <- softmax(lp_parts[i]);
    // Each Person Belongs to One of Two Latent Groups
    z[i] <- bernoulli_rng(prob[2]);
    // Correct Count
    pct[i] <- if_else(z[i] == truth[i], 1, 0);
  }
  pc <- sum(pct);
}"

cheat.dat  <- read.table("cheat.csv", header=F, sep=",")
cheatt.dat <- read.table("cheatt.csv", header=F, sep="")
truth <- cheatt.dat$V1  # truth = 1 if cheater
k <- apply(cheat.dat, 1, sum)  # total correct per participant
p <- length(k)  # number of people
n <- 40         # total trials

data <- list(p=p, k=k, n=n, truth=truth) # To be passed on to Stan

myinits <- list(
  list(mudiff=.1, phi=.5, mubon=.5, lambdabon=30, lambdache=25, 
       theta=matrix(rep(.5, 2 * p), 2, p)),
  list(mudiff=.15, phi=.5, mubon=.5, lambdabon=25, lambdache=30,
       theta=matrix(rep(.5, 2 * p), 2, p))) 

# Parameters to be monitored:
parameters <- c("theta", "z", "mubon", "lambdabon", "muche", "lambdache", 
                "mudiff", "phi", "alpha", "beta", "pc") 

# The following command calls Stan with specific options.
# For a detailed description type "?rstan".
samples <- stan(model_code=model,   
                data=data, 
                init=myinits,  # If not specified, gives random inits
                pars=parameters,
                iter=5000, 
                chains=2, 
                thin=1,
                # warmup = 100,  # Stands for burn-in; Default = iter/2
                # seed = 123  # Setting seed; Default is random seed
)
# Now the values for the monitored parameters are in the "samples" object, 
# ready for inspection.

print(samples)
pc <- extract(samples)$pc / p  # to get proportion correct
mean(pc)

# plot 6.9
#make the two panel plot:
windows(width=8,height=6) #this command works only under Windows!
layout(matrix(c(1,2),2,1))
layout.show(2)
par(cex.main = 1.5, mar = c(5, 6, 4, 5) + 0.1, mgp = c(3.5, 1, 0), cex.lab = 1.5,
    font.lab = 2, cex.axis = 1.3, bty = "n", las=1)
bins <- c(-1:n)+.5
bonafide <- hist(k[truth==0], breaks=bins, plot=F)$counts
cheat    <- hist(k[truth==1], breaks=bins, plot=F)$counts

counts <- rbind(bonafide, cheat)
barplot(counts, main=" ", xlab=" ", col=c("grey","white"),
  legend.text = c("Bona Fide","Cheater"), args.legend = list(x="topleft"),
  beside=TRUE, axes=F)
# bottom panel:
par(cex.main = 1.5, mar = c(5, 6, 4, 5) + 0.1, mgp = c(3.5, 1, 0), cex.lab = 1.5,
    font.lab = 2, cex.axis = 1.3, bty = "n", las=1)
pc.line <- array()
for (i in 1:41) {
  pc.line[i] <- mean((k>=(i-1))==truth)
}

dev.new() # so the plot below does not overwrite the plot above

plot(c(0:40), pc.line, type="l", lwd=2, xlim=c(0,40), ylim=c(0.4,1), 
     xlab="Number of Items Recalled Correctly", 
     ylab=" ", axes=F)
axis(1, at=c(0,seq(from=5,by=5,to=40)))
axis(2, at=c(.5,.75,1))
par(las=0)
mtext("Prop. Correct",side=2, line=2.5,cex=1.5)
# Now add the distribution:
pc.dens <- density(pc)
polygon(c(0,pc.dens$y,0,0), c(pc.dens$x[1]-.01,pc.dens$x,pc.dens$x[1]+.01,
                              pc.dens$x[1]-.01), col="green")

# plot 6.10
windows()
par(cex.main = 1.5, mar = c(5, 6, 4, 5) + 0.1, mgp = c(3.5, 1, 0), cex.lab = 1.5,
    font.lab = 2, cex.axis = 1.3, bty = "n", las=1)
plot(k,summary(samples)$summary[237:354, 1],ylim=c(0,1),xlim=c(0,n), 
     xlab= "Number of Items Recalled Correctly", ylab="Cheater Classification", 
     lwd=2, pch=4) 
# in the code, z=0 is bonafide and z=1 is cheating
# so z gives the prob of being assigned to cheating group

