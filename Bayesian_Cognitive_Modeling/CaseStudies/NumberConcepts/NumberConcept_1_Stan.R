# clears workspace: 
rm(list=ls()) 
library(rstan)

model <- "
// Knower Level Model Applied to Give-N Data
data { 
  int<lower=1> ns;
  int<lower=1> nz;
  int<lower=1> gn;
  int gnq[ns];
  int gq[ns,21];  // no. columns = max(gnq)
  int ga[ns,21];  // no. columns = max(gnq)
}
parameters {
  vector<lower=0,upper=1>[gn] pitmp;
  real<lower=1,upper=1000> v;
} 
transformed parameters {
  simplex[gn] pi;
  simplex[gn] npiprime[nz,gn];
  vector[nz] lp_parts[ns];

  // Base rate
  pi <- pitmp / sum(pitmp);
  
  // Model
  for (i in 1:nz) {
    for (j in 1:gn) {
      vector[gn] piprime;
      for (k in 1:gn) {
        real ind1;
        real ind2;
        real ind3;
        real ind4;
        real ind5;
        
        // Will be 1 if Knower-Level (i.e, i-1) is Same or Greater than Answer
        ind1 <- step((i - 1) - k);
        // Will be 1 for the Possible Answer that Matches the Question
        ind2 <- k == j;
        // Will be 1 for 0-Knowers
        ind3 <- i == 1;
        // Will be 1 for HN-Knowers
        ind4 <- i == nz;
        ind5 <- ind3 + ind4 * (2 + ind2) 
                + (1 - ind4) * (1 - ind3) * (ind1 * ind2 + ind1 + 1);
                
        if (ind5 == 1)
          piprime[k] <- pi[k];
        else if (ind5 == 2)  
          piprime[k] <- 1 / v * pi[k];
        else if (ind5 == 3)  
          piprime[k] <- v * pi[k];
      }  
      for (k in 1:gn)
        npiprime[i,j,k] <- piprime[k] / sum(piprime);
    }
  }
  
  for (i in 1:ns) {
    for (m in 1:nz) {
      real lp_parts_tmp;
      lp_parts_tmp <- 0;
      
      // Probability a z[i]-Knower Will Answer ga[i,j] to Question gq[i,j]
      // is a Categorical Draw From Their Distribution over the 1:gn Toys
      for (j in 1:gnq[i])
        lp_parts_tmp <- lp_parts_tmp
                         + categorical_log(ga[i,j], npiprime[m,gq[i,j]]);

      lp_parts[i,m] <- log(1.0 / nz) + lp_parts_tmp;
    }  
  }
}
model {
  for (i in 1:ns)
    increment_log_prob(log_sum_exp(lp_parts[i]));
}
generated quantities {
  vector[nz] prob;
  int z[ns];
  int predga[ns,gn];
  int predz[nz,gn];
  int predpi;
  
  for (i in 1:ns) {
    prob <- softmax(lp_parts[i]);
    z[i] <- categorical_rng(prob);
  }
  
  // Posterior Predictive
  for (i in 1:ns)
    for (j in 1:gn)
      predga[i,j] <- categorical_rng(npiprime[z[i],j]);
  
  // Posterior Prediction For Knower Levels
  for (i in 1:nz)
    for (j in 1:gn)
      predz[i,j] <- categorical_rng(npiprime[i,j]);
  
  predpi <- categorical_rng(pi);
}"

load("fc_given.RData")  # Load all data for the model

# to be passed on to Stan
data <- c("ns", "gnq", "gn", "ga", "gq", "nz") 

myinits <- list(
  list(v=2, pitmp=rep(1 / gn, gn)),
  list(v=2, pitmp=rep(1 / gn, gn)))

# parameters to be monitored:  
parameters <- c("predga", "predz", "predpi", "v", "z")

# The following command calls Stan with specific options.
# For a detailed description type "stan".
samples <- stan(model_code=model,   
                data=data, 
                init=myinits,  # If not specified, gives random inits
                pars=parameters,
                iter=600, 
                chains=2, 
                thin=1,
                warmup = 100,  # Stands for burn-in; Default = iter/2
                # seed = 123  # Setting seed; Default is random seed
)

zSamples <- extract(samples)$z
predpi <- extract(samples)$predpi
predz <- extract(samples)$predz
predga <- extract(samples)$predga

#### Figure 19.2 ####
windows(7,4)
par(mar=c(3, 2, 1, 1) + .1, mgp=c(1.3, 0.2, 0), cex.lab=1.2)
barplot(table(predpi), ylim=c(0, max(table(predpi)) * 1.2), col="black", 
        yaxt="n", xlab="Number", ylab="")
box()
title(ylab="Probability", line=.2)
axis(3, at=seq(.7, 17.5, by=1.2), label=FALSE, tck = 0.02)

#### Figure 19.3 ####
windows(9, 6)
par(mfrow=c(4, 5), mar=c(1, 0, 2, 2) + .1, oma=c(2.6, 2.8, 1, 0),
    mgp=c(1.5, 0.15, 0))
for (i in 1:ns) {
  zTable <- table(factor(zSamples[, i], levels=as.character(1:6)))
  
  barplot(zTable, col="black", xaxt="n", yaxt="n", main=paste("Child", i),
          ylim=c(0, length(zSamples[, 1]))) 
  
  if (i == 16)
    axis(1, at=seq(.7, 7.5, by=1.2), label=c("P", "1", "2", "3", "4", "C"), tck = 0.1)
  else
    axis(1, at=seq(.7, 7.5, by=1.2), label=FALSE, tck = 0.1)
  
  axis(3, at=seq(.7, 7.5, by=1.2), label=FALSE, tck = 0.1)
  box()
}
mtext("Knower", side=1, line=0.2, at=.055, adj=0, outer=TRUE)  
mtext("Level", side=1, line=1.2, at=.065, adj=0, outer=TRUE)  
mtext("Posterior", side=2, line=1.2, at=.07, adj=0, outer=TRUE)  
mtext("Mass", side=2, line=0.2, at=.09, adj=0, outer=TRUE)

#### Figure 19.4 ####
subjlist = c(15, 2, 4, 3, 10, 20)
sc=2
sc2=1
sc3=2
cm=-.4
cb=.99
shadedSize <- 2

windows(9, 6)
par(mfrow=c(2, 3), mar=c(2, 2, 2, 1) + .1, oma=c(2, 3, 0, 0), mgp=c(.2, .2, 0),
    cex.lab=1)
for (s in subjlist) {
  plot(NA, xlim=c(.5, 15.5), ylim=c(.5, 15.5), axes=FALSE, cex.main=1, ylab="", 
       main=paste("Child", s), xlab="")
  axis(1, at=c(1, 2, 3, 4, 5, 8, 10), tck=0.01)
  axis(2, at=1:15, las=1, tck=0.01)
  axis(3, at=c(1, 2, 3, 4, 5, 8, 10), labels=FALSE, tck=0.01)
  axis(4, at=1:15, labels=FALSE, tck=0.01)
  box()
  
  for (i in 1:gn) {
    count <- hist(predga[, s, i], plot=FALSE, breaks=seq(0.5, 15.5))$counts
    count <- count / max(count)
    
    for (j in 1:gn)
        points(i, j, pch=15, col=gray(min(1, max(0, (cm * count[j] + cb)))), 
               cex=shadedSize)
    
    tmp <- gq[s, ]
    match <- which(tmp == i)
    if (length(match) != 0) {
      for (j in 1:gn) {
        count <- sum(ga[s, match] == j)
        count <- count / length(match)
        if (count > 0)
          points(i, j, pch=22, cex=sc * sqrt(count), lwd=sc3 * count) 
      }
    }
  }
}
mtext("Question", side=1, line=0.2, outer=TRUE)   
mtext("Answer", side=2, line=0.2, outer=TRUE)   

#### Figure 19.5 ####
mainTitle <- c("PN-Knower", "One-Knower", "Two-Knower", "Three-Knower",
               "Four-Knower", "CP-Knower")
dm <- array(0, dim=c(gn, gn, nz))
mv <- c()
for (i in 1:ns) {  # computing mode of z for each child
  uz <- unique(zSamples[, i])
  mv[i] <- uz[which.max(tabulate(match(zSamples[, i], uz)))]
}

windows(9, 6)
par(mfrow=c(2, 3), mar=c(2, 2, 2, 1) + .1, oma=c(2, 3, 0, 0), mgp=c(.2, .2, 0),
    cex.lab=1)
for (z in 1:nz) {
  plot(NA, xlim=c(.5, 15.5), ylim=c(.5, 15.5), axes=FALSE, cex.main=1, ylab="", 
       main=mainTitle[z], xlab="")
  axis(1, at=c(1, 2, 3, 4, 5, 8, 10), tck=0.01)
  axis(2, at=1:15, las=1, tck=0.01)
  axis(3, at=c(1, 2, 3, 4, 5, 8, 10), labels=FALSE, tck=0.01)
  axis(4, at=1:15, labels=FALSE, tck=0.01)
  box()
  
  for (i in 1:gn) {
    count <- hist(predz[, z, i], plot=FALSE, breaks=seq(0.5, 15.5))$counts
    count <- count / max(count)
    
    for (j in 1:gn)
      points(i, j, pch=15, col=gray(min(1, max(0, (cm * count[j] + cb)))), 
             cex=shadedSize)
  }
  match <- which(mv == z)
  
  for (i in 1:length(match))
    for (j in 1:gnq[match[i]])
      dm[gq[match[i], j], ga[match[i], j], z] <- dm[gq[match[i], j], ga[match[i], j], z] + 1
  
  for (i in 1:gn) {
    if (sum(dm[i, , z]) == 0)
      count <- rep(0, 15)
    else
      count <- dm[i, , z] / sum(dm[i, , z])
    for (j in 1:gn)
      if (count[j] > 0)
        points(i, j, pch=22, cex=sc * sqrt(count[j]), lwd=sc3 * count[j])
  }
}
mtext("Question", side=1, line=0.2, outer=TRUE)   
mtext("Answer", side=2, line=0.2, outer=TRUE)   

