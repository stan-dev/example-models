# clears workspace: 
rm(list=ls()) 
library(rstan)

model <- "
// Knower Level Model Applied to Fast-Cards Data
data { 
  int<lower=1> ns;
  int<lower=1> nz;
  int<lower=1> gn;
  int<lower=1> fn;
  int fnq[ns];
  int fq[ns,21];  // no. of columns = max(fnq)
  int fa[ns,21];  // no. of columns = max(fnq)
}
parameters {
  vector<lower=0,upper=1>[fn] fpitmp;
  real<lower=1,upper=1000> v;
} 
transformed parameters {
  simplex[fn] fpi;
  simplex[fn] fnpiprime[nz,gn];
  vector[nz] lp_parts[ns];

  // Base rate
  fpi <- fpitmp / sum(fpitmp);
  
  // Model
  for (i in 1:nz) {
    for (j in 1:gn) {
      vector[fn] fpiprime;
      for (k in 1:fn) {
        real find1;
        real find2;
        real find3;
        real find4;
        real find5;
        
        // Will be 1 if Knower-Level (i.e, i-1) is Same or Greater than Answer
        find1 <- step((i - 1) - k);
        // Will be 1 for the Possible Answer that Matches the Question
        find2 <- k == j;
        // Will be 1 for 0-Knowers
        find3 <- i == 1;
        // Will be 1 for HN-Knowers
        find4 <- i == nz;
        find5 <- find3 + find4 * (2 + find2) 
                + (1 - find4) * (1 - find3) * (find1 * find2 + find1 + 1);
                
        if (find5 == 1)
          fpiprime[k] <- fpi[k];
        else if (find5 == 2)  
          fpiprime[k] <- 1 / v * fpi[k];
        else if (find5 == 3)  
          fpiprime[k] <- v * fpi[k];
      }  
      for (k in 1:fn)
        fnpiprime[i,j,k] <- fpiprime[k] / sum(fpiprime);
    }
  }
  
  for (i in 1:ns) {
    for (m in 1:nz) {
      real lp_parts_tmp;
      lp_parts_tmp <- 0;
      
      // Probability a z[i]-Knower Will Give fa[i,j] when Asked for fq[i,j]
      for (j in 1:fnq[i])
        lp_parts_tmp <- lp_parts_tmp
                         + categorical_log(fa[i,j], fnpiprime[m,fq[i,j]]);

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
  int predfa[ns,gn];
  int predz[nz,gn];
  int predfpi;
  
  for (i in 1:ns) {
    prob <- softmax(lp_parts[i]);
    z[i] <- categorical_rng(prob);
  }
  
  // Posterior Predictive
  for (i in 1:ns)
    for (j in 1:gn)
      predfa[i,j] <- categorical_rng(fnpiprime[z[i],j]);
  
  // Posterior Prediction For Knower Levels
  for (i in 1:nz)
    for (j in 1:gn)
      predz[i,j] <- categorical_rng(fnpiprime[i,j]);
  
  predfpi <- categorical_rng(fpi);
}"

load("fc_given.RData")  # Load all data for the model

# Answers in data
wlist <- sort(setdiff(unique(as.vector(fa)), 0))
fn <- length(wlist)
fqm <- matrix(0, nrow(fa), ncol(fa))
fam <- matrix(0, nrow(fa), ncol(fa))

for (i in 1:ns) {
  for (j in 1:max(fnq)) {
    if (fa[i, j] == 0)
      fam[i, j] <- 0
    else 
      fam[i, j] <- which(fa[i, j] == wlist)
    
    if (fq[i, j] == 0)
      fqm[i, j] <- 0
    else 
      fqm[i, j] <- which(fq[i, j] == wlist)
  }
}
fa <- fam
fq <- fqm

# to be passed on to Stan
data <- list(ns=ns, fnq=fnq, fn=fn, fa=fa, fq=fq, nz=nz, gn=gn) 

myinits <- list(
  list(v=2, fpitmp=rep(1 / fn, fn)),
  list(v=2, fpitmp=rep(1 / fn, fn)))

# parameters to be monitored:  
parameters <- c("predfa", "predz", "predfpi", "v", "z")

# The following command calls Stan with specific options.
# For a detailed description type "?stan".
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
predfpi <- extract(samples)$predfpi
predz <- extract(samples)$predz

#### Figure 19.7 ####
tab <- table(factor(0, levels=1:50))
tabtmp <- table(predfpi)
names(tabtmp) <- wlist

for (i in names(tabtmp))
  tab[i] <- tabtmp[i]

windows(8,5)
par(mar=c(3, 2, 1, 1) + .1, mgp=c(1.3, 0.2, 0), cex.lab=1.2)
barplot(as.vector(tab), ylim=c(0, max(tab) * 1.2), col="black", yaxt="n", 
        xlab="Number", ylab="")
title(ylab="Probability", line=.2)
axis(1, at=c(.7, 11.5, 23.5, 35.5, 47.5, 59.5), tck=0,
     labels=c("1", "10", "20", "30", "40", "50"))
box()

#### Figure 19.8 ####
windows(9, 6)
par(mfrow=c(4, 5), mar=c(1, 0, 2, 2) + .1, oma=c(2.6, 2.8, 1, 0),
    mgp=c(1.5, 0.15, 0))
for (i in 1:ns) {
  zTable <- table(factor(zSamples[, i], levels=as.character(1:6)))
  
  barplot(zTable, col="black", xaxt="n", yaxt="n", main=paste("Child", i),
          ylim=c(0, length(zSamples[, 1]))) 
  
  if (i == 16)
    axis(1, at=seq(.7, 7.5, by=1.2), label=c("P", "1", "2", "3", "4", "C"),
         tck = 0.1)
  else
    axis(1, at=seq(.7, 7.5, by=1.2), label=FALSE, tck = 0.1)
  
  axis(3, at=seq(.7, 7.5, by=1.2), label=FALSE, tck = 0.1)
  box()
}
mtext("Knower", side=1, line=0.2, at=.055, adj=0, outer=TRUE)  
mtext("Level", side=1, line=1.2, at=.065, adj=0, outer=TRUE)  
mtext("Posterior", side=2, line=1.2, at=.07, adj=0, outer=TRUE)  
mtext("Mass", side=2, line=0.2, at=.09, adj=0, outer=TRUE)

#### Figure 19.9 ####
sc=2
sc3=2
cm=-.4
cb=.99
shadedSize <- 2

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
  plot(NA, xlim=c(.5, 15.5), ylim=c(.5, 20.5), axes=FALSE, cex.main=1, ylab="", 
       main=mainTitle[z], xlab="")
  axis(1, at=c(1, 2, 3, 4, 5, 8, 10), tck=0)
  axis(2, at=c(1, 10, 20), las=1, tck=0)
  axis(3, at=c(1, 2, 3, 4, 5, 8, 10), labels=FALSE, tck=0)
  axis(4, at=1:15, labels=FALSE, tck=0)
  box()
  
  for (i in 1:gn) {
    count <- hist(predz[, z, i], plot=FALSE, breaks=seq(0.5, 15.5))$counts
    count <- count / max(count)
    
    for (j in 1:fn)
      points(i, wlist[j], pch=15, col=gray(min(1, max(0, (cm * count[j] + cb)))), 
             cex=shadedSize)
  }
  match <- which(mv == z)
  if (length(match) != 0)
    for (i in 1:length(match))
      for (j in 1:fnq[match[i]])
        dm[fqm[match[i], j], fam[match[i], j], z] <- dm[fqm[match[i], j],
                                                        fam[match[i], j], z] + 1
    
  for (i in 1:gn) {
    if (sum(dm[i, , z]) == 0)
      count <- rep(0, 15)
    else
      count <- dm[i, , z] / sum(dm[i, , z])
    for (j in 1:fn)
      if (count[j] > 0)
        points(i, wlist[j], pch=22, cex=sc * sqrt(count[j]), lwd=sc3 * count[j])
  }
}
mtext("Question", side=1, line=0.2, outer=TRUE)   
mtext("Answer", side=2, line=0.2, outer=TRUE)   

