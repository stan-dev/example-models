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
  int gnq[ns];
  int fnq[ns];
  int gq[ns,21];  // no. of columns = max(gnq)
  int ga[ns,21];  // no. of columns = max(gnq)
  int fq[ns,21];  // no. of columns = max(fnq)
  int fa[ns,21];  // no. of columns = max(fnq)
}
parameters {
  vector<lower=0,upper=1>[gn] pitmp;
  vector<lower=0,upper=1>[fn] fpitmp;
  real<lower=1,upper=1000> gv;
  real<lower=1,upper=1000> fv;
} 
transformed parameters {
  simplex[gn] pi;
  simplex[fn] fpi;
  simplex[gn] npiprime[nz,gn];
  simplex[fn] fnpiprime[nz,gn];
  vector[nz] lp_parts[ns];

  // Base rates
  pi <- pitmp / sum(pitmp);
  fpi <- fpitmp / sum(fpitmp);
  
  // Give-N Part Model
  for (i in 1:nz) {
    for (j in 1:gn) {
      vector[gn] piprime;
      for (k in 1:gn) {
        real ind1;
        real ind2;
        real ind3;
        real ind4;
        real ind5;
        
        ind1 <- step((i - 1) - k);
        ind2 <- k == j;
        ind3 <- i == 1;
        ind4 <- i == nz;
        ind5 <- ind3 + ind4 * (2 + ind2) 
                + (1 - ind4) * (1 - ind3) * (ind1 * ind2 + ind1 + 1);
                
        if (ind5 == 1)
          piprime[k] <- pi[k];
        else if (ind5 == 2)  
          piprime[k] <- 1 / gv * pi[k];
        else if (ind5 == 3)  
          piprime[k] <- gv * pi[k];
      }  
      for (k in 1:gn)
        npiprime[i,j,k] <- piprime[k] / sum(piprime);
    }
  }
  
  // Fast-Cards Model
  for (i in 1:nz) {
    for (j in 1:gn) {
      vector[fn] fpiprime;
      for (k in 1:fn) {
        real find1;
        real find2;
        real find3;
        real find4;
        real find5;
        
        find1 <- step((i - 1) - k);
        find2 <- k == j;
        find3 <- i == 1;
        find4 <- i == nz;
        find5 <- find3 + find4 * (2 + find2) 
                + (1 - find4) * (1 - find3) * (find1 * find2 + find1 + 1);
                
        if (find5 == 1)
          fpiprime[k] <- fpi[k];
        else if (find5 == 2)  
          fpiprime[k] <- 1 / fv * fpi[k];
        else if (find5 == 3)  
          fpiprime[k] <- fv * fpi[k];
      }  
      for (k in 1:fn)
        fnpiprime[i,j,k] <- fpiprime[k] / sum(fpiprime);
    }
  }
  
  for (i in 1:ns) {
    for (m in 1:nz) {
      real lp_parts_tmp;
      lp_parts_tmp <- 0;
      
      for (j in 1:gnq[i])
        lp_parts_tmp <- lp_parts_tmp
                         + categorical_log(ga[i,j], npiprime[m,gq[i,j]]);
      
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
  int predga[ns,gn];
  int predfa[ns,gn];
  int predgaz[nz,gn];
  int predfaz[nz,gn];
  int predpi;
  int predfpi;
  
  for (i in 1:ns) {
    prob <- softmax(lp_parts[i]);
    z[i] <- categorical_rng(prob);
  }
  
  // Posterior Predictive
  for (i in 1:ns) {
    for (j in 1:gn) {
      predga[i,j] <- categorical_rng(npiprime[z[i],j]);
      predfa[i,j] <- categorical_rng(fnpiprime[z[i],j]);
    }
  }
  
  // Posterior Prediction For Knower Levels
  for (i in 1:nz) {
    for (j in 1:gn) {
      predgaz[i,j] <- categorical_rng(npiprime[i,j]);
      predfaz[i,j] <- categorical_rng(fnpiprime[i,j]);
    }
  }
  predpi <- categorical_rng(pi);
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
data <- list(ns=ns, fnq=fnq, fn=fn, fa=fa, fq=fq, nz=nz, gnq=gnq, gn=gn, 
             ga=ga, gq=gq) 

myinits <- list(
  list(gv=5, fv=5, pitmp=rep(1 / gn, gn), fpitmp=rep(1 / fn, fn)),
  list(gv=5, fv=5, pitmp=rep(1 / gn, gn), fpitmp=rep(1 / fn, fn)))

# parameters to be monitored:  
parameters <- c("predga", "predfa", "predfaz", "predgaz", "predpi", "predfpi", 
                "fv", "gv", "z")

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
predpi <- extract(samples)$predpi

#### a Figure ####
tab <- table(factor(0, levels=1:50))
tabtmp <- table(predfpi)
names(tabtmp) <- wlist
for (i in names(tabtmp))
  tab[i] <- tabtmp[i]

windows(8,4)
par(mfrow=c(1, 2), mar=c(3, 2, 1, 1) + .1, mgp=c(1.3, 0.2, 0), cex.lab=1.2)

par(mar=c(3, 2, 1, 1) + .1, mgp=c(1.3, 0.2, 0), cex.lab=1.2)
barplot(table(predpi), ylim=c(0, max(table(predpi) * 1.2)), col="black", yaxt="n",
        xlab="Number", ylab="")
box()
title(ylab="Probability", line=.2)
axis(3, at=seq(.7, 17.5, by=1.2), label=FALSE, tck = 0.02)

barplot(as.vector(tab), ylim=c(0, max(tab) * 1.2), col="black", yaxt="n",
        xlab="Number", ylab="")
title(ylab="Probability", line=.2)
axis(1, at=c(.7, 11.5, 23.5, 35.5, 47.5, 59.5), tck=0,
     labels=c("1", "10", "20", "30", "40", "50"))
box()

#### Figure 19.11 ####
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
