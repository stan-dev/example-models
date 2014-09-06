### WARNING: Run PsychometricFunction1.R first and do not close the R window. 

# Set working directory!

library(rstan)

model <- "
// Logistic Psychophysical Function with Contaminants
data { 
  int nsubjs;
  int nstim[nsubjs];
  int n[nsubjs,28];
  int r[nsubjs,28];
  int x[nsubjs,28];
  vector[nsubjs] xmean; 
}
parameters {
  real mua;
  real mub;
  real mup;
  real<lower=0,upper=1000> sigmaa;
  real<lower=0,upper=1000> sigmab;
  real<lower=0,upper=3> sigmap;
  vector[nsubjs] raw[3];
  matrix<lower=0,upper=1>[nsubjs,28] pi;    
} 
transformed parameters {
  vector[2] lp_parts[nsubjs,28];
  vector<lower=0,upper=1>[nsubjs] phi;  
  vector[nsubjs] alpha;
  vector[nsubjs] beta;
  vector[nsubjs] probitphi;
  
  alpha <- mua + sigmaa * raw[1];  // Matt Trick
  beta <- mub + sigmab * raw[2];  // Matt Trick
  probitphi <- mup + sigmap * raw[3];  // Matt Trick

  for (i in 1:nsubjs)
    phi[i] <- Phi_approx(probitphi[i]);  // just approx. to improve performance
  
  for (i in 1:nsubjs) {
    for (j in 1:nstim[i]) {  
      real theta; 
      theta <- inv_logit(alpha[i] + beta[i] * (x[i,j] - xmean[i]));
      
      lp_parts[i,j,1] <- log1m(phi[i]) + binomial_log(r[i,j], n[i,j], theta);
      lp_parts[i,j,2] <- log(phi[i]) + binomial_log(r[i,j], n[i,j], pi[i,j]);
    }
  }
}
model {
  // Priors
  mua ~ normal(0, sqrt(1000));
  mub ~ normal(0, sqrt(1000));
  mup ~ normal(0, 1); 
  
  for (i in 1:3)
    raw[i] ~ normal(0, 1);  // for Matt Trick

  for (i in 1:nsubjs)
    for (j in 1:nstim[i])
      increment_log_prob(log_sum_exp(lp_parts[i,j]));
}
generated quantities {
  int<lower=0,upper=1> z[nsubjs,28];

  for (i in 1:nsubjs) {
    for (j in 1:nstim[i]) {  
      vector[2] prob;
      prob <- softmax(lp_parts[i,j]);
      z[i,j] <- bernoulli_rng(prob[2]);
    }
  }
}"

x <- as.matrix(read.table("data_x.txt", sep="\t"))
x[is.na(x)] = -99  # transforming because Stan won't accept NAs 

n <- as.matrix(read.table("data_n.txt", sep="\t"))
n[is.na(n)] = -99  # transforming because Stan won't accept NAs 

r <- as.matrix(read.table("data_r.txt", sep="\t"))
r[is.na(r)] = -99  # transforming because Stan won't accept NAs 

rprop <- as.matrix(read.table("data_rprop.txt", sep="\t"))

xmean <- c(318.888, 311.0417, 284.4444, 301.5909, 
           296.2000, 305.7692, 294.6429, 280.3571)
nstim <- c(27, 24, 27, 22, 25, 26, 28, 28)
nsubjs <- 8

# to be passed on to Stan
data <- list(x=x, xmean=xmean, n=n, r=r, nsubjs=nsubjs, nstim=nstim) 

myinits <- list( 
  list(raw=matrix(rep(0, 3 * nsubjs), 3, nsubjs),
       mua=0, mub=0, mup=0, sigmaa=1, sigmab=1, sigmap=1,  
       pi=matrix(runif(nsubjs * 28), nsubjs, 28)),
  list(raw=matrix(rep(0, 3 * nsubjs), 3, nsubjs),
       mua=0, mub=0, mup=0, sigmaa=1, sigmab=1, sigmap=1,  
       pi=matrix(runif(nsubjs * 28), nsubjs, 28)))

parameters <- c("alpha", "beta", "z")  # Parameters to be monitored

# The following command calls Stan with specific options.
# For a detailed description type "?rstan".
samples <- stan(model_code=model,
                data=data, 
                init=myinits,
                pars=parameters,
                iter=1200, 
                chains=2, 
                thin=1,
                warmup=200,  # Stands for burn-in; Default = iter/2
                # seed = 123  # Setting seed; Default is random seed
)
# Now the values for the monitored parameters are in the "samples" object, 
# ready for inspection.

# Transforming back to NAs
x[x == -99] = NA 
n[n == -99] = NA
r[r == -99] = NA

# Extracting the necessary parameters
zmean = list()
for (i in 1:8){
  ztemp = matrix(NA,nstim[i],1)
  for (j in 1:nstim[i]){
    ztemp[j,] = summary(samples)$summary[paste("z[",as.character(i),",",as.character(j),"]",sep=""),"mean"]
  }
  zmean[[i]] = ztemp
}
alpha2    = extract(samples)$alpha
beta2     = extract(samples)$beta
alphaMAP2  = c(rep(0,nsubjs))
betaMAP2  = c(rep(0,nsubjs))
alpha_sel2 = matrix(NA,20,8) 
beta_sel2 = matrix(NA,20,8) 

# Constructing MAP-estimates and alpha/beta range
for (i in 1:nsubjs)
{
	alphaMAP2[i]   <- density(alpha2[,i])$x[which(density(alpha2[,i])$y==max(density(alpha2[,i])$y))]
	betaMAP2[i]    <- density(beta2[,i])$x[which(density(beta2[,i])$y==max(density(beta2[,i])$y))]
	alpha_sel2[,i] <- sample(alpha2[,i],20)
	beta_sel2[,i]  <- sample(beta2[,i],20)
}
	
############################## PSYCHOMETRIC FUNCTIONS ##############################

F4 <- function(X,s) # only the MAP estimate; use this to plot psychometric functions
{
  exp(alphaMAP2[s] + betaMAP2[s]*(X - xmean[s]))/(1+exp(alphaMAP2[s] + betaMAP2[s]*(X - xmean[s])))
}

F4inv <- function(Y,s)
{
  (log(-Y/(Y-1))-alphaMAP2[s])/betaMAP2[s]
}

F5 <- function(X,s) # function for all the posterior alpha/beta values; use this to calculate JND posterior
{
  exp(alpha2[,s] + beta2[,s]*(X - xmean[s]))/(1+exp(alpha2[,s] + beta2[,s]*(X - xmean[s])))
}

F5inv <- function(Y,s)
{
  (log(-Y/(Y-1))-alpha2[,s])/beta2[,s]
}

F6 <- function(X,s,g) # function for 20 grabbed posterior alpha/beta values; use this to plot overlapping sigmoids to visualize variance
{
  exp(alpha_sel2[g,s] + beta_sel2[g,s]*(X - xmean[s]))/(1+exp(alpha_sel2[g,s] + beta_sel2[g,s]*(X - xmean[s])))
}

##################################### JND/PSE calculation ########################################

	JND2 	  <- F5inv(0.84,c(1:nsubjs))-F5inv(0.5,c(1:nsubjs))
	JNDmap2 <- F4inv(0.84,c(1:nsubjs))-F4inv(0.5,c(1:nsubjs))								  				             
	PSE2 	  <- F5inv(0.5,c(1:nsubjs))+xmean
	PSEmap2 <- F4inv(0.5,c(1:nsubjs))+xmean
		
################## PLOTS ####################

### Figure 12.6 
dev.new(width=10,height=5)
layout(matrix(1:nsubjs,2,4,byrow=T))
par(mar=c(1,2,2,0),oma=c(5,5,1,1))
for (i in 1:nsubjs)
{
	scale <- seq(x[i,1],x[i,nstim[i]], by=.1)
	plot(x[i,],rprop[i,],main=paste("Subject",as.character(i)),xlab="",ylab="",pch=22, col="black", bg=c(grey(1-zmean[[i]])), ylim=c(0,1), yaxt="n",xaxt="n")
	lines(scale,F1(scale,i),type="l", col="light grey",lwd=3)
	lines(scale,F4(scale,i),type="l")
	segments(x0=x[i,1],x1=PSEmap2[i]+JNDmap2[i],y0=0.84,lty=2)
	segments(x0=x[i,1],x1=PSEmap2[i],y0=0.5,lty=2)
	segments(y0=0,y1=0.84,x0=PSEmap2[i]+JNDmap2[i],lty=2)
	segments(y0=0,y1=0.5,x0=PSEmap2[i],lty=2)
	if (i==1 | i==5) 
  {
		axis(2,las=1,yaxp=c(0,1,2))
		axis(2,at=0.84,las=1)
	}
	if (i>4) 
    axis(1)
}
mtext("Proportion 'Long' Response",side=2,line=2,outer=T,cex=1.4)
mtext("Test Interval (ms)",side=1,outer=T,line=3,cex=1.4)

### Figure 12.7
	
dev.new(width=10,height=5)
layout(matrix(1:nsubjs,2,4,byrow=T))
par(mar=c(1,2,2,0),oma=c(5,5,1,1))
for (i in 1:nsubjs)
{
	plot(density(JND[,i]), col="dark grey", type="h", main=paste("Subject",as.character(i)), ylab="", xlab="", xlim=c(10,100), ylim=c(0,0.12), las=1,yaxt="n",xaxt="n")
	par(new=T)
	plot(density(JND2[,i]), col="light grey",type="h", main="",ylab="",xlab="",xlim=c(10,100),ylim=c(0,0.12),yaxt="n",xaxt="n",bty="n")
	par(new=F)
	if (i==1 | i==5) axis(2,las=1, yaxp=c(0,0.12,3))
	if (i>4) axis(1)
	if (i==4) legend("topright",c("non-contaminant model","contaminant model"),col=c("dark grey","light grey"),bty="n",pch=c(15,15))
}
mtext("Posterior Density",side=2,line=2,outer=T,cex=1.4)
mtext("JND (ms)",side=1,line=3,outer=T,cex=1.4)