# clears workspace: 
rm(list=ls(all=TRUE)) 

library(rstan)

model <- '
# SIMPLE Model
data { 
    int dsets;
    int y[40, dsets];
    int n[dsets];
    int listlength[dsets];
    int m[40, dsets];
}
transformed data {
	int maxlistlength;
	
	maxlistlength <- max(listlength);
}
parameters {
    vector<lower=0, upper=100>[dsets] c;
    vector<lower=0, upper=100>[dsets] s;
    vector<lower=0, upper=1>[dsets] t;
} 
transformed parameters {
    matrix[maxlistlength, dsets] theta;
    vector<lower=0, upper=1>[maxlistlength] sim;
    vector<lower=0, upper=1>[maxlistlength] disc;
    vector<lower=0, upper=1>[maxlistlength] resp;
    
    # Similarities, Discriminabilities, and Response Probabilities
    for (x in 1:dsets) {
        for (i in 1:listlength[x]) {
            for (j in 1:maxlistlength) {
                if (j <= listlength[x]) {
                    # Similarities
                    sim[j] <- exp(-c[x] * fabs(log(m[i, x]) - log(m[j, x])));
                } else {				
                    sim[j] <- 0.0;
                    resp[j] <- 0.0;
                }
            }
            # Discriminabilities
            disc <- sim / sum(sim);

            for (j in 1:listlength[x]) { 
                # Response Probabilities
                resp[j] <- 1.0 / (1.0 + exp(-s[x] * (disc[j] - t[x])));
            }

            # Free Recall Overall Response Probability
            theta[i, x] <- fmin(.999, sum(resp));
        }
    }
}
model {
    # Prior
    t ~ beta(1, 1);
    
    # Observed Data
    for (x in 1:dsets) {
        for (i in 1:listlength[x]) {
            y[i, x] ~ binomial(n[x], theta[i, x]);
        }
    }
}
generated quantities {
    real predy[maxlistlength, dsets];
    real predpc[maxlistlength, dsets];
    
    # Predicted Data
    for (x in 1:dsets) {
        for (i in 1:listlength[x]) {
            predy[i, x] <- binomial_rng(n[x], theta[i, x]);
            predpc[i, x] <- predy[i, x] / n[x];
        }
    } 	
}'

y          <- matrix(scan("k_M.txt", sep=","), ncol=40, nrow=6, byrow=T) 
n          <- c(1440,1280,1520,1520,1200,1280)
listlength <- c(10,15,20,20,30,40)
pc         <- matrix(scan("pc_M.txt", sep=","), ncol=40, nrow=6, byrow=T) 

dsets <- 6 

# Loop over conditions
m <- y*0
for (dset in 1:dsets)
{
    if (dset==1)
    {
        nwords <- 10
        lag    <- 2
        offset <- 15
    } 
    if (dset==2)
    {
        nwords <- 15
        lag    <- 2
        offset <- 20
    } 
    if (dset==3)
    {
        nwords <- 20
        lag    <- 2
        offset <- 25
    } 
    if (dset==4)
    {
        nwords <- 20
        lag    <- 1
        offset <- 10
    } 
    if (dset==5)
    {
        nwords <- 30
        lag    <- 1
        offset <- 15
    } 
    if (dset==6)
    {
        nwords <- 40
        lag    <- 1
        offset <- 20
    } 
    # Temporal Offset For Free Recall
    m[dset,1:nwords] <- offset+seq(from=(nwords-1)*lag, by=-lag, to=0)
}

y <- t(y)
m <- t(m)

# To be passed on to Stan
data <- list(y=y, n=n, listlength=listlength, m=m, dsets=dsets) 

myinits <- list(
    list(c=rep(15, dsets), s=rep(10, dsets), t=rep(.5, dsets)),
    list(c=rep(15, dsets), s=rep(10, dsets), t=rep(.5, dsets)),
    list(c=rep(15, dsets), s=rep(10, dsets), t=rep(.5, dsets)))

parameters <- c("c", "t", "s", "predpc")  # Parameters to be monitored

# The following command calls Stan with specific options.
# For a detailed description type "?rstan".
samples <- stan(model_code=model,   
                data=data, 
                init=myinits,  # If not specified, gives random inits
                pars=parameters,
                iter=1000, 
                chains=3, 
                thin=1,
                warmup = 200,  # Stands for burn-in; Don't set to 0 or 
                # ... low values, it can malfunction; Default = iter/2
                # seed = 123  # Setting seed; Default is random seed
)
# Now the values for the monitored parameters are in the "samples" object, 
# ready for inspection.

# NB. The predpc values of -999 are obviously dummy values;
# they should not appear in any of the plots.
# NB2. For better convergence, you might want to take a lunch break and run
# multiple chains, for more iterations, and do some thinning.

#Figure 15.2
layout(matrix(c(
    7,1,2,3,
    7,4,5,6,
    8,8,8,8
),3,4,byrow=T), c(1,2,2,2), c(2,2,.5))
layout.show(8)
hm <- 20
ll <- listlength

for (dset in 1:dsets) {
    plot(-1,-1,xlim=c(0,40),ylim=c(0,1),xlab="",ylab="",las=1)
    for (i in 1:ll[dset]) { 
        data <- extract(samples)$predpc[, i, dset]
        points(i+runif(hm,0,1)*.1,data[ceiling(runif(hm,0,1)*length(data))],col="grey")
    }
    points(1:ll[dset],pc[dset,1:ll[dset]],xlim=c(0,40),ylim=c(0,1))
    lines(1:ll[dset],pc[dset,1:ll[dset]])
    
    box("plot")
}
par(mar=c(rep(0,4)))
plot.new()
text(.45,.5,"Probability Correct",cex=2.5,srt=90)
plot.new()
text(.5,.5,"Serial Position",cex=2.5,mar=c(rep(0,4)))
