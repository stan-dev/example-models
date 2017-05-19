## TwoCptModel.stan
## Run two compartment model using built-in analytical solution 
## Heavily anotated to help new users

data{
  int<lower = 1> nt; # number of events
  int<lower = 1> nObs; # number of observation
  int<lower = 1> iObs[nObs]; # index of observation
  
  # NONMEM data
  int<lower = 1> cmt[nt];
  int evid[nt];
  int addl[nt];
  int ss[nt];
  real amt[nt];
  real time[nt];
  real rate[nt];
  real ii[nt];
  
  vector<lower = 0>[nObs] cObs; # observed concentration (Dependent Variable)
}

transformed data{
  vector[nObs] logCObs;
  int nTheta;
  
  logCObs = log(cObs);
  nTheta = 11; # number of parameters in Two Compartment Model
}

parameters{
  real<lower = 0> CL;
  real<lower = 0> Q;
  real<lower = 0> V1;
  real<lower = 0> V2;
  real<lower = 0> ka;
  real<lower = 0> sigma;

}

transformed parameters{
  vector<lower = 0>[nTheta] theta[1]; # The [1] indicates the parameters are constant
  vector<lower = 0>[nt] cHat;
  vector<lower = 0>[nObs] cHatObs;
  matrix<lower = 0>[nt, 3] x; 

  theta[1][1] = CL;
  theta[1][2] = Q;
  theta[1][3] = V1;
  theta[1][4] = V2;
  theta[1][5] = ka;
  theta[1][6] = 1; # F1
  theta[1][7] = 1; # F2
  theta[1][8] = 1; # F3
  theta[1][9] = 0; # tlag1
  theta[1][10] = 0; # tlag2
  theta[1][11] = 0; # tlag3

  # PKModelTwoCpt takes in the parameters matrix and the NONMEM data,
  # and returns a matrix with the predicted amount in each compartment 
  # at each event.
  x = PKModelTwoCpt(theta, time, amt, rate, ii, evid, cmt, addl, ss);

  cHat = col(x, 2) ./ V1; # we're interested in the amount in the second compartment 

  for(i in 1:nObs){
    cHatObs[i] = cHat[iObs[i]]; ## predictions for observed data records
  }
}

model{
  # informative prior
  CL ~ lognormal(log(10), 0.25);
  Q ~ lognormal(log(15), 0.5);
  V1 ~ lognormal(log(35), 0.25);
  V2 ~ lognormal(log(105), 0.5);
  ka ~ lognormal(log(2.5), 1);
  sigma ~ cauchy(0, 1);

  logCObs ~ normal(log(cHatObs), sigma);
}

generated quantities{
  real cObsPred[nObs];

  for(i in 1:nObs){
      cObsPred[i] = exp(normal_rng(log(cHatObs[i]), sigma));
    }
			 
}


