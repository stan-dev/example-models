## LinTwoCptModelExample.stan
## Run two compartment model using matrix exponential solution
## Heavily anotated to help new users

data{
  int<lower = 1> nt; # number of events
  int<lower = 1> nObs; # number of observations
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
  
  vector<lower = 0>[nObs] cObs; # observed concentration (dependent variable)
}

transformed data{
  vector[nObs] logCObs;

  logCObs = log(cObs);
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
  matrix[3, 3] K;
  real k10;
  real k12;
  real k21;
  vector<lower = 0>[6] theta[1]; # The [1] indicates the parameters are constant
  vector<lower = 0>[nt] cHat;
  vector<lower = 0>[nObs] cHatObs;
  matrix<lower = 0>[nt, 3] x;

  k10 = CL / V1;
  k12 = Q / V1;
  k21 = Q / V2;
  
  K = rep_matrix(0, 3, 3);
  
  K[1, 1] = -ka;
  K[2, 1] = ka;
  K[2, 2] = -(k10 + k12);
  K[2, 3] = k21;
  K[3, 2] = k12;
  K[3, 3] = -k21;

  theta[1][1] = 1; # F1
  theta[1][2] = 1; # F2
  theta[1][3] = 1; # F3
  theta[1][4] = 0; # tlag1
  theta[1][5] = 0; # tlag2
  theta[1][6] = 0; # tlag3

  # linCptModel takes in the constant rate matrix, the object theta which
  # contains the biovariability fraction and the lag time of each compartment,
  # and the NONMEM data.
  x = linCptModel(K, theta, time, amt, rate, ii, evid, cmt, addl, ss);

  cHat = col(x, 2) ./ V1;

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
