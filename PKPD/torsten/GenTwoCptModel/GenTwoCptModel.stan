## GenTwoCptModelExample.stan
## v1.3 - Example
## Test model for Stan Pmetrics function PKModelTwoCpt.
## Heavily anotated to help new users

functions{
  
  # define ODE system for two compartmnt model
  real[] twoCptModelODE(real t,
			real[] x,
			real[] parms,
			real[] rate, # in this example, rate is treated as data
			int[] dummy){
			  
    # Parameters
    real Q;
    real CL;
    real V1;
    real V2;
    real ka;
    
    # Re-parametrization
    real k12;
    real k21;
    real k10;
    
    # Return object (derivative)
    real y[3]; # 1 element per compartment of
               # the model

    CL = parms[1];
    Q = parms[2];
    V1 = parms[3];
    V2 = parms[4];
    ka = parms[5];
   
    # Reparametrization (for PK model)
    k10 = CL / V1;
    k12 = Q / V1;
    k21 = Q / V2;

    # PK component of the ODE system
    y[1] = -ka*x[1];
    y[2] = ka*x[1] - (k10 + k12)*x[2] + k21*x[3];
    y[3] = k12*x[2] - k21*x[3];

    return y;
  }
  
}
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
  int nTheta;
  
  logCObs = log(cObs);
  nTheta = 11; # number of parameters 
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

  # generalCptModel takes in the ODE system, the number of compartment 
  # (here we have a two compartment model with first order absorption, so
  # three compartments), the parameters matrix, the NONEM data, and the tuning
  # parameters (relative tolerance, absolute tolerance, and maximum number of steps)
  # of the ODE integrator. The user can choose between the bdf and the rk45 integrator. 
  # Returns a matrix with the predicted amount in each compartment 
  # at each event.

//  x = generalCptModel_bdf(twoCptModelODE, 3,
//                          theta, time, amt, rate, ii, evid, cmt, addl, ss,
//                          1e-8, 1e-8, 1e8);

   x = generalCptModel_rk45(twoCptModelODE, 3,
                           theta, time, amt, rate, ii, evid, cmt, addl, ss,
                           1e-8, 1e-8, 1e8);

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


