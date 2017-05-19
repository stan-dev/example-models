functions{

    real[] twoCptNeutModelODE(real t,
			real[] x,
			real[] parms,
			real[] rdummy,
			int[] idummy){
    real k10;
    real k12;
    real k21;
    real CL;
    real Q;
    real V1;
    real V2;
    real ka;
    real mtt;
    real circ0;
    real gamma;
    real alpha;
    real ktr;
    real dxdt[8];
    real conc;
    real EDrug;
    real transit1;
    real transit2;
    real transit3;
    real circ;
    real prol;

    CL = parms[1];
    Q = parms[2];
    V1 = parms[3];
    V2 = parms[4];
    ka = parms[5];
    mtt = parms[6];	
    circ0 = parms[7];
    gamma = parms[8];
    alpha = parms[9];

    k10 = CL / V1;
    k12 = Q / V1;
    k21 = Q / V2;

    ktr = 4 / mtt;
  
    dxdt[1] = -ka * x[1];
    dxdt[2] = ka * x[1] - (k10 + k12) * x[2] + k21 * x[3];
    dxdt[3] = k12 * x[2] - k21 * x[3];
    conc = x[2]/V1;
    EDrug = alpha * conc;
    // x[4], x[5], x[6], x[7] and x[8] are differences from circ0.
    prol = x[4] + circ0;
    transit1 = x[5] + circ0;
    transit2 = x[6] + circ0;
    transit3 = x[7] + circ0;
    circ = fmax(machine_precision(), x[8] + circ0); // Device for implementing a modeled 
                                                    // initial condition
    dxdt[4] = ktr * prol * ((1 - EDrug) * ((circ0 / circ)^gamma) - 1);
    dxdt[5] = ktr * (prol - transit1);
    dxdt[6] = ktr * (transit1 - transit2);
    dxdt[7] = ktr * (transit2 - transit3);
    dxdt[8] = ktr * (transit3 - circ);

    return dxdt;

  }
}

data{
  int<lower = 1> nt;
  int<lower = 1> nObsPK;
  int<lower = 1> nObsPD;
  int<lower = 1> iObsPK[nObsPK];
  int<lower = 1> iObsPD[nObsPD];
  real<lower = 0> amt[nt];
  int<lower = 1> cmt[nt];
  int<lower = 0> evid[nt];
  real<lower = 0> time[nt];
  real<lower = 0> ii[nt];
  int<lower = 0> addl[nt];
  int<lower = 0> ss[nt];
  real rate[nt];
  vector<lower = 0>[nObsPK] cObs;
  vector<lower = 0>[nObsPD] neutObs;
  
  ## data for population model
  int<lower = 1> nSubjects;
  int<lower = 1> start[nSubjects];
  int<lower = 1> end[nSubjects];
  real<lower = 0> weight[nSubjects];
  
  ## data for priors
  // real<lower = 0> CLHatPrior;
  // real<lower = 0> QHatPrior;
  // real<lower = 0> V1HatPrior;
  // real<lower = 0> V2HatPrior;
  // real<lower = 0> kaHatPrior;
  // real<lower = 0> CLHatPriorCV;
  // real<lower = 0> QHatPriorCV;
  // real<lower = 0> V1HatPriorCV;
  // real<lower = 0> V2HatPriorCV;
  // real<lower = 0> kaHatPriorCV;
  real<lower = 0> circ0HatPrior;
  real<lower = 0> circ0HatPriorCV;
  real<lower = 0> mttHatPrior;
  real<lower = 0> mttHatPriorCV;
  real<lower = 0> gammaPrior;
  real<lower = 0> gammaPriorCV;
  real<lower = 0> alphaHatPrior;
  real<lower = 0> alphaHatPriorCV;
}

transformed data{
  vector[nObsPK] logCObs;
  vector[nObsPD] logNeutObs;
##  int idummy[0];
##  real rdummy[0];

  int nTheta;
  int nIIV;

  logCObs = log(cObs);
  logNeutObs = log(neutObs);
  
  nIIV = 7; ## parameters with IIV
  nTheta = 25; ## number of parameters

}

parameters{

  real<lower = 0> CLHat;
  real<lower = 0> QHat;
  real<lower = 0> V1Hat;
  real<lower = 0> V2Hat;
  real<lower = 0> kaHat;
  real<lower = 0> mttHat;
  real<lower = 0> circ0Hat;
  real<lower = 0> alphaHat;
  real<lower = 0> gamma;
  real<lower = 0> sigma;
  real<lower = 0> sigmaNeut;
  
  ## IIV parameters
  cholesky_factor_corr[nIIV] L;
  vector<lower = 0>[nIIV] omega;
  matrix[nIIV, nSubjects] etaStd;
  
}

transformed parameters{
  vector[nt] cHat;
  vector[nObsPK] cHatObs;
  vector[nt] neutHat;
  vector[nObsPD] neutHatObs;
  matrix[nt, 8] x;
  vector<lower = 0>[nTheta] parms[1]; # The [1] indicates the parameters are constant
  
  ## variables for Matt's trick
  vector<lower = 0>[nIIV] thetaHat;
  matrix<lower = 0>[nSubjects, nIIV] thetaM; 

  ## Matt's trick to use unit scale
  thetaHat[1] = CLHat; 
  thetaHat[2] = QHat;
  thetaHat[3] = V1Hat;
  thetaHat[4] = V2Hat;
  thetaHat[5] = mttHat;
  thetaHat[6] = circ0Hat;
  thetaHat[7] = alphaHat;
  thetaM = (rep_matrix(thetaHat, nSubjects) .* 
             exp(diag_pre_multiply(omega, L * etaStd)))';
  
  for(i in 1:8) {
    parms[1][9 + i] = 1; # F
    parms[1][17 + i] = 0; # tlag
  }

  for(i in 1:nSubjects) {

    parms[1][1] = thetaM[i, 1] * (weight[i] / 70)^0.75; # CL
    parms[1][2] = thetaM[i, 2] * (weight[i] / 70)^0.75; # Q
    parms[1][3] = thetaM[i, 3] * (weight[i] / 70); # V1
    parms[1][4] = thetaM[i, 4] * (weight[i] / 70); # V2
    parms[1][5] = kaHat; # ka
    parms[1][6] = thetaM[i, 5]; # mtt
    parms[1][7] = thetaM[i, 6]; # circ0
    parms[1][8] = gamma;
    parms[1][9] = thetaM[i, 7]; # alpha

    x[start[i]:end[i]] = generalCptModel_rk45(twoCptNeutModelODE, 8,
                                              parms, 
                                              time[start[i]:end[i]], 
                                              amt[start[i]:end[i]], 
                                              rate[start[i]:end[i]], 
                                              ii[start[i]:end[i]], 
                                              evid[start[i]:end[i]], 
                                              cmt[start[i]:end[i]], 
                                              addl[start[i]:end[i]], 
                                              ss[start[i]:end[i]],
                                              1e-6, 1e-6, 1e8);
                             
    cHat[start[i]:end[i]] = x[start[i]:end[i], 2] / parms[1][3]; # divide by V1
    neutHat[start[i]:end[i]] = x[start[i]:end[i], 8] + parms[1][7]; # Add baseline
    
  }
  
  for(i in 1:nObsPK) cHatObs[i] = cHat[iObsPK[i]];
  for(i in 1:nObsPD) neutHatObs[i] = neutHat[iObsPD[i]];

}

model{
  
  ## Priors
  CLHat ~ normal(0, 20);
  QHat ~ normal(0, 20);
  V1Hat ~ normal(0, 100);
  V2Hat ~ normal(0, 1000);
  kaHat ~ normal(0, 5);
  sigma ~ cauchy(0, 1);
  
  mttHat ~ lognormal(log(mttHatPrior), mttHatPriorCV);
  circ0Hat ~ lognormal(log(circ0HatPrior), circ0HatPriorCV);
  alphaHat ~ lognormal(log(alphaHatPrior), alphaHatPriorCV);
  gamma ~ lognormal(log(gammaPrior), gammaPriorCV);
  sigmaNeut ~ cauchy(0, 1);

  ## Parameters for Matt's trick
  L ~ lkj_corr_cholesky(1);
  to_vector(etaStd) ~ normal(0, 1);
  omega ~ cauchy(0, 1);

  ## observed data likelihood
  logCObs ~ normal(log(cHatObs), sigma);
  logNeutObs ~ normal(log(neutHatObs), sigmaNeut);
}

generated quantities{

  matrix[nt, 8] xPred;
  vector<lower = 0>[nTheta] parmsPred[1]; # [1] indicates the parameters are constant
  vector[nt] cHatPred;
  vector[nt] neutHatPred;
  vector<lower = 0>[nObsPK] cHatObsCond;
  vector<lower = 0>[nObsPK] cHatObsPred;
  vector<lower = 0>[nObsPD] neutHatObsCond;
  vector<lower = 0>[nObsPD] neutHatObsPred;

  ## Variables for IIV  
  matrix[nIIV, nSubjects] etaStdPred;
  matrix<lower = 0>[nSubjects, nIIV] thetaPredM;
  corr_matrix[nIIV] rho;
  
  rho = L * L';
  for(i in 1:nSubjects) {
    for(j in 1:nIIV) {
      etaStdPred[j, i] = normal_rng(0, 1);
    }
  }
  thetaPredM = (rep_matrix(thetaHat, nSubjects) .* 
                exp(diag_pre_multiply(omega, L * etaStdPred)))';
                
  for(i in 1:8) {
    parmsPred[1][9 + i] = 1; # F
    parmsPred[1][17 + i] = 0; # tlag
  }
  
  for(i in 1:nSubjects) {
    parmsPred[1][1] = thetaPredM[i, 1] * (weight[i] / 70)^0.75; # CL
    parmsPred[1][2] = thetaPredM[i, 2] * (weight[i] / 70)^0.75; # Q
    parmsPred[1][3] = thetaPredM[i, 3] * (weight[i] / 70); # V1
    parmsPred[1][4] = thetaPredM[i, 4] * (weight[i] / 70); # V2
    parmsPred[1][5] = kaHat; # ka
    parmsPred[1][6] = thetaPredM[i, 5]; # mtt
    parmsPred[1][7] = thetaPredM[i, 6]; # circ0
    parmsPred[1][8] = gamma; # gamma
    parmsPred[1][9] = thetaPredM[i, 7]; # alpha
    
    xPred[start[i]:end[i]] = generalCptModel_rk45(twoCptNeutModelODE, 8,
                                                    parmsPred, 
                                                    time[start[i]:end[i]], 
                                                    amt[start[i]:end[i]],
                                                    rate[start[i]:end[i]],
                                                    ii[start[i]:end[i]],
                                                    evid[start[i]:end[i]],
                                                    cmt[start[i]:end[i]],
                                                    addl[start[i]:end[i]],
                                                    ss[start[i]:end[i]],
                                                    1e-6, 1e-6, 1e8);
    
    cHatPred[start[i]:end[i]] = xPred[start[i]:end[i], 2] / parmsPred[1][3]; # divide by V1
    neutHatPred[start[i]:end[i]] = xPred[start[i]:end[i], 8] + parmsPred[1][7]; # Add baseline
  }

  ## predictions for observed data records
  for(i in 1:nObsPK) cHatObsPred[i] = cHatPred[iObsPK[i]];
  for(i in 1:nObsPD) neutHatObsPred[i] = neutHatPred[iObsPD[i]];
  
  for(i in 1:nObsPK) {
    cHatObsCond[i] = exp(normal_rng(log(fmax(machine_precision(), cHatObs[i])), sigma));
    cHatObsPred[i] = exp(normal_rng(log(fmax(machine_precision(), cHatObsPred[i])), sigma));
  }
  
  for(i in 1:nObsPD) {
    neutHatObsCond[i] = exp(normal_rng(log(fmax(machine_precision(), neutHatObs[i])), sigmaNeut));
    neutHatObsPred[i] = exp(normal_rng(log(fmax(machine_precision(), neutHatObsPred[i])), sigmaNeut));
  }
}
