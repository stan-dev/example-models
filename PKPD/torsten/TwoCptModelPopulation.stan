## TwoCptModelPopulation.stan
## v1.3 - example



data{
  int<lower = 1> nt;
  int<lower = 1> nObs;
  int<lower = 1> nSubjects;
  int nIIV;
  int<lower = 1> iObs[nObs];
  int<lower = 1> start[nSubjects];
  int<lower = 1> end[nSubjects];
  int<lower = 1> cmt[nt];
  int evid[nt];
  int addl[nt];
  int ss[nt];
  real amt[nt];
  real time[nt];
  real rate[nt];
  real ii[nt];
  vector<lower = 0>[nObs] cObs;
}

transformed data{
  vector[nObs] logCObs;
  int nTheta;
  int nti[nSubjects];
  int ntMax;

  logCObs = log(cObs);
  nTheta = 11;
  
  for(i in 1:nSubjects) nti[i] = end[i] - start[i] + 1;
}

parameters{
  real<lower = 0, upper = 500> CLHat;
  real<lower = 0, upper = 500> QHat;
  real<lower = 0, upper = 3500> V1Hat;
  real<lower = 0, upper = 3500> V2Hat;
  real<lower = 0, upper = 100> kaHat;
  real<lower = 0> sigma;
  
  # Inter-Individual variability
  cholesky_factor_corr[nIIV] L;
  vector<lower = 0.01, upper = 2>[nIIV] omega;
  matrix[nIIV, nSubjects] etaStd;
}

transformed parameters{
  vector<lower = 0>[nIIV] thetaHat;
  matrix<lower = 0>[nSubjects, nIIV] thetaM; // variable required for Matt's trick
  vector<lower = 0>[nTheta] theta[1];
  vector<lower = 0>[nt] cHat;
  vector<lower = 0>[nObs] cHatObs;
  matrix<lower = 0>[nt, 3] x;

  thetaHat[1] = CLHat;
  thetaHat[2] = QHat;
  thetaHat[3] = V1Hat;
  thetaHat[4] = V2Hat;
  thetaHat[5] = kaHat;

  ## Matt's trick to use unit scale 
  thetaM = (rep_matrix(thetaHat, nSubjects) .* exp(diag_pre_multiply(omega, L * etaStd)))'; 
  
  for(j in 1:nSubjects)
  {
    theta[1][1] = thetaM[j, 1]; # CL
    theta[1][2] = thetaM[j, 2]; # Q
    theta[1][3] = thetaM[j, 3]; # V1
    theta[1][4] = thetaM[j, 4]; # V2
    theta[1][5] = thetaM[j, 5]; # ka
    theta[1][6] = 1; # F1
    theta[1][7] = 1; # F2
    theta[1][8] = 1; # F3
    theta[1][9] = 0; # tlag1
    theta[1][10] = 0; # tlag2
    theta[1][11] = 0; # tlag3

    x[start[j]:end[j]] = PKModelTwoCpt(theta, 
                                       time[start[j]:end[j]], 
                                       amt[start[j]:end[j]],
                                       rate[start[j]:end[j]],
                                       ii[start[j]:end[j]],
                                       evid[start[j]:end[j]],
                                       cmt[start[j]:end[j]],
                                       addl[start[j]:end[j]],
                                       ss[start[j]:end[j]]);
                                       
    cHat[start[j]:end[j]] = col(x[start[j]:end[j]], 2) ./ theta[1][3]; ## divide by V1
  }

  for(i in 1:nObs){
    cHatObs[i] = cHat[iObs[i]]; ## predictions for observed data records
  }
}

model{
  ## Prior
  CLHat ~ lognormal(log(10), 0.25);
  QHat ~ lognormal(log(15), 0.5);
  V1Hat ~ lognormal(log(35), 0.25);
  V2Hat ~ lognormal(log(105), 0.5);
  kaHat ~ lognormal(log(2.5), 1);
  
  L ~ lkj_corr_cholesky(1);
  
  ## Inter-individual variability (see transformed parameters block
  ## for translation to PK parameters)
  to_vector(etaStd) ~ normal(0, 1);
  
  sigma ~ cauchy(0, 5);
  logCObs ~ normal(log(cHatObs), sigma);
}

generated quantities{
    vector[nt] cHatPred;
    real logCObsPred[nt];
    matrix[nt, 3] xPred;
    matrix[nIIV, nSubjects] etaStdPred;
    matrix<lower=0>[nSubjects, nIIV] thetaPredM;
    corr_matrix[nIIV] rho;
    vector<lower = 0>[nTheta] thetaPred[1];
    real logCObsCond[nt];

    rho = L * L';
    
    for(i in 1:nSubjects){
      for(j in 1:nIIV){ 
        etaStdPred[j, i] = normal_rng(0, 1);
      }
    }
    
    thetaPredM = (rep_matrix(thetaHat, nSubjects) .* exp(diag_pre_multiply(omega, L * etaStdPred)))';

    for(j in 1:nSubjects){
      
      thetaPred[1][1] = thetaPredM[j,1]; # CL
      thetaPred[1][2] = thetaPredM[j,2]; # Q 
      thetaPred[1][3] = thetaPredM[j,3]; # V1
      thetaPred[1][4] = thetaPredM[j,4]; # V2
      thetaPred[1][5] = thetaPredM[j,5]; # ka 
      thetaPred[1][6] = 1; # F1
      thetaPred[1][7] = 1; # F2
      thetaPred[1][8] = 1; # F3
      thetaPred[1][9] = 0; # tlag1
      thetaPred[1][10] = 0; # tlag2
      thetaPred[1][11] = 0; # tlag3
    
      xPred[start[j]:end[j],] = PKModelTwoCpt(thetaPred,
                                        time[start[j]:end[j]],
                                        amt[start[j]:end[j]],
                                        rate[start[j]:end[j]],
                                        ii[start[j]:end[j]],
                                        evid[start[j]:end[j]],
                                        cmt[start[j]:end[j]],
                                        addl[start[j]:end[j]],
                                        ss[start[j]:end[j]]);
      for(i in start[j]:end[j]){
        cHatPred[i] = xPred[i, 2] / thetaPred[1][3];
      }
    }
    
    // Need to check the following lines: individual predictions do not work
    for(i in 1:nt){
      if(evid[i] == 1){
        logCObsCond[i] = -99;
        logCObsPred[i] = -99;
      }
      else{
        logCObsCond[i] = normal_rng(log(cHat[i]), sigma);
        logCObsPred[i] = normal_rng(log(cHatPred[i]), sigma);
      }
    }   

}