data{
  int<lower = 1> nSubjects;
  int<lower = 1> nt;
  int<lower = 1> nObs;
  int<lower = 1> iObs[nObs];
  real<lower = 0> amt[nt];
  real<lower = 0> rate[nt];
  real<lower = 0> ii[nt];
  int<lower = 1> cmt[nt];
  int<lower = 0> evid[nt];
  int<lower = 0> addl[nt];
  int<lower = 0> ss[nt];
  int<lower = 1> start[nSubjects];
  int<lower = 1> end[nSubjects];
  real<lower = 0> weight[nSubjects];
  real<lower = 0> time[nt];
  vector<lower = 0>[nObs] cObs;
}

transformed data{
  vector[nObs] logCObs;
  int<lower = 1> nRandom;

  logCObs = log(cObs);

  nRandom = 5;
}

parameters{
  real<lower = 0> CLHat;
  real<lower = 0> QHat;
  real<lower = 0> V1Hat;
  real<lower = 0> V2Hat;
  real<lower = 0> kaHat;
  corr_matrix[nRandom] rho;
  vector<lower = 0>[nRandom] omega;
  real<lower = 0> sigma;
  vector[nRandom] logtheta[nSubjects];
}

transformed parameters{
  vector<lower = 0>[nRandom] thetaHat;
  cov_matrix[nRandom] Omega;
  real<lower = 0> CL[nSubjects];
  real<lower = 0> Q[nSubjects];
  real<lower = 0> V1[nSubjects];
  real<lower = 0> V2[nSubjects];
  real<lower = 0> ka[nSubjects];
  vector<lower = 0>[nt] cHat;
  vector<lower = 0>[nObs] cHatObs;
  matrix<lower = 0>[nt, 3] x;
  vector[11] parms[nSubjects,1];

  thetaHat[1] = CLHat;
  thetaHat[2] = QHat;
  thetaHat[3] = V1Hat;
  thetaHat[4] = V2Hat;
  thetaHat[5] = kaHat;

  Omega = quad_form_diag(rho, omega); ## diag_matrix(omega) * rho * diag_matrix(omega)

  for(j in 1:nSubjects){
    CL[j] = exp(logtheta[j, 1]) * (weight[j] / 70)^0.75;
    Q[j] = exp(logtheta[j, 2]) * (weight[j] / 70)^0.75;
    V1[j] = exp(logtheta[j, 3]) * weight[j] / 70;
    V2[j] = exp(logtheta[j, 4]) * weight[j] / 70;
    ka[j] = exp(logtheta[j, 5]);
    
    parms[j, 1][1] = CL[j];
    parms[j, 1][2] = Q[j];
    parms[j, 1][3] = V1[j];
    parms[j, 1][4] = V2[j];
    parms[j, 1][5] = ka[j];
    parms[j, 1][6] = 1; # F1
    parms[j, 1][7] = 1; # F2
    parms[j, 1][8] = 1; # F3
    parms[j, 1][9] = 0; # tlag1
    parms[j, 1][10] = 0; # tlag2
    parms[j, 1][11] = 0; # tlag3

    x[start[j]:end[j],] = PKModelTwoCpt(parms[j],
					 time[start[j]:end[j]],
					 amt[start[j]:end[j]],
					 rate[start[j]:end[j]],
					 ii[start[j]:end[j]],
					 evid[start[j]:end[j]],
					 cmt[start[j]:end[j]],
					 addl[start[j]:end[j]],
					 ss[start[j]:end[j]]);
    for(i in start[j]:end[j])
      cHat[i] = x[i, 2] / V1[j];
  }

  cHatObs = cHat[iObs]; ## predictions for observed data records

}

model{
    CLHat ~ normal(0, 20);
    QHat ~ normal(0, 20);
    V1Hat ~ normal(0, 100);
    V2Hat ~ normal(0, 1000);
    kaHat ~ normal(0, 5);
    omega ~ cauchy(0, 2);
    rho ~ lkj_corr(1); 
    sigma ~ cauchy(0, 5);

    ## Inter-individual variability
    logtheta ~ multi_normal(log(thetaHat), Omega);

    logCObs ~ normal(log(cHatObs), sigma); ## observed data likelihood
}

generated quantities{
  vector[nRandom] logthetaPred[nSubjects];
  vector<lower = 0>[nt] cHatPred;
  real cObsCond[nt];
  real cObsPred[nt];
  real<lower = 0> CLPred[nSubjects];
  real<lower = 0> QPred[nSubjects];
  real<lower = 0> V1Pred[nSubjects];
  real<lower = 0> V2Pred[nSubjects];
  real<lower = 0> kaPred[nSubjects];
  matrix[nt, 3] xPred;
  vector[11] parmsPred[nSubjects,1];

  for(j in 1:nSubjects){
    logthetaPred[j] = multi_normal_rng(log(thetaHat), Omega);
    CLPred[j] = exp(logthetaPred[j, 1]) * (weight[j] / 70)^0.75;
    QPred[j] = exp(logthetaPred[j, 2]) * (weight[j] / 70)^0.75;
    V1Pred[j] = exp(logthetaPred[j, 3]) * weight[j] / 70;
    V2Pred[j] = exp(logthetaPred[j, 4]) * weight[j] / 70;
    kaPred[j] = exp(logthetaPred[j, 5]);

    parmsPred[j, 1][1] = CLPred[j];
    parmsPred[j, 1][2] = QPred[j];
    parmsPred[j, 1][3] = V1Pred[j];
    parmsPred[j, 1][4] = V2Pred[j];
    parmsPred[j, 1][5] = kaPred[j];
    parmsPred[j, 1][6] = 1; # F1
    parmsPred[j, 1][7] = 1; # F2
    parmsPred[j, 1][8] = 1; # F3
    parmsPred[j, 1][9] = 0; # tlag1
    parmsPred[j, 1][10] = 0; # tlag2
    parmsPred[j, 1][11] = 0; # tlag3

    xPred[start[j]:end[j],] = PKModelTwoCpt(parmsPred[j],
					 time[start[j]:end[j]],
					 amt[start[j]:end[j]],
					 rate[start[j]:end[j]],
					 ii[start[j]:end[j]],
					 evid[start[j]:end[j]],
					 cmt[start[j]:end[j]],
					 addl[start[j]:end[j]],
					 ss[start[j]:end[j]]);

    for(i in start[j]:end[j])
      cHatPred[i] = xPred[i, 2] / V1Pred[j];
  }

  for(i in 1:nt){
    if(time[i] == 0){
      cObsCond[i] = -99;
      cObsPred[i] = -99;
    }else{
      cObsCond[i] = exp(normal_rng(log(cHat[i]), sigma));
      cObsPred[i] = exp(normal_rng(log(cHatPred[i]), sigma));
    }
  }

}
