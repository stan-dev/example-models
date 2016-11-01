
data{
  int<lower = 1> nSubjects;
  int<lower = 1> nt;
  int<lower = 1> nObs;
  int<lower = 1> iObs[nObs];
  real<lower = 0> amt[nt];
  int<lower = 1> cmt[nt];
  int<lower = 0> evid[nt];
  int<lower = 1> start[nSubjects];
  int<lower = 1> end[nSubjects];
  real<lower = 0> time[nt];
  vector<lower = 0>[nObs] cObs;
  vector[nObs] respObs;
  real<lower = 0> weight[nSubjects];
  real<lower = 0> rate[nt];
  real<lower = 0> ii[nt];
  int<lower = 0> addl[nt];
  int<lower = 0> ss[nt];
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
  real<lower = 0> ke0Hat;
  real<lower = 0> EC50Hat;
  vector<lower = 0>[nRandom] omega;
  corr_matrix[nRandom] rho;
  real<lower = 0> omegaKe0;
  real<lower = 0> omegaEC50;
  real<lower = 0> sigma;
  real<lower = 0> sigmaResp;
  vector[nRandom] logtheta[nSubjects];
  real logKe0[nSubjects];
  real logEC50[nSubjects];
}

transformed parameters{
  vector<lower = 0>[nRandom] thetaHat;
  cov_matrix[nRandom] Omega;
  real<lower = 0> CL[nSubjects];
  real<lower = 0> Q[nSubjects];
  real<lower = 0> V1[nSubjects];
  real<lower = 0> V2[nSubjects];
  real<lower = 0> ka[nSubjects];
  real<lower = 0> ke0[nSubjects];
  real<lower = 0> EC50[nSubjects];
  vector<lower = 0>[8] parms[1];
  matrix[4, 4] K;
  real k10;
  real k12;
  real k21;
  vector<lower = 0>[nt] cHat;
  vector<lower = 0>[nObs] cHatObs;
  vector<lower = 0>[nt] respHat;
  vector<lower = 0>[nObs] respHatObs;
  vector<lower = 0>[nt] ceHat;
  matrix[nt, 4] x;
  
  thetaHat[1] = CLHat;
  thetaHat[2] = QHat;
  thetaHat[3] = V1Hat;
  thetaHat[4] = V2Hat;
  thetaHat[5] = kaHat;

  Omega = quad_form_diag(rho, omega); ## diag_matrix(omega) * rho * diag_matrix(omega)

  for(i in 1:4){
    parms[1][i] = 1; # F
    parms[1][4 + i] = 0; # tlag
  }

  for(j in 1:nSubjects){
    CL[j] = exp(logtheta[j, 1]) * (weight[j] / 70)^0.75;
    Q[j] = exp(logtheta[j, 2]) * (weight[j] / 70)^0.75;
    V1[j] = exp(logtheta[j, 3]) * weight[j] / 70;
    V2[j] = exp(logtheta[j, 4]) * weight[j] / 70;
    ka[j] = exp(logtheta[j, 5]);
    ke0[j] = exp(logKe0[j]);
    EC50[j] = exp(logEC50[j]);
    
    k10 = CL[j] / V1[j];
    k12 = Q[j] / V1[j];
    k21 = Q[j] / V2[j];

    K = rep_matrix(0, 4, 4);
    
    K[1, 1] = -ka[j];
    K[2, 1] = ka[j];
    K[2, 2] = -(k10 + k12);
    K[2, 3] = k21;
    K[3, 2] = k12;
    K[3, 3] = -k21;
    K[4, 2] = ke0[j];
    K[4, 4] = -ke0[j];
    
    x[start[j]:end[j],] = linCptModel(K, parms,
                                      time[start[j]:end[j]],
                                      amt[start[j]:end[j]],
                                      rate[start[j]:end[j]],
                                      ii[start[j]:end[j]],
                                      evid[start[j]:end[j]],
                                      cmt[start[j]:end[j]],
                                      addl[start[j]:end[j]],
                                      ss[start[j]:end[j]]);

    cHat[start[j]:end[j]] = 1000 * x[start[j]:end[j], 2] ./ V1[j];
    ceHat[start[j]:end[j]] = 1000 * x[start[j]:end[j], 4] ./ V1[j];
    respHat[start[j]:end[j]] = 100 * ceHat[start[j]:end[j]] ./ 
       (EC50[j] + ceHat[start[j]:end[j]]);
  }

  for(i in 1:nObs){
    cHatObs[i] = cHat[iObs[i]];
    respHatObs[i] = respHat[iObs[i]];
  }
}

model{
    CLHat ~ normal(0, 20);
    QHat ~ normal(0, 40);
    V1Hat ~ normal(0, 150);
    V2Hat ~ normal(0, 150);
    kaHat ~ normal(0, 5);
    ke0Hat ~ normal(0, 2);
    EC50Hat ~ normal(0, 200);
    omega ~ cauchy(0, 2);
    rho ~ lkj_corr(1); 
    omegaKe0 ~ cauchy(0, 2);
    omegaEC50 ~ cauchy(0, 2);
    sigma ~ cauchy(0, 2);
    sigmaResp ~ cauchy(0, 5);

    ## Inter-individual variability
    logtheta ~ multi_normal(log(thetaHat), Omega);
    logKe0 ~ normal(log(ke0Hat), omegaKe0);
    logEC50 ~ normal(log(EC50Hat), omegaEC50);

    logCObs ~ normal(log(cHatObs), sigma); 
    respObs ~ normal(respHatObs, sigmaResp); 
}

generated quantities{
  vector[nRandom] logthetaPred[nSubjects];
  real logKe0Pred[nSubjects];
  real logEC50Pred[nSubjects];
  vector<lower = 0>[nt] cHatPred;
  vector<lower = 0>[nt] cObsCond;
  vector<lower = 0>[nt] cObsPred;
  vector<lower = 0>[nt] respHatPred;
  vector<lower = 0>[nt] ceHatPred;
  vector[nt] respObsCond;
  vector[nt] respObsPred;
  real<lower = 0> CLPred[nSubjects];
  real<lower = 0> QPred[nSubjects];
  real<lower = 0> V1Pred[nSubjects];
  real<lower = 0> V2Pred[nSubjects];
  real<lower = 0> kaPred[nSubjects];
  real<lower = 0> ke0Pred[nSubjects];
  real<lower = 0> EC50Pred[nSubjects];
  matrix[nt, 4] xPred;
  vector<lower = 0>[8] parmsPred[1];
  matrix[4, 4] KPred;
  real k10Pred;
  real k12Pred;
  real k21Pred;

  for(i in 1:4){
    parmsPred[1][i] = 1; # F
    parmsPred[1][4 + i] = 0; # tlag
  }

  for(j in 1:nSubjects){
    logthetaPred[j] = multi_normal_rng(log(thetaHat), Omega);
    logKe0Pred[j] = normal_rng(log(ke0Hat), omegaKe0);
    logEC50Pred[j] = normal_rng(log(EC50Hat), omegaEC50);
    CLPred[j] = exp(logthetaPred[j, 1]) * (weight[j] / 70)^0.75;
    QPred[j] = exp(logthetaPred[j, 2]) * (weight[j] / 70)^0.75;
    V1Pred[j] = exp(logthetaPred[j, 3]) * weight[j] / 70;
    V2Pred[j] = exp(logthetaPred[j, 4]) * weight[j] / 70;
    kaPred[j] = exp(logthetaPred[j, 5]);
    ke0Pred[j] = exp(logKe0Pred[j]);
    EC50Pred[j] = exp(logEC50Pred[j]);
    
    k10Pred = CLPred[j] / V1Pred[j];
    k12Pred = QPred[j] / V1Pred[j];
    k21Pred = QPred[j] / V2Pred[j];

    KPred = rep_matrix(0, 4, 4);
    
    KPred[1, 1] = -kaPred[j];
    KPred[2, 1] = kaPred[j];
    KPred[2, 2] = -(k10Pred + k12Pred);
    KPred[2, 3] = k21Pred;
    KPred[3, 2] = k12Pred;
    KPred[3, 3] = -k21Pred;
    KPred[4, 2] = ke0Pred[j];
    KPred[4, 4] = -ke0Pred[j];
    
    xPred[start[j]:end[j],] = linCptModel(KPred, parmsPred,
                                      time[start[j]:end[j]],
                                      amt[start[j]:end[j]],
                                      rate[start[j]:end[j]],
                                      ii[start[j]:end[j]],
                                      evid[start[j]:end[j]],
                                      cmt[start[j]:end[j]],
                                      addl[start[j]:end[j]],
                                      ss[start[j]:end[j]]);    

    cHatPred[start[j]:end[j]] = 1000 * xPred[start[j]:end[j], 2] ./ V1Pred[j];
    ceHatPred[start[j]:end[j]] = 1000 * xPred[start[j]:end[j], 4] ./ V1Pred[j];
    respHatPred[start[j]:end[j]] = 100 * ceHatPred[start[j]:end[j]] ./
      (EC50Pred[j] + ceHatPred[start[j]:end[j]]);
  }
  
  for(i in 1:nt){
    if(time[i] == 0){
      cObsCond[i] = 0;
      cObsPred[i] = 0;
      respObsCond[i] = 0;
      respObsPred[i] = 0;
    }else{
      cObsCond[i] = exp(normal_rng(log(cHat[i]), sigma));
      cObsPred[i] = exp(normal_rng(log(cHatPred[i]), sigma));
      respObsCond[i] = normal_rng(respHat[i], sigmaResp);
      respObsPred[i] = normal_rng(respHatPred[i], sigmaResp);
    }
  }

}

