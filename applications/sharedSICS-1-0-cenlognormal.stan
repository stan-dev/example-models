data {
  int<lower=2> N; //Num datapoints
  int<lower=0, upper=N> Nobs;
  int<lower=0, upper=N> Nmiss;

  int<lower=1, upper=N> obsInd[Nobs];
  int<lower=1, upper=N> missInd[Nmiss];

  int<lower=1> nFeatures;
  int<lower=1> nProteins;
  int<lower=1> nPeptides;
  int<lower=1> nPeptideRuns;
  int<lower=2> nDigests;
  int<lower=2> nDigestEffects;
  int<lower=1> nRuns;
  int<lower=2> nSamples;
  int<lower=1> nPeptideVarianceLevels;
  int<lower=1> nResidualVarianceLevels;

  int<lower=1> nNormalisationLevels;
  vector[nNormalisationLevels] normalisation;
  int<lower=1,upper=nNormalisationLevels> normalisationMap[N];

  //Parameters for scaled-inv-chi-squared prior on peptide variance
  vector<lower=0>[nPeptideVarianceLevels] peptideNu;
  vector<lower=0>[nPeptideVarianceLevels] peptideTau;

  //Parameters for scaled-inv-chi-squared prior on residual variance  
  vector<lower=0>[nResidualVarianceLevels] residualNu;
  vector<lower=0>[nResidualVarianceLevels] residualTau;

  real proteinRefPriorMean;
  real<lower=0> proteinRefPriorSD;
  real<lower=0> samplePriorSD;

  vector[nFeatures+nPeptideRuns] concentrationPrior; //Prior for spectrum effects
  int<lower=0,upper=nFeatures+nPeptideRuns> csAllFeatures[nPeptideRuns+1];
  int<lower=0,upper=nFeatures> csObsFeatures[nPeptideRuns+1];
  int<lower=0,upper=nFeatures> nFeaturesInPeptide[nPeptideRuns];

  int<lower=1,upper=nProteins> proteinMap[nProteins*nSamples];

  //Selection of effect levels for each data point
  int<lower=1,upper=nFeatures+nPeptideRuns> featureMap[N];
  int<lower=1,upper=nDigestEffects> peptideDigestMap[N];

  int<lower=1,upper=nPeptides*nSamples*nProteins*nSamples> nNonZeroElts;
  int<lower=1,upper = nPeptides*nSamples> lhs[nNonZeroElts];
  int<lower=1,upper = nProteins*nSamples> rhs[nNonZeroElts];

  matrix[nSamples*nProteins, (nSamples-1)*nProteins] sampleQR;
  vector<lower=0>[(nSamples-1)*nProteins] sampleRawPriorSD;

  matrix[nDigestEffects,nDigestEffects-nPeptides] digestQR;
  vector[nDigestEffects-nPeptides] digestRawPriorSD;
  
  vector[N] residualRawPriorSD;

  int<lower=1,upper=nPeptides*nSamples> peptideSampleMap[N];

  int<lower=1,upper=nResidualVarianceLevels> residualVarianceMap[N]; //residual variance selection

  int<lower=1,upper=nPeptideVarianceLevels> digestEffectMap[nDigestEffects]; //used only for generated quantities

  int debug;

  vector[Nobs] obs_count; //Log-Counts
  vector[Nmiss] miss_upper; //Missing Log-Counts Upper bounds
}

parameters {
  vector[nProteins] logProteinReferenceAbundance_raw;
  vector[nFeatures] beta_Feature_raw;

  vector[nProteins*(nSamples-1)] beta_Sample_raw;

  vector[nDigestEffects-nPeptides] epsilon_Digest_raw;
  vector<lower=0>[nPeptideVarianceLevels] sigmaSq_Peptide;

  vector<lower=0>[nResidualVarianceLevels] sigmaSq_Residual;

}

transformed parameters {
  vector[nProteins] logProteinReferenceAbundance;
  vector[nProteins*nSamples] logProteinSampleAbundance;
  vector[nPeptides*nSamples] logPeptideSampleAbundance;
  vector[nFeatures+nPeptideRuns] beta_Feature;
  vector[nPeptideRuns] jacAdj;

  vector[nDigestEffects] epsilon_Digest;

  vector[nProteins*nSamples] beta_Sample;

  vector<lower=0>[nPeptideVarianceLevels] sigma_Peptide;
  vector<lower=0>[nResidualVarianceLevels] sigma_Residual;

  logProteinReferenceAbundance = proteinRefPriorMean + proteinRefPriorSD * logProteinReferenceAbundance_raw;

  //Simplex transformed onto K-1 Hypercube, then logit transformed into R^(K-1)
  //Re-implementing Stan's stick breaking for each subvector
  for (p in 1:nPeptideRuns) {
    real log_stick_len = 0.0;
    jacAdj[p] = 0.0;

    for (f in 1:nFeaturesInPeptide[p]){
      //beta_Feature[csAllFeatures[p] + 1 + f] = beta_Feature_raw[csObsFeatures[p]+f] / sumBetaFeature[p];
      real eq_share = -log(nFeaturesInPeptide[p] - f + 1);
      real bfRawAdj = beta_Feature_raw[csObsFeatures[p]+f] + eq_share;
      beta_Feature[csAllFeatures[p] + f] = log_stick_len + log_inv_logit(bfRawAdj);
      jacAdj[p] += log_stick_len;
      jacAdj[p] -= log1p_exp(-bfRawAdj);
      jacAdj[p] -= log1p_exp(bfRawAdj);
      //stick_len -= exp(beta_Feature[csAllFeatures[p] + f]);
      log_stick_len = log_diff_exp(log_stick_len, beta_Feature[csAllFeatures[p] + f]);
    }
    //What's left of the stick has no contribution to Jacobian
    beta_Feature[csAllFeatures[p+1]] = log_stick_len;

  }

  sigma_Peptide = sqrt(sigmaSq_Peptide);
  sigma_Residual = sqrt(sigmaSq_Residual);
  
  epsilon_Digest = sigma_Peptide[digestEffectMap] .* (digestQR * epsilon_Digest_raw);

  beta_Sample = samplePriorSD * (sampleQR * beta_Sample_raw);

  logProteinSampleAbundance = logProteinReferenceAbundance[proteinMap] + beta_Sample;

  //Initialize logPeptideSampleAbundance at -Inf
  logPeptideSampleAbundance = rep_vector(negative_infinity(), nPeptides*nSamples);
  for (i in 1:nNonZeroElts) {
    logPeptideSampleAbundance[lhs[i]] = log_sum_exp(logPeptideSampleAbundance[lhs[i]], logProteinSampleAbundance[rhs[i]]);
  }

}

model {
  logProteinReferenceAbundance_raw ~ normal(0,1);
  for (p in 1:nPeptideRuns) {
    //target += dirichlet_lpdf(beta_Feature[csAllFeatures[p]+1:csAllFeatures[p+1]] | concentrationPrior[csAllFeatures[p]+1:csAllFeatures[p+1]]);
    target += sum((concentrationPrior[csAllFeatures[p]+1:csAllFeatures[p+1]] - 1) .* beta_Feature[csAllFeatures[p]+1:csAllFeatures[p+1]]);
    target += jacAdj[p]; //Jacobian Adjustment
  }

  epsilon_Digest_raw ~ normal(0,digestRawPriorSD);
  sigmaSq_Peptide ~ scaled_inv_chi_square(peptideNu, peptideTau);

  sigmaSq_Residual ~ scaled_inv_chi_square(residualNu, residualTau);


  beta_Sample_raw ~ normal(0,sampleRawPriorSD);

  {
    vector[N] C;
    C = beta_Feature[featureMap] + logPeptideSampleAbundance[peptideSampleMap] + epsilon_Digest[peptideDigestMap];
    C += normalisation[normalisationMap];
    if (debug) {
      print("sigmaResidual: ", sigma_Residual);
      print("sigma_Peptide: ", sigma_Peptide);
      print("epsilon_Digest_raw: ", epsilon_Digest_raw);
      print("logPeptideSampleAbundance: ", logPeptideSampleAbundance);
      print("logProteinReferenceAbundance: ", logProteinReferenceAbundance);
      print("logProteinSampleAbundance: ", logProteinSampleAbundance);
      print("beta_Sample: ", beta_Sample);
      print("beta_Feature (raw): ", beta_Feature_raw);
      print("beta_Feature: ", exp(beta_Feature));
      print("jacAdj: ", jacAdj);
      print("C: ", C);
    }
    obs_count ~ normal(C[obsInd],sigma_Residual[residualVarianceMap[obsInd]]);
    if (Nmiss > 0) {
      target += normal_lcdf(miss_upper | C[missInd], sigma_Residual[residualVarianceMap[missInd]]);
    }
  }

}


generated quantities {
}
