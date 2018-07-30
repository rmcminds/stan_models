functions {
    real log_mean_prod_exp_csr_mat_vect(int m, vector w, int[] v, int[] u, vector d) {
      vector[rows(w)] logX;
      real logMeanLength;
      for (i in 1:m) {
        logX[u[i]:(u[i + 1] - 1)] = w[u[i]:(u[i + 1] - 1)] + d[v[u[i]:(u[i + 1] - 1)]];
      }
      logMeanLength = log_sum_exp(logX) - log(m);
      return logMeanLength;
    }
    vector log_prod_exp_mat_vect(matrix logA, vector logB) {
      vector[rows(logA)] logX;
      row_vector[rows(logB)] logBPrime = logB';
      for (i in 1:rows(logA)) {
        logX[i] = log_sum_exp(logA[i, ] + logBPrime);
      }
      return logX;
    }
    matrix log_outer_exp(vector logA, vector logB) {
      return rep_matrix(logA, rows(logB)) + rep_matrix(logB', rows(logA));
    }
}
data {
    int NSamples;
    int NObs;
    int NMicrobeNodes;
    int NMicrobeTips;
    int NFactors;
    int NEffects;
    int NHostNodes;
    int NHostTips;
    int NTimeBins;
    int present[NObs];
    int sampleNames[NObs];
    int microbeTipNames[NObs];
    real<lower=0> aveStDPriorExpect;
    real<lower=0> aveStDMetaPriorExpect;
    matrix[NEffects, NFactors] factLevelMat;
    matrix<lower=0>[NHostNodes, NTimeBins] edgeToBin;
    vector<lower=0>[NMicrobeNodes] microbeEdges;
    int nnzMicrobeAncestors;
    int nuMicrobeAncestors;
    vector<lower=0>[nnzMicrobeAncestors] wMicrobeAncestors;
    int vMicrobeAncestors[nnzMicrobeAncestors];
    int uMicrobeAncestors[nuMicrobeAncestors];
    int nnzMicrobeTipAncestors;
    int nuMicrobeTipAncestors;
    vector<lower=0>[nnzMicrobeTipAncestors] wMicrobeTipAncestors;
    int vMicrobeTipAncestors[nnzMicrobeTipAncestors];
    int uMicrobeTipAncestors[nuMicrobeTipAncestors];
    int nnzHostAncestors;
    int nuHostAncestors;
    vector<lower=0>[nnzHostAncestors] wHostAncestors;
    int vHostAncestors[nnzHostAncestors];
    int uHostAncestors[nuHostAncestors];
    int nnzHostTipAncestors;
    int nuHostTipAncestors;
    vector<lower=0>[nnzHostTipAncestors] wHostTipAncestors;
    int vHostTipAncestors[nnzHostTipAncestors];
    int uHostTipAncestors[nuHostTipAncestors];
    int nnzHostKronMicrobeAncestors;
    int nuHostKronMicrobeAncestors;
    vector<lower=0>[nnzHostKronMicrobeAncestors] wHostKronMicrobeAncestors;
    int vHostKronMicrobeAncestors[nnzHostKronMicrobeAncestors];
    int uHostKronMicrobeAncestors[nuHostKronMicrobeAncestors];
    int nnzHostKronMicrobeTipAncestors;
    int nuHostKronMicrobeTipAncestors;
    vector<lower=0>[nnzHostKronMicrobeTipAncestors] wHostKronMicrobeTipAncestors;
    int vHostKronMicrobeTipAncestors[nnzHostKronMicrobeTipAncestors];
    int uHostKronMicrobeTipAncestors[nuHostKronMicrobeTipAncestors];
    int nnzFullModelMat;
    int nuFullModelMat;
    vector[nnzFullModelMat] wFullModelMat;
    int vFullModelMat[nnzFullModelMat];
    int uFullModelMat[nuFullModelMat];
}
transformed data {
    vector[NMicrobeNodes] logMicrobeEdges = log(microbeEdges);
    vector[nnzMicrobeTipAncestors] wMicrobeTipAncestorsLog = log(wMicrobeTipAncestors);
    matrix[NHostNodes, NTimeBins] logEdgeToBin = log(edgeToBin);
    vector[nnzHostTipAncestors] wHostTipAncestorsLog = log(wHostTipAncestors);
    vector[nnzHostKronMicrobeTipAncestors] wHostKronMicrobeTipAncestorsLog = log(wHostKronMicrobeTipAncestors);
}
parameters {
    real<lower=0> aveStD;
    simplex[2 * NFactors + 3] stDProps;
    vector[NTimeBins - 1] timeBinMetaVar;
    real<lower=0> aveStDMeta;
    simplex[3] metaVarProps;
    simplex[2] hostMetaVarProps;
    vector[NMicrobeNodes] phyloLogVarMultPrev;
    vector[NHostNodes] phyloLogVarMultADiv;
    vector[NHostNodes * NMicrobeNodes] phyloLogVarMultRaw;
    real globalIntercept;
    vector[NEffects + NHostNodes] rawAlphaDivEffects;
    matrix[NEffects + NHostNodes + 1, NMicrobeNodes] rawMicrobeNodeEffects;
}
transformed parameters {
    vector<lower=0>[2 * NFactors + 3] scales;
    vector<lower=0>[3] metaScales;
    vector[NMicrobeNodes] microbeVarRaw;
    row_vector<lower=0>[NMicrobeNodes] microbeScales;
    vector[NTimeBins] logRelativeEvolRates;
    vector[NHostNodes] hostVarRaw;
    vector<lower=0>[NHostNodes] hostScales;
    matrix[NHostNodes, NMicrobeNodes] phyloVarRaw;
    matrix<lower=0>[NHostNodes, NMicrobeNodes] phyloScales;
    vector[NEffects + NHostNodes + 1] scaledAlphaDivEffects;
    matrix[NEffects + NHostNodes + 1, NMicrobeNodes] scaledMicrobeNodeEffects;
    scales
        = sqrt((2 * NFactors + 3) * stDProps)
          * aveStD;
    metaScales
        = sqrt(3 * metaVarProps)
          * aveStDMeta;
    microbeVarRaw
        = csr_matrix_times_vector(NMicrobeNodes,
                                  NMicrobeNodes,
                                  wMicrobeAncestors,
                                  vMicrobeAncestors,
                                  uMicrobeAncestors,
                                  phyloLogVarMultPrev
                                    * metaScales[1])
          + logMicrobeEdges;
    microbeScales
        = exp(0.5 * (microbeVarRaw
                     - log_mean_prod_exp_csr_mat_vect(NMicrobeTips,
                                                      wMicrobeTipAncestorsLog,
                                                      vMicrobeTipAncestors,
                                                      uMicrobeTipAncestors,
                                                      microbeVarRaw)))';
    logRelativeEvolRates
        = append_row(0,
            cumulative_sum(timeBinMetaVar)
                * metaScales[2]
                * sqrt(hostMetaVarProps[1]));
    hostVarRaw
        = csr_matrix_times_vector(NHostNodes,
                                  NHostNodes,
                                  wHostAncestors,
                                  vHostAncestors,
                                  uHostAncestors,
                                  phyloLogVarMultADiv
                                    * metaScales[2]
                                    * sqrt(hostMetaVarProps[2]))
          + log_prod_exp_mat_vect(logEdgeToBin, logRelativeEvolRates);
    hostScales
        = exp(log(scales[2 * NFactors + 1])
              + 0.5 * (hostVarRaw
                       - log_mean_prod_exp_csr_mat_vect(NHostTips,
                                                        wHostTipAncestorsLog,
                                                        vHostTipAncestors,
                                                        uHostTipAncestors,
                                                        hostVarRaw)));
    phyloVarRaw
        = to_matrix(csr_matrix_times_vector(NHostNodes * NMicrobeNodes,
                                            NHostNodes * NMicrobeNodes,
                                            wHostKronMicrobeAncestors,
                                            vHostKronMicrobeAncestors,
                                            uHostKronMicrobeAncestors,
                                            phyloLogVarMultRaw
                                              * metaScales[3]),
                    NHostNodes,
                    NMicrobeNodes)
          + log_outer_exp(hostVarRaw, microbeVarRaw);
    phyloScales
        = exp(log(scales[2 * NFactors + 2])
              + 0.5 * (phyloVarRaw
                       - log_mean_prod_exp_csr_mat_vect(NHostTips * NMicrobeTips,
                                                        wHostKronMicrobeTipAncestorsLog,
                                                        vHostKronMicrobeTipAncestors,
                                                        uHostKronMicrobeTipAncestors,
                                                        to_vector(phyloVarRaw))));
    scaledAlphaDivEffects
        = append_row(globalIntercept,
            append_row(factLevelMat * segment(scales, 1, NFactors),
                hostScales)
            .* rawAlphaDivEffects);
    scaledMicrobeNodeEffects
        = append_row(
            append_row(
                scales[2 * NFactors + 3],
                factLevelMat * segment(scales, NFactors + 1, NFactors))
            * microbeScales,
            phyloScales)
          .* rawMicrobeNodeEffects;
}
model {
    matrix[NSamples, NMicrobeTips] sampleTipEffects;
    vector[NObs] logit_ratios;
    aveStD ~ exponential(1.0 / aveStDPriorExpect);
    stDProps ~ dirichlet(rep_vector(1, 2 * NFactors + 3));
    timeBinMetaVar ~ normal(0,1);
    aveStDMeta ~ exponential(1.0 / aveStDMetaPriorExpect);
    metaVarProps ~ dirichlet(rep_vector(1, 3));
    hostMetaVarProps ~ dirichlet(rep_vector(1, 2));
    phyloLogVarMultPrev ~ normal(0,1);
    phyloLogVarMultADiv ~ normal(0,1);
    phyloLogVarMultRaw ~ normal(0,1);
    globalIntercept ~ normal(0, 50);
    rawAlphaDivEffects ~ normal(0,1);
    to_vector(rawMicrobeNodeEffects) ~ normal(0,1);
    sampleTipEffects = to_matrix(csr_matrix_times_vector(NSamples * NMicrobeTips,
                                                         (1 + NEffects + NHostNodes) * NMicrobeNodes,
                                                         wFullModelMat,
                                                         vFullModelMat,
                                                         uFullModelMat,
                                                         to_vector(scaledMicrobeNodeEffects
                                                                   + rep_matrix(scaledAlphaDivEffects, NMicrobeNodes))),
                                 NSamples, NMicrobeTips);
    for (n in 1:NObs)
        logit_ratios[n] = sampleTipEffects[sampleNames[n], microbeTipNames[n]];
    present ~ bernoulli_logit(logit_ratios);
}
