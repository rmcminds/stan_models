functions {
    matrix log_prod_exp_mats(matrix logA, matrix logB) {
      matrix[rows(logA), cols(logB)] logX;
      for (i in 1:rows(logX)) {
        for (j in 1:cols(logX)) {
          logX[i, j] = log_sum_exp(logA[i, ] + logB[, j]');
        }
      }
      return logX;
    }
    vector log_prod_exp_mat_vect(matrix logA, vector logB) {
      vector[rows(logA)] logX;
      row_vector[rows(logB)] logBPrime = logB';
      for (i in 1:rows(logA)) {
        logX[i] = log_sum_exp(logA[i, ] + logBPrime);
      }
      return logX;
    }
    matrix log_outer_exp(vector logA, row_vector logB) {
      return rep_matrix(logA, cols(logB)) + rep_matrix(logB, rows(logA));
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
    matrix[NMicrobeNodes, NMicrobeNodes] microbeAncestors;
    matrix[NEffects, NFactors] factLevelMat;
    matrix[NSamples, NEffects] modelMat;
    matrix[NHostNodes, NTimeBins] edgeToBin;
    matrix[NHostNodes, NHostNodes] hostAncestors;
    matrix[NSamples, NHostNodes] hostAncestorsExpanded;
    row_vector[NMicrobeNodes] microbeEdges;
}
transformed data {
    matrix[NHostNodes, NTimeBins] logEdgeToBin = log(edgeToBin);
    row_vector[NMicrobeNodes] logMicrobeEdges = log(microbeEdges);
    matrix[NHostTips, NHostNodes] logHostAncestors = log(hostAncestors[1:NHostTips, ]);
    matrix[NMicrobeNodes, NMicrobeTips] logMicrobeAncestors = log(microbeAncestors[, 1:NMicrobeTips]);
    matrix[NSamples, 1 + NEffects + NHostNodes] fullModelMat = append_col(rep_vector(1, NSamples),
                                                                append_col(modelMat,
                                                                 hostAncestorsExpanded));
    vector[rows(csr_extract_w(fullModelMat))] wFullModelMat = csr_extract_w(fullModelMat);
    int vFullModelMat[size(csr_extract_v(fullModelMat))] = csr_extract_v(fullModelMat);
    int uFullModelMat[size(csr_extract_u(fullModelMat))] = csr_extract_u(fullModelMat);
    matrix[NMicrobeTips, NMicrobeNodes] microbeAncestorsForTipsT = microbeAncestors[, 1:NMicrobeTips]';
    vector[rows(csr_extract_w(microbeAncestorsForTipsT))] wMicrobeAncestorsForTipsT = csr_extract_w(microbeAncestorsForTipsT);
    int vMicrobeAncestorsForTipsT[size(csr_extract_v(microbeAncestorsForTipsT))] = csr_extract_v(microbeAncestorsForTipsT);
    int uMicrobeAncestorsForTipsT[size(csr_extract_u(microbeAncestorsForTipsT))] = csr_extract_u(microbeAncestorsForTipsT);
}
parameters {
    real<lower=0> aveStD;
    simplex[2 * NFactors + 3] stDProps;
    vector[NTimeBins - 1] timeBinMetaVar;
    real<lower=0> aveStDMeta;
    simplex[3] metaVarProps;
    simplex[2] hostMetaVarProps;
    row_vector[NMicrobeNodes - NMicrobeTips] phyloLogVarMultPrev;
    vector[NHostNodes - NHostTips] phyloLogVarMultADiv;
    matrix[NHostNodes - NHostTips, NMicrobeNodes - NMicrobeTips] phyloLogVarMultRaw;
    real globalIntercept;
    vector[NEffects + NHostNodes] rawAlphaDivEffects;
    matrix[NEffects + NHostNodes + 1, NMicrobeNodes] rawMicrobeNodeEffects;
}
transformed parameters {
    vector<lower=0>[2 * NFactors + 3] scales;
    vector<lower=0>[3] metaScales;
    vector[NTimeBins] logRelativeEvolRates;
    vector<lower=0>[NTimeBins] relativeEvolRates;
    row_vector[NMicrobeNodes] logMicrobeVarRaw;
    row_vector<lower=0>[NMicrobeNodes] microbeVarRaw;
    row_vector<lower=0>[NMicrobeNodes] microbeScales;
    vector[NHostNodes] logHostVarRaw;
    vector<lower=0>[NHostNodes] hostScales;
    matrix[NHostNodes, NMicrobeNodes] logPhyloScalesRaw;
    matrix<lower=0>[NHostNodes, NMicrobeNodes] phyloScales;
    vector[NEffects + NHostNodes + 1] scaledAlphaDivEffects;
    matrix[NEffects + NHostNodes + 1, NMicrobeNodes] scaledMicrobeNodeEffects;
    scales
        = sqrt((2 * NFactors + 3) * stDProps)
          * aveStD;
    metaScales
        = sqrt(3 * metaVarProps)
          * aveStDMeta;
    logMicrobeVarRaw
        = metaScales[1] * phyloLogVarMultPrev
              * microbeAncestors[(NMicrobeTips + 1):, ]
          + logMicrobeEdges;
    microbeVarRaw
        = exp(logMicrobeVarRaw);
    microbeScales
        = sqrt(microbeVarRaw
               / mean(csr_matrix_times_vector(NMicrobeTips,
                                              NMicrobeNodes,
                                              wMicrobeAncestorsForTipsT,
                                              vMicrobeAncestorsForTipsT,
                                              uMicrobeAncestorsForTipsT,
                                              microbeVarRaw')));
    logRelativeEvolRates
        = append_row(0,
            cumulative_sum(timeBinMetaVar)
                * metaScales[2]
                * sqrt(hostMetaVarProps[1]));
    relativeEvolRates
        = exp(logRelativeEvolRates);
    logHostVarRaw
        = hostAncestors[, (NHostTips + 1):]
              * phyloLogVarMultADiv
                * metaScales[2]
                * sqrt(hostMetaVarProps[2])
          + log_prod_exp_mat_vect(logEdgeToBin,
                                  logRelativeEvolRates);
    hostScales
        = exp(log(scales[2 * NFactors + 1])
              + (logHostVarRaw
                  - log_sum_exp(log_prod_exp_mat_vect(logHostAncestors,
                                                      logHostVarRaw))
                  + log(NHostTips))
                 * 0.5);
    logPhyloScalesRaw
        = hostAncestors[, (NHostTips + 1):]
              * (phyloLogVarMultRaw * metaScales[3])
              * microbeAncestors[(NMicrobeTips + 1):, ]
          + log_outer_exp(logHostVarRaw, logMicrobeVarRaw);
    phyloScales
        = exp(log(scales[2 * NFactors + 2])
              + (logPhyloScalesRaw
                 - log_sum_exp(log_prod_exp_mats(log_prod_exp_mats(logHostAncestors,
                                                                   logPhyloScalesRaw),
                                                 logMicrobeAncestors))
                 + log(NHostTips * NMicrobeTips))
                * 0.5); // standardization as for hostScales, on log scale for stability
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
    matrix[NEffects + NHostNodes + 1, NMicrobeTips] microbePreMult;
    matrix[NSamples, NMicrobeTips] sampleTipEffects;
    vector[NObs] logit_ratios;
    aveStD ~ exponential(1.0 / aveStDPriorExpect);
    stDProps ~ dirichlet(rep_vector(1, 2 * NFactors + 3));
    timeBinMetaVar ~ normal(0,1);
    aveStDMeta ~ exponential(1.0);
    metaVarProps ~ dirichlet(rep_vector(1, 3));
    hostMetaVarProps ~ dirichlet(rep_vector(1, 2));
    phyloLogVarMultPrev ~ normal(0,1);
    phyloLogVarMultADiv ~ normal(0,1);
    to_vector(phyloLogVarMultRaw) ~ normal(0,1);
    globalIntercept ~ normal(0, 50);
    rawAlphaDivEffects ~ normal(0,1);
    to_vector(rawMicrobeNodeEffects) ~ normal(0,1);
    for (e in 1:(NEffects + NHostNodes + 1)) {
        microbePreMult[e, ]
            = csr_matrix_times_vector(NMicrobeTips,
                                      NMicrobeNodes,
                                      wMicrobeAncestorsForTipsT,
                                      vMicrobeAncestorsForTipsT,
                                      uMicrobeAncestorsForTipsT,
                                      scaledMicrobeNodeEffects[e, ]')';
    }
    for (t in 1:NMicrobeTips) {
        sampleTipEffects[, t]
            = csr_matrix_times_vector(NSamples,
                                      1 + NEffects + NHostNodes,
                                      wFullModelMat,
                                      vFullModelMat,
                                      uFullModelMat,
                                      scaledAlphaDivEffects
                                       + microbePreMult[, t]);
    }
    for (n in 1:NObs)
        logit_ratios[n] = sampleTipEffects[sampleNames[n], microbeTipNames[n]];
    present ~ bernoulli_logit(logit_ratios);
}
