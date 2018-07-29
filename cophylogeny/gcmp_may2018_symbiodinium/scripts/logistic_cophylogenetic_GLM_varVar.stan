functions {
    vector log_prod_exp_csr_mat_vect(int m, vector w, int[] v, int[] u, vector d) {
      vector[m] logX;
      for (i in 1:m) {
        logX[i] = log_sum_exp(w[u[i]:(u[i+1]-1)] + d[v[u[i]:(u[i+1]-1)]]);
      }
      return logX;
    }
    real mean_of_log_prod_exp_csr_mat_vect(int m, vector w, int[] v, int[] u, vector d) {
      vector[rows(w)] logX;
      real meanLength;
      for (i in 1:m) {
        logX[u[i]:(u[i+1]-1)] = w[u[i]:(u[i+1]-1)] + d[v[u[i]:(u[i+1]-1)]];
      }
      logMeanLength = log_sum_exp(logX) - log(m);
      return logMeanLength;
    }
    matrix log_prod_exp_csr_mat_dense_mat(int m, vector w, int[] v, int[] u, matrix d) {
      matrix[m, cols(d)] logX;
      for (j in 1:cols(d)) {
        for (i in 1:m) {
            logX[i, j] = log_sum_exp(w[u[i]:(u[i+1]-1)] + d[v[u[i]:(u[i+1]-1)], j]);
        }
      }
      return logX;
    }
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
    real<lower=0> aveStDMetaPriorExpect;
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
    matrix[NMicrobeNodes, NMicrobeTips] logMicrobeAncestors = log(microbeAncestors[, 1:NMicrobeTips]);
    matrix[NSamples, 1 + NEffects + NHostNodes] fullModelMat = append_col(rep_vector(1, NSamples),
                                                                append_col(modelMat,
                                                                 hostAncestorsExpanded));
    vector[rows(csr_extract_w(fullModelMat))] wFullModelMat = csr_extract_w(fullModelMat);
    int vFullModelMat[size(csr_extract_v(fullModelMat))] = csr_extract_v(fullModelMat);
    int uFullModelMat[size(csr_extract_u(fullModelMat))] = csr_extract_u(fullModelMat);
    vector[rows(csr_extract_w(hostAncestors))] wHostAncestors = csr_extract_w(hostAncestors[1:NHostTips, ]);
    vector[rows(wHostAncestors)] wLogHostAncestors = log(wHostAncestors);
    int vHostAncestors[size(csr_extract_v(hostAncestors[1:NHostTips, ]))] = csr_extract_v(hostAncestors[1:NHostTips, ]);
    int uHostAncestors[size(csr_extract_u(hostAncestors[1:NHostTips, ]))] = csr_extract_u(hostAncestors[1:NHostTips, ]);
    vector[rows(csr_extract_w(microbeAncestors[, 1:NMicrobeTips]'))] wMicrobeAncestorsT = csr_extract_w(microbeAncestors[, 1:NMicrobeTips]');
    vector[rows(wMicrobeAncestorsT)] wLogMicrobeAncestorsT = log(wMicrobeAncestorsT);
    int vMicrobeAncestorsT[size(csr_extract_v(microbeAncestors[, 1:NMicrobeTips]'))] = csr_extract_v(microbeAncestors[, 1:NMicrobeTips]');
    int uMicrobeAncestorsT[size(csr_extract_u(microbeAncestors[, 1:NMicrobeTips]'))] = csr_extract_u(microbeAncestors[, 1:NMicrobeTips]');
}
parameters {
    real<lower=0> aveStD;
    simplex[2 * NFactors + 3] stDProps;
    vector[NTimeBins - 1] timeBinMetaVar;
    real<lower=0> aveStDMeta;
    simplex[3] metaVarProps;
    simplex[2] hostMetaVarProps;
    row_vector[NMicrobeNodes] phyloLogVarMultPrev;
    vector[NHostNodes] phyloLogVarMultADiv;
    matrix[NHostNodes, NMicrobeNodes] phyloLogVarMultRaw;
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
              * microbeAncestors
          + logMicrobeEdges;
    microbeScales
        = exp((logMicrobeVarRaw'
              - mean_of_log_prod_exp_csr_mat_vect(NMicrobeTips,
                                                  wLogMicrobeAncestorsT,
                                                  vMicrobeAncestorsT,
                                                  uMicrobeAncestorsT,
                                                  logMicrobeVarRaw'))
              * 0.5)';
    logRelativeEvolRates
        = append_row(0,
            cumulative_sum(timeBinMetaVar)
                * metaScales[2]
                * sqrt(hostMetaVarProps[1]));
    relativeEvolRates
        = exp(logRelativeEvolRates);
    logHostVarRaw
        = hostAncestors
              * phyloLogVarMultADiv
                * metaScales[2]
                * sqrt(hostMetaVarProps[2])
          + log_prod_exp_mat_vect(logEdgeToBin,
                                  logRelativeEvolRates);
    hostScales
        = exp(log(scales[2 * NFactors + 1])
              + (logHostVarRaw
                 - mean_of_log_prod_exp_csr_mat_vect(NHostTips,
                                                     wLogHostAncestors,
                                                     vHostAncestors,
                                                     uHostAncestors,
                                                     logHostVarRaw))
                 * 0.5);
    logPhyloScalesRaw
        = hostAncestors
              * (phyloLogVarMultRaw * metaScales[3])
              * microbeAncestors
          + log_outer_exp(logHostVarRaw, logMicrobeVarRaw);
    phyloScales
        = exp(log(scales[2 * NFactors + 2])
              + (logPhyloScalesRaw
                 - log_sum_exp(log_prod_exp_mats(log_prod_exp_csr_mat_dense_mat(NHostTips,
                                                                                wLogHostAncestors,
                                                                                vHostAncestors,
                                                                                uHostAncestors,
                                                                                logPhyloScalesRaw),
                                                 logMicrobeAncestors))
                 + log(NHostTips * NMicrobeTips))
                * 0.5); // should be able to optimize, maybe by precalculating the kronecker of host and microbe ancestry and using a single csr x vect multiplication
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
    aveStDMeta ~ exponential(1.0 / aveStDMetaPriorExpect);
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
                                      wMicrobeAncestorsT,
                                      vMicrobeAncestorsT,
                                      uMicrobeAncestorsT,
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
