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
    matrix[NSamples, 1 + NEffects + NHostNodes] fullModelMat = append_col(rep_vector(1, NSamples),
                                                                append_col(modelMat,
                                                                 hostAncestorsExpanded));
    vector[rows(csr_extract_w(fullModelMat))] wFullModelMat = csr_extract_w(fullModelMat);
    int vFullModelMat[size(csr_extract_v(fullModelMat))] = csr_extract_v(fullModelMat);
    int uFullModelMat[size(csr_extract_u(fullModelMat))] = csr_extract_u(fullModelMat);

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
    vector[NTimeBins] relativeEvolRates;
    row_vector<lower=0>[NMicrobeNodes] microbeVarRaw;
    row_vector<lower=0>[NMicrobeNodes] microbeScales;
    vector<lower=0>[NHostNodes] hostVarRaw;
    vector<lower=0>[NHostNodes] hostScales;
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
        = exp(metaScales[1] * phyloLogVarMultPrev
              * microbeAncestors[(NMicrobeTips + 1):, ])
          .* microbeEdges;
    microbeScales
        = sqrt(microbeVarRaw
               / mean(microbeVarRaw * microbeAncestors[, 1:NMicrobeTips]));
    relativeEvolRates
        = append_row(1,
            exp(cumulative_sum(timeBinMetaVar)
                * metaScales[2]
                * sqrt(hostMetaVarProps[1])));
    hostVarRaw
        = exp(hostAncestors[, (NHostTips + 1):]
              * phyloLogVarMultADiv
                * metaScales[2]
                * sqrt(hostMetaVarProps[2]))
          .* (edgeToBin
              * relativeEvolRates);
    hostScales
        = scales[2 * NFactors + 1]
          * sqrt(hostVarRaw
                 / mean(hostAncestors[1:NHostTips, ] * hostVarRaw));
    phyloScales
        = exp(hostAncestors[, (NHostTips + 1):]
              * (phyloLogVarMultRaw * metaScales[3])
              * microbeAncestors[(NMicrobeTips + 1):, ])
          .* (hostVarRaw * microbeVarRaw);
    phyloScales
        = scales[2 * NFactors + 2]
          * sqrt(phyloScales
                 / mean(hostAncestors[1:NHostTips, ]
                        * phyloScales
                        * microbeAncestors[, 1:NMicrobeTips]));
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
    aveStDMeta ~ exponential(1.0);
    metaVarProps ~ dirichlet(rep_vector(1, 3));
    hostMetaVarProps ~ dirichlet(rep_vector(1, 2));
    phyloLogVarMultPrev ~ normal(0,1);
    phyloLogVarMultADiv ~ normal(0,1);
    to_vector(phyloLogVarMultRaw) ~ normal(0,1);
    globalIntercept ~ normal(0, 50);
    rawAlphaDivEffects ~ normal(0,1);
    to_vector(rawMicrobeNodeEffects) ~ normal(0,1);

    microbePreMult = (microbeAncestors[, 1:NMicrobeTips]' * scaledMicrobeNodeEffects')';

    for (t in 1:NMicrobeTips) {
        sampleTipEffects[, t]
            = csr_matrix_times_vector(NSamples,
                                      1 + NEffects + NHostNodes,
                                      wFullModelMat,
                                      vFullModelMat,
                                      uFullModelMat,
                                      scaledAlphaDivEffects
                                       + scaledMicrobeNodeEffects * microbeAncestors[, t]);
    }
    for (n in 1:NObs)
        logit_ratios[n] = sampleTipEffects[sampleNames[n], microbeTipNames[n]];
    present ~ bernoulli_logit(logit_ratios);
}
