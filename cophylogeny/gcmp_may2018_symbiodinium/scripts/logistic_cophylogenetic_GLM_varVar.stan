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
    vector[nnzMicrobeTipAncestors] vMicrobeTipAncestors;
    vector[nuMicrobeTipAncestors] uMicrobeTipAncestors;
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
    matrix[NHostNodes, NTimeBins] logEdgeToBin = log(edgeToBin);
    vector[nnzHostAncestors] wHostAncestorsLog = log(wHostAncestors);
    vector[nnzHostTipAncestors] wHostTipAncestorsLog = log(wHostTipAncestors);
    vector[nnzHostKronMicrobeAncestors] wHostKronMicrobeAncestorsLog = log(wHostKronMicrobeAncestors);
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
    vector<lower=0>[NMicrobeNodes] microbeVarRaw;
    row_vector<lower=0>[NMicrobeNodes] microbeScales;
    vector[NTimeBins] logRelativeEvolRates;
    vector<lower=0>[NTimeBins] relativeEvolRates;
    vector<lower=0>[NHostNodes] hostVarRaw;
    vector<lower=0>[NHostNodes] hostScales;
    matrix<lower=0>[NHostNodes, NMicrobeNodes] phyloVarRaw;
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
                                  wMicrobeAncestors,
                                  vMicrobeAncestors,
                                  uMicrobeAncestors,
                                  phyloLogVarMultPrev
                                    * metaScales[1])
          .* microbeEdges;
    microbeScales
        = sqrt(microbeVarRaw
               / mean(csr_matrix_times_vector(NMicrobeTips,
                                              wMicrobeTipAncestors,
                                              vMicrobeTipAncestors,
                                              uMicrobeTipAncestors,
                                              microbeVarRaw)))';
    logRelativeEvolRates
        = append_row(0,
            cumulative_sum(timeBinMetaVar)
                * metaScales[2]
                * sqrt(hostMetaVarProps[1]));
    relativeEvolRates
        = exp(logRelativeEvolRates);
    hostVarRaw
        = csr_matrix_times_vector(NHostNodes,
                                  wHostAncestors,
                                  vHostAncestors,
                                  uHostAncestors,
                                  phyloLogVarMultADiv
                                    * metaScales[2]
                                    * sqrt(hostMetaVarProps[2]))
          .* (edgeToBin * relativeEvolRates);
    hostScales
        = scales[2 * NFactors + 1]
          * sqrt(hostVarRaw
                 / mean(csr_matrix_times_vector(NHostTips,
                                                wHostTipAncestors,
                                                vHostTipAncestors,
                                                uHostTipAncestors,
                                                hostVarRaw)));
    phyloVarRaw
        = to_matrix(csr_matrix_times_vector(NHostNodes * NMicrobeNodes,
                                            wHostKronMicrobeAncestors,
                                            vHostKronMicrobeAncestors,
                                            uHostKronMicrobeAncestors,
                                            phyloLogVarMultRaw
                                              * metaScales[3]),
                    NHostNodes,
                    NMicrobeNodes)
          .* (hostVarRaw * microbeVarRaw');
    phyloScales
        = scales[2 * NFactors + 2]
          * sqrt(phyloVarRaw
                 / mean(csr_matrix_times_vector(NHostTips * NMicrobeTips,
                                                wHostKronMicrobeTipAncestors,
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
