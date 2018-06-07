data {
    int NSamples;
    int NObs;
    int NMicrobeFactors;
    int NMicrobeEffects;
    int NMicrobeNodes;
    int NMicrobeTips;
    int NSampleFactors;
    int NSampleEffects;
    int NHostNodes;
    int NHostTips;
    int NTimeBins;
    int present[NObs];
    int sampleNames[NObs];
    int microbeTipNames[NObs];
    real<lower=0> aveStDPriorExpect;
    matrix[NMicrobeNodes, NMicrobeNodes] microbeAncestors;
    matrix[NMicrobeEffects, NMicrobeFactors] sampleFactLevelMat;
    matrix[NMicrobeTips, NMicrobeEffects] microbeModelMat;
    matrix[NSampleEffects, NSampleFactors] sampleFactLevelMat;
    matrix[NSamples, NSampleEffects] sampleModelMat;
    vector[NTimeBins] timeBinSizes;
    matrix[NHostNodes, NTimeBins] edgeToBin;
    matrix[NHostNodes, NHostNodes] hostAncestors;
    matrix[NSamples, NHostNodes] hostAncestorsExpanded;
    row_vector[NMicrobeNodes] microbeEdges;
}
transformed data {
    int NVarPar = (NMicrobeFactors + 1) * (NSampleFactors + 1) + 3;
    // alpha div and prevalence effects, specificity effects, and 3 phylogeny effects
}
parameters {
    real<lower=0> aveStD;
    simplex[NVarPar] stDProps;
    simplex[NTimeBins] timeBinProps;
    row_vector[NMicrobeNodes - NMicrobeTips] phyloLogVarMultPrev;
    vector[NHostNodes - NHostTips] phyloLogVarMultADiv;
    matrix[NHostNodes - NHostTips, NMicrobeNodes - NMicrobeTips] phyloLogVarMultRaw;
    real globalIntercept;
    vector[NSampleEffects + NHostNodes] rawAlphaDivEffects;
    matrix[NSampleEffects + NHostNodes + 1, NMicrobeEffects + NMicrobeNodes] rawMicrobeNodeEffects;
}
transformed parameters {
    vector<lower=0>[NVarPar] scales;
    row_vector<lower=0>[NMicrobeNodes] microbeVarRaw;
    row_vector<lower=0>[NMicrobeNodes] microbeScales;
    vector<lower=0>[NHostNodes] hostVarRaw;
    vector<lower=0>[NHostNodes] hostScales;
    matrix<lower=0>[NHostNodes, NMicrobeNodes] phyloScales;
    vector[NSampleEffects + NHostNodes + 1] scaledAlphaDivEffects;
    matrix[NSampleEffects + NHostNodes + 1, NMicrobeEffects + NMicrobeNodes] scaledMicrobeNodeEffects;
    scales
        = sqrt(NVarPar * stDProps)
          * aveStD;
    microbeVarRaw
        = exp(phyloLogVarMultPrev
              * microbeAncestors[(NMicrobeTips + 1):, ])
          .* microbeEdges;
    microbeScales
        = sqrt(microbeVarRaw
               / mean(microbeVarRaw * microbeAncestors[, 1:NMicrobeTips]));
    hostVarRaw
        = exp(hostAncestors[, (NHostTips + 1):]
              * phyloLogVarMultADiv)
          .* (edgeToBin * timeBinProps);
    hostScales
        = scales[2 * NSampleFactors + 1]
          * sqrt(hostVarRaw
                 / mean(hostAncestors[1:NHostTips, ] * hostVarRaw));
    phyloScales
        = exp(hostAncestors[, (NHostTips + 1):]
              * phyloLogVarMultRaw
              * microbeAncestors[(NMicrobeTips + 1):, ])
          .* (hostVarRaw * microbeVarRaw);
    phyloScales
        = scales[2 * NSampleFactors + 2]
          * sqrt(phyloScales
                 / mean(hostAncestors[1:NHostTips, ]
                        * phyloScales
                        * microbeAncestors[, 1:NMicrobeTips]));
    scaledAlphaDivEffects
        = append_row(globalIntercept,
            append_row(sampleFactLevelMat * segment(scales, 1, NSampleFactors),
                hostScales)
          .* rawAlphaDivEffects);
    scaledMicrobeNodeEffects
        = append_row(
            append_row(
                scales[2 * NSampleFactors + 3],
                sampleFactLevelMat * segment(scales, NSampleFactors + 1, NSampleFactors))
            * microbeScales,
            phyloScales)
          .* rawMicrobeNodeEffects;
}
model {
    matrix[NSamples, NMicrobeTips] sampleTipEffects;
    vector[NObs] logit_ratios;
    aveStD ~ exponential(1.0 / aveStDPriorExpect);
    stDProps ~ dirichlet(rep_vector(1, NVarPar));
    timeBinProps ~ dirichlet(NTimeBins * timeBinSizes);
    phyloLogVarMultPrev ~ normal(0,1);
    phyloLogVarMultADiv ~ normal(0,1);
    to_vector(phyloLogVarMultRaw) ~ normal(0,1);
    globalIntercept ~ normal(0, 50);
    rawAlphaDivEffects ~ normal(0,1);
    to_vector(rawMicrobeNodeEffects) ~ normal(0,1);
    sampleTipEffects
        = append_col(rep_vector(1, NSamples),
            append_col(sampleModelMat,
                hostAncestorsExpanded))
          * (rep_matrix(scaledAlphaDivEffects, NMicrobeTips)
             + scaledMicrobeNodeEffects * microbeAncestors[, 1:NMicrobeTips]);
    for (n in 1:NObs)
        logit_ratios[n] = sampleTipEffects[sampleNames[n], microbeTipNames[n]];
    present ~ bernoulli_logit(logit_ratios);
}
generated quantities {
    vector[NTimeBins] relativeEvolRates = timeBinProps ./ timeBinSizes;
}
