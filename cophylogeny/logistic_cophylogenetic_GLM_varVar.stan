data {
    int NEstSamples;
    int NObs;
    int NEstTaxNodes;
    int NEstTips;
    int NFactors;
    int NEffects;
    int NHostNodes;
    int NHostTips;
    int NTimeBins;
    int present[NObs];
    int samplenames[NObs];
    int tipnames[NObs];
    real<lower=0> taxAveStDPriorExpect;
    vector[NTimeBins] timeBinSizes;
    matrix[NEstTaxNodes, NEstTaxNodes] ancestors;
    matrix[NEffects, NFactors] factLevelMat;
    matrix[NEstSamples, NEffects] modelMat;
    matrix[NHostNodes, NTimeBins] edgetobin;
    matrix[NHostNodes, NHostNodes] hostAncestors;
    matrix[NEstSamples, NHostNodes] hostAncestorsExpanded;
    row_vector[NEstTaxNodes] taxEdges;
}
parameters {
    real<lower=0> aveStD;
    simplex[2 * NFactors + 2] stDProps;
    simplex[NTimeBins] timeBinProps;
    row_vector[NEstTaxNodes - NEstTips] phyloLogVarMultPrev;
    vector[NHostNodes - NHostTips] phyloLogVarMultADiv;
    matrix[NHostNodes - NHostTips, NEstTaxNodes - NEstTips] phyloLogVarMultRaw;
    real globalIntercept;
    vector[NEffects + NHostNodes] rawAlphaDivEffects;
    vector[NEstTaxNodes] rawPrevalenceEffects;
    matrix[NEffects + NHostNodes, NEstTaxNodes] rawTaxNodeEffects;
}
transformed parameters {
    vector<lower=0>[2 * NFactors + 2] scales;
    row_vector<lower=0>[NEstTaxNodes] taxVarRaw;
    row_vector<lower=0>[NEstTaxNodes] taxScales;
    vector<lower=0>[NHostNodes] hostVarRaw;
    vector<lower=0>[NHostNodes] hostScales;
    matrix<lower=0>[NHostNodes, NEstTaxNodes] phyloScales;
    vector[NEffects + NHostNodes] scaledAlphaDivEffects;
    matrix[NEffects + NHostNodes, NEstTaxNodes] scaledTaxNodeEffects;
    matrix[NEffects + NHostNodes, NEstTips] scaledTipEffects;
    scales
        = sqrt((2 * NFactors + 2) * stDProps)
          * aveStD;
    taxVarRaw
        = exp(phyloLogVarMultPrev
              * ancestors[(NEstTips + 1):, ])
          .* taxEdges;
    taxScales
        = sqrt(taxVarRaw
               / mean(taxVarRaw * ancestors[, 1:NEstTips]));
    hostVarRaw
        = exp(hostAncestors[, (NHostTips + 1):]
              * phyloLogVarMultADiv)
          .* (edgetobin * timeBinProps);
    hostScales
        = scales[2 * NFactors + 1]
          * sqrt(hostVarRaw
                 / mean(hostAncestors[1:NHostTips, ] * hostVarRaw));
    phyloScales
        = exp(hostAncestors[, (NHostTips + 1):]
              * phyloLogVarMultRaw
              * ancestors[(NEstTips + 1):, ])
          .* (hostVarRaw * taxVarRaw);
    phyloScales
        = scales[2 * NFactors + 2]
          * sqrt(phyloScales
                 / mean(hostAncestors[1:NHostTips, ]
                        * phyloScales
                        * ancestors[, 1:NEstTips]));
    scaledAlphaDivEffects
        = append_row(globalIntercept,
            append_row(factLevelMat * segment(scales, 1, NFactors),
                hostScales)
          .* rawAlphaDivEffects);
    scaledTaxNodeEffects
        = append_row(factLevelMat * segment(scales, NFactors + 1, NFactors) * taxScales,
                phyloScales)
          .* rawTaxNodeEffects;
    scaledTipEffects
        = scaledTaxNodeEffects * ancestors[, 1:NEstTips];
}
model {
    matrix[NEstSamples, NEstTips] sampleTipEffects;
    vector[NObs] logit_ratios;
    aveStD ~ exponential(1.0 / taxAveStDPriorExpect);
    stDProps ~ dirichlet(rep_vector(1, 2 * NFactors + 2));
    timeBinProps ~ dirichlet(NTimeBins * timeBinSizes);
    phyloLogVarMultPrev ~ normal(0,1);
    phyloLogVarMultADiv ~ normal(0,1);
    to_vector(phyloLogVarMultRaw) ~ normal(0,1);
    globalIntercept ~ normal(0, 50);
    rawAlphaDivEffects ~ normal(0,1);
    to_vector(rawTaxNodeEffects) ~ normal(0,1);
    sampleTipEffects
        = append_col(modelMat, hostAncestorsExpanded)
          * (rep_matrix(scaledAlphaDivEffects, NEstTips) + scaledTipEffects);
    for (n in 1:NObs)
        logit_ratios[n] = sampleTipEffects[samplenames[n], tipnames[n]];
    present ~ bernoulli_logit(logit_ratios);
}
generated quantities {
    vector[NTimeBins] relativeEvolRates = timeBinProps ./ timeBinSizes;
}
