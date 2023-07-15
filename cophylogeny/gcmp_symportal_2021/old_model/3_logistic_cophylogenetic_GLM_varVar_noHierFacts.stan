functions {
    vector rescaleOU(matrix nhs, real alpha) {
        return (exp(-2.0 * alpha * (1 - nhs[,2]))
                - exp(-2.0 * alpha * (1 - nhs[,1])));
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
    int present[NObs];
    int sampleNames[NObs];
    int microbeTipNames[NObs];
    real<lower=0> aveStDPriorExpect;
    real<lower=0> aveStDMetaPriorExpect;
    real<lower=0> hostOUAlphaPriorExpect;
    real<lower=0> microbeOUAlphaPriorExpect;
    matrix[NEffects, NFactors] factLevelMat;
    matrix[NSamples, NEffects + NHostNodes + 1] modelMat; //matrix will be used the same in the model but is not necessarily binary
    int NSumTo0;
    matrix[NSumTo0, NEffects] baseLevelMat; //categorical factors are better modeled with sum-to-zero contrasts, where there is one fewer parameter than the nubmer of levels in the factor. The expectation for the missing level is the sum of all the other effects, and this matrix is used to calculate it
    matrix[NMicrobeNodes, NMicrobeNodes] microbeAncestorsT;
    matrix[NMicrobeNodes + 1, NMicrobeTips] microbeTipAncestorsT;
    matrix[NHostNodes, NHostNodes] hostAncestors;
    matrix[NHostTips, NHostNodes] hostTipAncestors;
    matrix[NHostNodes, 2] hostNH;
    matrix[NMicrobeNodes, 2] microbeNH;
}
parameters {
    real<lower=0> aveStD;
    simplex[2 * NFactors + 3] stDProps;
    real<lower=0> hostOUAlpha;
    real<lower=0> microbeOUAlpha;
    real<lower=0> aveStDMeta;
    simplex[3] metaVarProps;
    row_vector[NMicrobeNodes] phyloLogVarMultPrev;
    vector[NHostNodes] phyloLogVarMultADiv;
    matrix[NHostNodes, NMicrobeNodes] phyloLogVarMultRaw;
    matrix[NEffects + NHostNodes + 1, NMicrobeNodes + 1] rawMicrobeNodeEffects;
}
transformed parameters {
    vector<lower=0>[2 * NFactors + 3] scales;
    vector<lower=0>[3] metaScales;
    row_vector<lower=0>[NMicrobeNodes] microbeVarRaw;
    row_vector<lower=0>[NMicrobeNodes] microbeScales;
    vector<lower=0>[NHostNodes] hostVarRaw;
    vector<lower=0>[NHostNodes] hostScales;
    matrix<lower=0>[NHostNodes, NMicrobeNodes] phyloVarRaw;
    matrix<lower=0>[NHostNodes, NMicrobeNodes] phyloScales;
    matrix[NEffects + NHostNodes + 1, NMicrobeNodes + 1] scaledMicrobeNodeEffects;
    scales
        = sqrt((2 * NFactors + 3) * stDProps)
          * aveStD;
    metaScales
        = sqrt(3 * metaVarProps)
          * aveStDMetaPriorExpect
          * aveStDMeta;
    microbeVarRaw
        = rescaleOU(microbeNH, microbeOUAlpha)'
          .* exp((phyloLogVarMultPrev
                  * metaScales[1])
                 * microbeAncestorsT);
    microbeScales
        = sqrt(microbeVarRaw
               / mean(microbeVarRaw * microbeTipAncestorsT[2:,]));
    hostVarRaw
        = rescaleOU(hostNH, hostOUAlpha)
          .* exp(hostAncestors
              * (phyloLogVarMultADiv
                 * metaScales[2]));
    hostScales
        = scales[2 * NFactors + 1]
          * sqrt(hostVarRaw
                 / mean(hostTipAncestors * hostVarRaw));
    phyloVarRaw
        = exp(hostAncestors
              * (phyloLogVarMultRaw
                 * metaScales[3])
              * microbeAncestorsT)
          .* (hostVarRaw * microbeVarRaw);
    phyloScales
        = scales[2 * NFactors + 2]
          * sqrt(phyloVarRaw
                 / mean(hostTipAncestors
                        * (phyloVarRaw
                           * microbeTipAncestorsT[2:,])));
    scaledMicrobeNodeEffects
        = append_col(
                append_row(1.0,
                           append_row(factLevelMat * segment(scales, 1, NFactors),
                                      hostScales)),
                append_row(
                    append_row(
                        scales[2 * NFactors + 3],
                        factLevelMat * segment(scales, NFactors + 1, NFactors))
                    * microbeScales,
                    phyloScales))
          .* rawMicrobeNodeEffects;
}
model {
    matrix[NSamples, NMicrobeTips] sampleTipEffects;
    vector[NObs] logit_ratios;
    target += exponential_lpdf(aveStD | 1.0 / aveStDPriorExpect);
    target += dirichlet_lpdf(stDProps | rep_vector(1, 2 * NFactors + 3));
    target += exponential_lpdf(hostOUAlpha | 1.0 / hostOUAlphaPriorExpect);
    target += exponential_lpdf(microbeOUAlpha | 1.0 / microbeOUAlphaPriorExpect);
    target += exponential_lpdf(aveStDMeta | 1.0);
    target += dirichlet_lpdf(metaVarProps | rep_vector(1, 3));
    target += std_normal_lpdf(phyloLogVarMultPrev);
    target += std_normal_lpdf(phyloLogVarMultADiv);
    target += std_normal_lpdf(to_vector(phyloLogVarMultRaw));
    target += std_normal_lpdf(to_vector(rawMicrobeNodeEffects)[2:]);
    target += logistic_lpdf(rawMicrobeNodeEffects[1,1] | 0,1);
    target += std_normal_lpdf(to_vector(baseLevelMat * rawMicrobeNodeEffects[2:(NEffects + 1),])); //for symmetrical prior expectations, shrinkage needs to occur on the 'missing' level the same as for the others. This causes the marginal variance of all effects to be lower than expected; an adjustment for this is included in the model matrix prior to model fitting
    sampleTipEffects = modelMat * (scaledMicrobeNodeEffects * microbeTipAncestorsT);
    for (n in 1:NObs)
        logit_ratios[n] = sampleTipEffects[sampleNames[n], microbeTipNames[n]];
    target += bernoulli_logit_lpmf(present | logit_ratios);
}
generated quantities {
    matrix[NSumTo0, NMicrobeNodes + 1] baseLevelEffects;
    baseLevelEffects
        = baseLevelMat * scaledMicrobeNodeEffects[2:(NEffects + 1),];
}
