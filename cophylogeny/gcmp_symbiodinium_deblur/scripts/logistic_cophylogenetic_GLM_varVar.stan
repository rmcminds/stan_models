functions {
    vector rescaleOU(matrix nhs, real alpha) {
        return (exp(-2.0 * alpha * (1 - nhs[,2]))
                - exp(-2.0 * alpha * (1 - nhs[,1])));
    }
}
data {
    int NSamples;
    int NObs;
    int NMicrobeFactors;
    int NMicrobeEffects;
    int NMicrobeSumTo0;
    int NMicrobeNodes;
    int NMicrobeTips;
    int NHostFactors;
    int NHostEffects;
    int NHostSumTo0;
    int NHostNodes;
    int NHostTips;
    int present[NObs];
    int sampleNames[NObs];
    int microbeTipNames[NObs];
    real<lower=0> aveStDPriorExpect;
    real<lower=0> aveStDMetaPriorExpect;
    real<lower=0> hostOUAlphaPriorExpect;
    real<lower=0> microbeOUAlphaPriorExpect;
    matrix[NHostEffects, NHostFactors] hostFactLevelMat;
    matrix[NSamples, NHostEffects + NHostNodes + 1] hostModelMat;
    matrix[NHostSumTo0, NHostEffects] hostBaseLevelMat;
    matrix[NMicrobeEffects, NMicrobeFactors] microbeFactLevelMat;
    matrix[NMicrobeNodes + NMicrobeEffects + 1, NMicrobeTips] microbeModelMat;
    matrix[NMicrobeSumTo0, NMicrobeEffects] hostBaseLevelMat;
    matrix[NMicrobeNodes, NMicrobeNodes] microbeAncestorsT;
    matrix[NMicrobeNodes, NMicrobeTips] microbeTipAncestorsT;
    matrix[NHostNodes, NHostNodes] hostAncestors;
    matrix[NHostTips, NHostNodes] hostTipAncestors;
    matrix[NHostNodes, 2] hostNodeHeights;
    matrix[NMicrobeNodes, 2] microbeNodeHeights;
    real<lower=0> globalScale;
}
parameters {
    real<lower=0> aveStD;
    simplex[4 * NHostFactors * NMicrobeFactors - 1] stDProps;
    real<lower=0> hostOUAlpha;
    real<lower=0> microbeOUAlpha;
    real<lower=0> aveStDMeta;
    simplex[3] metaVarProps;
    row_vector[NMicrobeNodes] phyloLogVarMultPrev;
    vector[NHostNodes] phyloLogVarMultADiv;
    matrix[NHostNodes, NMicrobeNodes] phyloLogVarMultRaw;
    matrix[NHostEffects + NHostNodes + 1, NMicrobeNodes + NMicrobeEffects + 1] rawMicrobeNodeEffects;
}
transformed parameters {
    vector<lower=0>[4 * NHostFactors * NMicrobeFactors - 1] scales;
    vector<lower=0>[3] metaScales;
    row_vector<lower=0>[NMicrobeNodes] microbeVarRaw;
    row_vector<lower=0>[NMicrobeNodes] microbeScales;
    vector<lower=0>[NHostNodes] hostVarRaw;
    vector<lower=0>[NHostNodes] hostScales;
    matrix<lower=0>[NHostNodes, NMicrobeNodes] phyloVarRaw;
    matrix<lower=0>[NHostNodes, NMicrobeNodes] phyloScales;
    matrix[NHostEffects + NHostNodes + 1, NMicrobeNodes + NMicrobeEffects + 1] scaledMicrobeNodeEffects;
    scales
        = sqrt((4 * NHostFactors * NMicrobeFactors - 1) * stDProps)
          * aveStD;
    metaScales
        = sqrt(3 * metaVarProps)
          * aveStDMeta;
    microbeVarRaw
        = rescaleOU(microbeNodeHeights, microbeOUAlpha)'
          .* exp((phyloLogVarMultPrev
                  * metaScales[1])
                 * microbeAncestorsT);
    microbeScales
        = sqrt(microbeVarRaw
               / mean(microbeVarRaw * microbeTipAncestorsT));
    hostVarRaw
        = rescaleOU(hostNodeHeights, hostOUAlpha)
          .* exp(hostAncestors
              * (phyloLogVarMultADiv
                 * metaScales[2]));
    hostScales
        = scales[2 * NHostFactors + 1]
          * sqrt(hostVarRaw
                 / mean(hostTipAncestors * hostVarRaw));

    phyloVarRaw
        = exp(hostAncestors
              * (phyloLogVarMultRaw
                 * metaScales[3])
              * microbeAncestorsT)
          .* (hostVarRaw * microbeVarRaw);
    phyloScales
        = scales[2 * NHostFactors + 2]
          * sqrt(phyloVarRaw
                 / mean(hostTipAncestors
                        * (phyloVarRaw
                           * microbeTipAncestorsT)));
    scaledMicrobeNodeEffects
        = append_col(
                append_row(globalScale,
                           append_row(hostFactLevelMat * segment(scales, 1, NHostFactors),
                                      hostScales)),
                append_row(
                    append_row(
                        scales[2 * NHostFactors + 3],
                        hostFactLevelMat * segment(scales, NHostFactors + 1, NHostFactors))
                    * microbeScales,
                    phyloScales))
          .* rawMicrobeNodeEffects;
}
model {
    matrix[NSamples, NMicrobeTips] sampleTipEffects;
    vector[NObs] logit_ratios;
    aveStD ~ exponential(1.0 / aveStDPriorExpect);
    stDProps ~ dirichlet(rep_vector(1, 2 * NHostFactors + 3));
    hostOUAlpha ~ exponential(1.0 / hostOUAlphaPriorExpect);
    microbeOUAlpha ~ exponential(1.0 / microbeOUAlphaPriorExpect);
    aveStDMeta ~ exponential(1.0 / aveStDMetaPriorExpect);
    metaVarProps ~ dirichlet(rep_vector(1, 3));
    phyloLogVarMultPrev ~ normal(0,1);
    phyloLogVarMultADiv ~ normal(0,1);
    to_vector(phyloLogVarMultRaw) ~ normal(0,1);
    to_vector(rawMicrobeNodeEffects) ~ normal(0,1);
    to_vector(hostBaseLevelMat * rawMicrobeNodeEffects[2:(NHostEffects + 1),]) ~ normal(0,1);
    sampleTipEffects = hostModelMat * scaledMicrobeNodeEffects * microbeModelMat;
    for (n in 1:NObs)
        logit_ratios[n] = sampleTipEffects[sampleNames[n], microbeTipNames[n]];
    present ~ bernoulli_logit(logit_ratios);
}
generated quantities {
    matrix[NHostSumTo0, NMicrobeNodes + 1] hostBaseLevelEffects
        = hostBaseLevelMat * scaledMicrobeNodeEffects[2:(NHostEffects + 1),];
}
