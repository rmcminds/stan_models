functions {
    vector rescaleOU(matrix nhs, real alpha) {
        return (exp(-2.0 * alpha * (1 - nhs[,2]))
                - exp(-2.0 * alpha * (1 - nhs[,1]))); //the history of evolution is erased as it gets farther away in time
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
    real<lower=0> hostOUAlphaPriorExpect; //expectation for the strength of the ornstein-uhlenbeck effect in hosts
    real<lower=0> microbeOUAlphaPriorExpect; //expectation for the strength of the ornstein-uhlenbeck effect in microbes
    matrix[NEffects, NFactors] factLevelMat;
    matrix[NSamples, NEffects + NHostNodes + 1] modelMat;
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
    real<lower=0> hostOUAlpha; //the strength of the ornstein-uhlenbeck effect in hosts
    real<lower=0> microbeOUAlpha; //the strength of the ornstein-uhlenbeck effect in microbes
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
        = rescaleOU(microbeNH, microbeOUAlpha)' //rescaled edges are a drop-in replacement for the times between divergences. **All factors are affected by a single parameter**
          .* exp((phyloLogVarMultPrev
                  * metaScales[1])
                 * microbeAncestorsT);
    microbeScales
        = sqrt(microbeVarRaw
               / mean(microbeVarRaw * microbeTipAncestorsT[2:,]));
    hostVarRaw
        = rescaleOU(hostNH, hostOUAlpha) //rescaled edges are a drop-in replacement
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
    target += exponential_lpdf(hostOUAlpha | 1.0 / hostOUAlphaPriorExpect); //the strength of the ornstein-uhlenbeck attraction is shrunk toward 0
    target += exponential_lpdf(microbeOUAlpha | 1.0 / microbeOUAlphaPriorExpect);
    target += exponential_lpdf(aveStDMeta | 1.0);
    target += dirichlet_lpdf(metaVarProps | rep_vector(1, 3));
    target += std_normal_lpdf(phyloLogVarMultPrev);
    target += std_normal_lpdf(phyloLogVarMultADiv);
    target += std_normal_lpdf(to_vector(phyloLogVarMultRaw));
    target += std_normal_lpdf(to_vector(rawMicrobeNodeEffects)[2:]);
    target += logistic_lpdf(rawMicrobeNodeEffects[1,1] | 0,1);
    sampleTipEffects = modelMat * (scaledMicrobeNodeEffects * microbeTipAncestorsT);
    for (n in 1:NObs)
        logit_ratios[n] = sampleTipEffects[sampleNames[n], microbeTipNames[n]];
    target += bernoulli_logit_lpmf(present | logit_ratios);
}
