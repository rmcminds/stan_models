functions {
    matrix makeNHMat(matrix pm, int NTips, vector logitNH) {
        matrix[rows(pm), 2] newNHs = rep_matrix(0, rows(pm), 2);
        for (row in 1:(rows(pm) - NTips)) {
            newNHs[NTips + row, 1] = pm[NTips + row,] * newNHs[,2];
            newNHs[NTips + row, 2] = exp(log_inv_logit(logitNH[row])
                                         + log(1 - newNHs[NTips + row, 1]))
                                     + newNHs[NTips + row, 1];
        }
        for (row in 1:NTips) {
            newNHs[row, 1] = pm[row,] * newNHs[,2];
            newNHs[row, 2] = 1;
        }
        return newNHs;
    }
    vector rescaleOU(matrix nhs, real alpha) {
        return (exp(-2.0 * alpha * (1 - nhs[,2]))
                - exp(-2.0 * alpha * (1 - nhs[,1])));
    }
}
data {
    int NSamples;
    int NMicrobeNodes;
    int NMicrobeTips;
    int NFactors;
    int NSubPerFactor[NFactors];
    int NEffects;
    int NHostNodes;
    int NHostTips;
    int NSeqs;
    int present[NSamples, NSeqs];
    matrix[NMicrobeTips, NSeqs] profileMat;
    real<lower=0> aveStDPriorExpect;
    real<lower=0> aveStDMetaPriorExpect;
    real<lower=0> hostOUAlphaPriorExpect;
    real<lower=0> microbeOUAlphaPriorExpect;
    real<lower=0> stDLogitMicrobePriorExpect;
    real<lower=0> stDLogitHostPriorExpect;
    matrix[NEffects, sum(NSubPerFactor)] subfactLevelMat;
    matrix[NSamples, NEffects + NHostNodes + 1] modelMat;
    int NSumTo0;
    matrix[NSumTo0, NEffects] baseLevelMat;
    matrix[NMicrobeNodes, NMicrobeNodes] microbeAncestorsT;
    matrix[NMicrobeNodes + 1, NMicrobeTips] microbeTipAncestorsT;
    matrix[NHostNodes, NHostNodes] hostAncestors;
    matrix[NHostTips, NHostNodes] hostTipAncestors;
    matrix[NMicrobeNodes, NMicrobeNodes] microbeParents;
    matrix[NHostNodes, NHostNodes] hostParents;
    vector[NHostNodes - NHostTips] hostLogitNH;
    vector[NMicrobeNodes - NMicrobeTips] microbeLogitNH;
    real<lower=0> globalScale;
}
transformed data {
    int NSubfactorGammas = 0;
    int NSubfactors = sum(NSubPerFactor);
    for(i in 1:NFactors) {
        if(NSubPerFactor[i] > 1) {
            NSubfactorGammas += NSubPerFactor[i] - 1;
        }
    }
}
parameters {
    real<lower=0> aveStD;
    vector<lower=0>[2 * NSubfactorGammas] subfactPropsRaw;
    simplex[2 * NFactors + 3] stDProps;
    real<lower=0> hostOUAlpha;
    real<lower=0> microbeOUAlpha;
    real<lower=0> aveStDMeta;
    simplex[3] metaVarProps;
    real<lower=0> stDLogitMicrobe;
    real<lower=0> stDLogitHost;
    vector[NMicrobeNodes - NMicrobeTips] phyloLogitVarMicrobe;
    vector[NHostNodes - NHostTips] phyloLogitVarHost;
    row_vector[NMicrobeNodes] phyloLogVarMultPrev;
    vector[NHostNodes] phyloLogVarMultADiv;
    matrix[NHostNodes, NMicrobeNodes] phyloLogVarMultRaw;
    matrix[NEffects + NHostNodes + 1, NMicrobeNodes + 1] rawMicrobeNodeEffects;
}
transformed parameters {
    simplex[2 * NSubfactors + 3] subfactProps;
    vector<lower=0>[2 * NSubfactors + 3] scales;
    vector<lower=0>[3] metaScales;
    matrix[NMicrobeNodes, 2] newMicrobeNHs;
    matrix[NHostNodes, 2] newHostNHs;
    row_vector<lower=0>[NMicrobeNodes] microbeVarRaw;
    row_vector<lower=0>[NMicrobeNodes] microbeScales;
    vector<lower=0>[NHostNodes] hostVarRaw;
    vector<lower=0>[NHostNodes] hostScales;
    matrix<lower=0>[NHostNodes, NMicrobeNodes] phyloVarRaw;
    matrix<lower=0>[NHostNodes, NMicrobeNodes] phyloScales;
    matrix[NEffects + NHostNodes + 1, NMicrobeNodes + 1] scaledMicrobeNodeEffects;
    matrix[NSamples, NMicrobeTips] sampleTipEffects;
    real dirichSubFact_lpdf = 0;
    {
        int rawStart = 1;
        int normStart = 1;
        for (i in 1:NFactors) {
            if(NSubPerFactor[i] > 1) {
                real sum_gamma = 1 + sum(segment(subfactPropsRaw,
                                                 rawStart,
                                                 NSubPerFactor[i] - 1));
                subfactProps[normStart:(normStart - 1 + NSubPerFactor[i])]
                    = append_row(1, segment(subfactPropsRaw, rawStart, NSubPerFactor[i] - 1))
                      / sum_gamma;
                dirichSubFact_lpdf += -NSubPerFactor[i] * log(sum_gamma)
                                      + dirichlet_lpdf(subfactProps[normStart:(normStart - 1 + NSubPerFactor[i])] | rep_vector(1, NSubPerFactor[i]));
                subfactProps[normStart:(normStart - 1 + NSubPerFactor[i])]
                    = subfactProps[normStart:(normStart - 1 + NSubPerFactor[i])]
                      * stDProps[i];
                sum_gamma = 1 + sum(segment(subfactPropsRaw,
                                            NSubfactorGammas + rawStart,
                                            NSubPerFactor[i] - 1));
                subfactProps[(NSubfactors + normStart):(NSubfactors + normStart - 1 + NSubPerFactor[i])]
                    = append_row(1, segment(subfactPropsRaw, NSubfactorGammas + rawStart, NSubPerFactor[i] - 1))
                      / sum_gamma;
                dirichSubFact_lpdf += -NSubPerFactor[i] * log(sum_gamma)
                                      + dirichlet_lpdf(subfactProps[(NSubfactors + normStart):(NSubfactors + normStart - 1 + NSubPerFactor[i])] | rep_vector(1, NSubPerFactor[i]));
                subfactProps[(NSubfactors + normStart):(NSubfactors + normStart - 1 + NSubPerFactor[i])]
                    = subfactProps[(NSubfactors + normStart):(NSubfactors + normStart - 1 + NSubPerFactor[i])]
                      * stDProps[NFactors + i];
                rawStart += NSubPerFactor[i] - 1;
            } else {
                subfactProps[normStart]
                    = stDProps[i];
                subfactProps[NSubfactors + normStart]
                    = stDProps[NFactors + i];
            }
            normStart += NSubPerFactor[i];
        }
    }
    subfactProps[(2 * NSubfactors + 1):(2 * NSubfactors + 3)]
        = stDProps[(2 * NFactors + 1):(2 * NFactors + 3)];
    scales
        = sqrt((2 * NSubfactors + 3) * subfactProps)
          * aveStD;
    metaScales
        = sqrt(3 * metaVarProps)
          * aveStDMetaPriorExpect
          * aveStDMeta;
    newMicrobeNHs
        = makeNHMat(microbeParents,
                    NMicrobeTips,
                    microbeLogitNH + stDLogitMicrobe
                                     * stDLogitMicrobePriorExpect
                                     * phyloLogitVarMicrobe);
    microbeVarRaw
        = rescaleOU(newMicrobeNHs, microbeOUAlpha)'
          .* exp((phyloLogVarMultPrev
                  * metaScales[1])
                 * microbeAncestorsT);
    microbeScales
        = sqrt(microbeVarRaw
               / mean(microbeVarRaw * microbeTipAncestorsT[2:,]));
    newHostNHs
        = makeNHMat(hostParents,
                    NHostTips,
                    hostLogitNH + stDLogitHost
                                  * stDLogitHostPriorExpect
                                  * phyloLogitVarHost);
    hostVarRaw
        = rescaleOU(newHostNHs, hostOUAlpha)
          .* exp(hostAncestors
              * (phyloLogVarMultADiv
                 * metaScales[2]));
    hostScales
        = scales[2 * NSubfactors + 1]
          * sqrt(hostVarRaw
                 / mean(hostTipAncestors * hostVarRaw));

    phyloVarRaw
        = exp(hostAncestors
              * (phyloLogVarMultRaw
                 * metaScales[3])
              * microbeAncestorsT)
          .* (hostVarRaw * microbeVarRaw);
    phyloScales
        = scales[2 * NSubfactors + 2]
          * sqrt(phyloVarRaw
                 / mean(hostTipAncestors
                        * (phyloVarRaw
                           * microbeTipAncestorsT[2:,])));
    scaledMicrobeNodeEffects
        = append_col(
                append_row(globalScale,
                           append_row(subfactLevelMat * segment(scales, 1, NSubfactors),
                                      hostScales)),
                append_row(
                    append_row(
                        scales[2 * NSubfactors + 3],
                        subfactLevelMat * segment(scales, NSubfactors + 1, NSubfactors))
                    * microbeScales,
                    phyloScales))
          .* rawMicrobeNodeEffects;
    sampleTipEffects
        = modelMat
          * (scaledMicrobeNodeEffects
             * microbeTipAncestorsT);
}
model {
    aveStD ~ exponential(1.0 / aveStDPriorExpect);
    target += dirichSubFact_lpdf;
    stDProps ~ dirichlet(rep_vector(1, 2 * NFactors + 3));
    hostOUAlpha ~ exponential(1.0 / hostOUAlphaPriorExpect);
    microbeOUAlpha ~ exponential(1.0 / microbeOUAlphaPriorExpect);
    aveStDMeta ~ exponential(1.0);
    metaVarProps ~ dirichlet(rep_vector(1, 3));
    stDLogitMicrobe ~ exponential(1.0);
    stDLogitHost ~ exponential(1.0);
    phyloLogitVarMicrobe ~ normal(0,1);
    phyloLogitVarHost ~ normal(0,1);
    phyloLogVarMultPrev ~ normal(0,1);
    phyloLogVarMultADiv ~ normal(0,1);
    to_vector(phyloLogVarMultRaw) ~ normal(0,1);
    to_vector(rawMicrobeNodeEffects) ~ normal(0,1);
    to_vector(baseLevelMat * rawMicrobeNodeEffects[2:(NEffects + 1),]) ~ normal(0,1);
    to_array_1d(present) ~ bernoulli_logit(to_array_1d(to_array_2d(sampleTipEffects * profileMat)));
}
generated quantities {
    matrix[NSamples, NMicrobeTips] profilePresence;
    matrix[NSumTo0, NMicrobeNodes + 1] baseLevelEffects;
    for (i  in 1:NSamples) {
        for (j in 1:NMicrobeTips) {
            profilePresence[i,j] = bernoulli_logit_rng(sampleTipEffects[i,j]);
        }
    }
    baseLevelEffects
        = baseLevelMat * scaledMicrobeNodeEffects[2:(NEffects + 1),];
}
