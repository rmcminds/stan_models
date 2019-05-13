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
    int NObs;
    int NMicrobeNodes;
    int NMicrobeTips;
    int NFactors;
    int NSubPerFactor[NFactors];
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
    simplex[2] codivVsCophyMetaVarProps[3];
    simplex[2] codivVsCophyVarProps[3];
    real<lower=0> stDLogitMicrobe;
    real<lower=0> stDLogitHost;
    vector[NMicrobeNodes - NMicrobeTips] phyloLogitVarMicrobe;
    vector[NHostNodes - NHostTips] phyloLogitVarHost;
    row_vector[NMicrobeNodes] phyloLogVarMultPrev;
    row_vector[NMicrobeNodes] phyloLogVarDivPrev;
    vector[NHostNodes] phyloLogVarMultADiv;
    vector[NHostNodes] phyloLogVarDivADiv;
    matrix[NHostNodes, NMicrobeNodes] phyloLogVarMultRaw;
    matrix[NHostNodes, NMicrobeNodes] phyloLogVarCodivRaw;
    matrix[NEffects + NHostNodes + 1, NMicrobeNodes + 1] rawMicrobeNodeEffects;
}
transformed parameters {
    simplex[2 * NSubfactors + 3] subfactProps;
    vector<lower=0>[2 * NSubfactors + 3] scales;
    vector<lower=0>[3] metaScales;
    matrix[NMicrobeNodes, 2] newMicrobeNHs;
    row_vector<lower=0>[NMicrobeNodes] newMicrobeEdges;
    matrix[NHostNodes, 2] newHostNHs;
    vector<lower=0>[NHostNodes] newHostEdges;
    row_vector<lower=0>[NMicrobeNodes] microbeOUExpectedVariance;
    row_vector[NMicrobeNodes] microbeRateShifts;
    row_vector<lower=0>[NMicrobeNodes] microbeRates;
    row_vector[NMicrobeNodes] microbeDivergenceVariance;
    row_vector<lower=0>[NMicrobeNodes] microbeDivergence;
    row_vector<lower=0>[NMicrobeNodes] microbeScales;
    vector<lower=0>[NHostNodes] hostOUExpectedVariance;
    vector[NHostNodes] hostRateShifts;
    vector<lower=0>[NHostNodes] hostRates;
    vector[NHostNodes] hostDivergenceVariance;
    vector<lower=0>[NHostNodes] hostDivergence;
    vector<lower=0>[NHostNodes] hostScales;
    matrix<lower=0>[NHostNodes, NMicrobeNodes] cophyloExpectedVariance;
    matrix[NHostNodes, NMicrobeNodes] cophyloRateShifts;
    matrix<lower=0>[NHostNodes, NMicrobeNodes] cophyloRates;
    matrix[NHostNodes, NMicrobeNodes] coDivergenceVariance;
    matrix<lower=0>[NHostNodes, NMicrobeNodes] coDivergence;
    matrix<lower=0>[NHostNodes, NMicrobeNodes] coPhyloScales;
    matrix[NEffects + NHostNodes + 1, NMicrobeNodes + 1] scaledMicrobeNodeEffects;
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
    newMicrobeEdges
        = (newMicrobeNHs[,2] - newMicrobeNHs[,1])';
    microbeOUExpectedVariance
        = rescaleOU(newMicrobeNHs, microbeOUAlpha)';
    microbeRateShifts
        = sqrt(newMicrobeEdges)
          .* phyloLogVarMultPrev
             * metaScales[1]
               * codivVsCophyMetaVarProps[1,1];
    microbeRates
        = microbeOUExpectedVariance
          .* exp(microbeRateShifts
                 * microbeAncestorsT);
    microbeRates
        = codivVsCophyVarProps[1,1]
          * microbeRates
          / mean(microbeRates * microbeTipAncestorsT[2:,]);
    microbeDivergenceVariance
        = sqrt(newMicrobeEdges)
          .* phyloLogVarDivPrev
             * metaScales[1]
               * codivVsCophyMetaVarProps[1,2];
    microbeDivergence
        = exp(microbeDivergenceVariance
              * microbeAncestorsT);
    microbeDivergence
        = codivVsCophyVarProps[1,2]
          * microbeDivergence
          / mean(microbeDivergence * microbeTipAncestorsT[2:,]);
    microbeScales
        = sqrt(microbeRates + microbeDivergence);
    newHostNHs
        = makeNHMat(hostParents,
                    NHostTips,
                    hostLogitNH + stDLogitHost
                                  * stDLogitHostPriorExpect
                                  * phyloLogitVarHost);
    newHostEdges
        = newHostNHs[,2] - newHostNHs[,1];
    hostOUExpectedVariance
        = rescaleOU(newHostNHs, hostOUAlpha);
    hostRateShifts
        = sqrt(newHostEdges)
          .* phyloLogVarMultADiv
             * metaScales[2]
               * codivVsCophyMetaVarProps[2,1];
    hostRates
        = hostOUExpectedVariance
          .* exp(hostAncestors
                 * hostRateShifts);
    hostRates
        = codivVsCophyVarProps[2,1]
          * hostRates
          / mean(hostTipAncestors * hostRates);
    hostDivergenceVariance
        = sqrt(newHostEdges)
          .* phyloLogVarDivADiv
             * metaScales[2]
               * codivVsCophyMetaVarProps[2,2];
    hostDivergence
        = exp(hostAncestors
              * hostDivergenceVariance);
    hostDivergence
        = codivVsCophyVarProps[2,2]
          * hostDivergence
          / mean(hostTipAncestors * hostDivergence);
    hostScales
        = scales[2 * NSubfactors + 1]
          * sqrt(hostRates + hostDivergence);
    cophyloExpectedVariance
        = hostRates
          * microbeRates;
    cophyloRateShifts
        = phyloLogVarMultRaw
          * metaScales[3]
            * codivVsCophyMetaVarProps[3,1];
    cophyloRates
        = cophyloExpectedVariance
          .* exp(hostAncestors
                 * cophyloRateShifts
                   * microbeAncestorsT);
    cophyloRates
        = codivVsCophyVarProps[3,1]
          * cophyloRates
          / mean(hostTipAncestors
                 * cophyloRates
                   * microbeTipAncestorsT[2:,]);
    coDivergenceVariance
        = phyloLogVarCodivRaw
          * metaScales[3]
            * codivVsCophyMetaVarProps[3,2];
    coDivergence
        = exp(hostAncestors
              * coDivergenceVariance
                * microbeAncestorsT);
    coDivergence
        = codivVsCophyVarProps[3,2]
          * coDivergence
          / mean(hostTipAncestors
                 * coDivergence
                   * microbeTipAncestorsT[2:,]);
    coPhyloScales
        = scales[2 * NSubfactors + 2]
          * sqrt(cophyloRates + coDivergence);
    scaledMicrobeNodeEffects
        = append_col(
                append_row(1.0,
                           append_row(subfactLevelMat * segment(scales, 1, NSubfactors),
                                      hostScales)),
                append_row(
                    append_row(
                        scales[2 * NSubfactors + 3],
                        subfactLevelMat * segment(scales, NSubfactors + 1, NSubfactors))
                    * microbeScales,
                    coPhyloScales))
          .* rawMicrobeNodeEffects;
}
model {
    matrix[NSamples, NMicrobeTips] sampleTipEffects;
    vector[NObs] logit_ratios;
    aveStD ~ exponential(1.0 / aveStDPriorExpect);
    target += dirichSubFact_lpdf;
    stDProps ~ dirichlet(rep_vector(1, 2 * NFactors + 3));
    for (i in 1:3)
        codivVsCophyVarProps[i] ~ dirichlet(rep_vector(1, 2));
    hostOUAlpha ~ exponential(1.0 / hostOUAlphaPriorExpect);
    microbeOUAlpha ~ exponential(1.0 / microbeOUAlphaPriorExpect);
    aveStDMeta ~ exponential(1.0);
    metaVarProps ~ dirichlet(rep_vector(1, 3));
    for (i in 1:3)
        codivVsCophyMetaVarProps[i] ~ dirichlet(rep_vector(1, 2));
    stDLogitMicrobe ~ exponential(1.0);
    stDLogitHost ~ exponential(1.0);
    phyloLogitVarMicrobe ~ normal(0,1);
    phyloLogitVarHost ~ normal(0,1);
    phyloLogVarMultPrev ~ normal(0,1);
    phyloLogVarMultADiv ~ normal(0,1);
    to_vector(phyloLogVarMultRaw) ~ normal(0,1);
    phyloLogVarDivPrev ~ normal(0,1);
    phyloLogVarDivADiv ~ normal(0,1);
    to_vector(phyloLogVarCodivRaw) ~ normal(0,1);
    to_vector(rawMicrobeNodeEffects)[2:] ~ normal(0,1);
    rawMicrobeNodeEffects[1,1] ~ logistic(0,1);
    to_vector(baseLevelMat * rawMicrobeNodeEffects[2:(NEffects + 1),]) ~ normal(0,1);
    sampleTipEffects = modelMat * (scaledMicrobeNodeEffects * microbeTipAncestorsT);
    for (n in 1:NObs)
        logit_ratios[n] = sampleTipEffects[sampleNames[n], microbeTipNames[n]];
    present ~ bernoulli_logit(logit_ratios);
}
generated quantities {
    matrix[NSumTo0, NMicrobeNodes + 1] baseLevelEffects;
    baseLevelEffects
        = baseLevelMat * scaledMicrobeNodeEffects[2:(NEffects + 1),];
}
