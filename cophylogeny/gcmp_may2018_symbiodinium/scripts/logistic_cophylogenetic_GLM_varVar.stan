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
    int NEffects;
    int NHostNodes;
    int NHostTips;
    int NLatent;
    int present[NObs];
    int sampleNames[NObs];
    int microbeTipNames[NObs];
    real<lower=0> aveStDPriorExpect;
    real<lower=0> aveStDMetaPriorExpect;
    real<lower=0> hostOUAlphaPriorExpect;
    real<lower=0> microbeOUAlphaPriorExpect;
    real<lower=0> stDLogitMicrobePriorExpect;
    real<lower=0> stDLogitHostPriorExpect;
    matrix[NEffects, NFactors] factLevelMat;
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
parameters {
    real<lower=0> aveStD;
    vector<lower=0>[2 * NFactors + NLatent + 3] stDPropsRaw;
    positive_ordered[NLatent] latentStDPropsRaw;
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
    matrix<lower=0>[NSamples, NLatent] latentVars;
    matrix[NEffects + NHostNodes + NLatent + 1, NMicrobeNodes + 1] rawMicrobeNodeEffects;
}
transformed parameters {
    simplex[2 * (NFactors + NLatent) + 3] stDProps;
    vector<lower=0>[2 * (NFactors + NLatent) + 3] scales;
    vector<lower=0>[3] metaScales;
    matrix[NMicrobeNodes, 2] newMicrobeNHs;
    matrix[NHostNodes, 2] newHostNHs;
    row_vector<lower=0>[NMicrobeNodes] microbeVarRaw;
    row_vector<lower=0>[NMicrobeNodes] microbeScales;
    vector<lower=0>[NHostNodes] hostVarRaw;
    vector<lower=0>[NHostNodes] hostScales;
    matrix<lower=0>[NHostNodes, NMicrobeNodes] phyloVarRaw;
    matrix<lower=0>[NHostNodes, NMicrobeNodes] phyloScales;
    matrix[NEffects + NHostNodes + NLatent + 1, NMicrobeNodes + 1] scaledMicrobeNodeEffects;
    matrix[NSamples, NMicrobeTips] sampleTipEffects;
    stDProps
        = append_row(
              append_row(stDPropsRaw[1:NLatent],
                         latentStDPropsRaw),
              stDPropsRaw[(NLatent + 1):(2 * NFactors + NLatent + 3)]);
    stDProps
        = stDProps / sum(stDProps);
    scales
        = sqrt((2 * (NFactors + NLatent) + 3) * stDProps)
          * aveStD;
    metaScales
        = sqrt(3 * metaVarProps)
          * aveStDMeta
          * aveStDMetaPriorExpect;
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
        = scales[2 * (NFactors + NLatent) + 1]
          * sqrt(hostVarRaw
                 / mean(hostTipAncestors * hostVarRaw));

    phyloVarRaw
        = exp(hostAncestors
              * (phyloLogVarMultRaw
                 * metaScales[3])
              * microbeAncestorsT)
          .* (hostVarRaw * microbeVarRaw);
    phyloScales
        = scales[2 * (NFactors + NLatent) + 2]
          * sqrt(phyloVarRaw
                 / mean(hostTipAncestors
                        * (phyloVarRaw
                           * microbeTipAncestorsT[2:,])));
    scaledMicrobeNodeEffects
        = append_col(
              append_row(segment(scales, 1, NLatent),
                  append_row(globalScale,
                      append_row(factLevelMat * segment(scales, 2 * NLatent + 1, NFactors),
                          hostScales))),
              append_row(
                append_row(segment(scales, NLatent + 1, NLatent),
                    append_row(scales[2 * (NFactors + NLatent) + 3],
                        factLevelMat * segment(scales, NFactors + 2 * NLatent + 1, NFactors)))
                * microbeScales,
                phyloScales))
          .* rawMicrobeNodeEffects;
    sampleTipEffects
        = append_col(latentVars, modelMat) * (scaledMicrobeNodeEffects * microbeTipAncestorsT);
}
model {
    vector[NObs] logit_ratios;
    aveStD ~ exponential(1.0 / aveStDPriorExpect);
    stDPropsRaw ~ exponential(1.0);
    latentStDPropsRaw ~ exponential(1.0);
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
    to_vector(latentVars) ~ normal(0,2);
    to_vector(rawMicrobeNodeEffects) ~ normal(0,1);
    to_vector(baseLevelMat * rawMicrobeNodeEffects[(2 * NLatent + 2):(2 * NLatent + NEffects + 1),]) ~ normal(0,1);
    for (n in 1:NObs)
        logit_ratios[n] = sampleTipEffects[sampleNames[n], microbeTipNames[n]];
    present ~ bernoulli_logit(logit_ratios);
}
generated quantities {
    matrix[NSumTo0, NMicrobeNodes + 1] baseLevelEffects;
    int present_pred[NObs];
    baseLevelEffects
        = baseLevelMat * scaledMicrobeNodeEffects[(2 * NLatent + 2):(2 * NLatent + NEffects + 1),];
    for (n in 1:NObs)
        present_pred[n] = bernoulli_logit_rng(sampleTipEffects[sampleNames[n], microbeTipNames[n]]);
}
