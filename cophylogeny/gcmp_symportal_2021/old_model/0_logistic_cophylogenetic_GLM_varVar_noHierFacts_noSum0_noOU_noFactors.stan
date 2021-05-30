data {
    int NSamples;
    int NObs; //usually the number of samples * the number of microbial taxa (but should be possible to use subset if for some reason the data aren't conclusive about the presence of particular microbes in some samples)
    int NMicrobeNodes;
    int NMicrobeTips;
    int NHostNodes;
    int NHostTips;
    int present[NObs]; //binary data, whether a given microbial taxon existed in a given sample
    int sampleNames[NObs]; //maps each data point to the relevant sample
    int microbeTipNames[NObs]; //maps each data point to the relevant microbial taxon
    real<lower=0> aveStDPriorExpect; //prior expectation for the rate of evolution
    real<lower=0> aveStDMetaPriorExpect; //prior expectation for the rate of evolution of the rate of evolution
    matrix[NMicrobeNodes, NMicrobeNodes] microbeAncestorsT; //binary matrix for microbial tree topology (transposed)
    matrix[NMicrobeNodes + 1, NMicrobeTips] microbeTipAncestorsT; //subset of microbeAncestorsT plus a row for the root
    matrix[NHostNodes, NHostNodes] hostAncestors; //binary matrix for host tree topology
    matrix[NHostTips, NHostNodes] hostTipAncestors; //subset of hostAncestors
    matrix[NSamples, NHostNodes + 1] modelMat; //binary matrix mapping samples to the nodes of the host tree, plus a column for the root
    vector[NHostNodes] hostEdges; //host tree edge lengths
    row_vector[NMicrobeNodes] microbeEdges; //microbial tree edge lengths
}
parameters {
    real<lower=0> aveStD; //mean rate of evolution of associations
    simplex[3] stDProps; //variance partitions for sigma2[h], sigma2[m], and sigma2[pm]
    real<lower=0> aveStDMeta; //mean expectated scale of the evolution associated with a divergence event of the rate of evolution
    simplex[3] metaVarProps; //variance partitions for metaSigmas
    row_vector[NMicrobeNodes] phyloLogVarMultPrev; //unscaled edge-specific shifts in the rate of evolution of overal microbe prevalence/generalism
    vector[NHostNodes] phyloLogVarMultADiv; //unscaled edge-specific shifts in the rate of evolution of overal host alpha diversity/generalism
    matrix[NHostNodes, NMicrobeNodes] phyloLogVarMultRaw; //unscaled shifts in the rate of evolution of associations between particular microbial clades with particular host clades
    matrix[NHostNodes + 1, NMicrobeNodes + 1] rawMicrobeNodeEffects; //unscaled evolution of the probability of association between particular microbial clades and particular host clades
}
transformed parameters {
    vector<lower=0>[3] scales; //mean rate of evolution of host generalism, microbial generalism, and host-microbe specificity (sigma[h], sigma[m], sigma[pm])
    vector<lower=0>[3] metaScales; //mean expectated scale of the evolution associated with a divergence event of the rate of evolution of generalism and specificity
    row_vector<lower=0>[NMicrobeNodes] microbeVarRaw; //unnormalized variance of the contrast between the expected overall prevalence of each microbial node and its parent
    row_vector<lower=0>[NMicrobeNodes] microbeScales; //variance normalized by the mean root-to-tip distance and converted to standard deviation
    vector<lower=0>[NHostNodes] hostVarRaw; //unnormalized variance of the contrast between the expected overall alpha diversity of each host node and its parent
    vector<lower=0>[NHostNodes] hostScales; //variance normalized by the mean root-to-tip distance and converted to standard deviation
    matrix<lower=0>[NHostNodes, NMicrobeNodes] phyloVarRaw; //unnormalized variance of the contrast between the expected specificity of each host node / microbe combination and the specificity of  their parents
    matrix<lower=0>[NHostNodes, NMicrobeNodes] phyloScales; //variance normalized by the mean root-to-tip distance and converted to standard deviation
    matrix[NHostNodes + 1, NMicrobeNodes + 1] scaledMicrobeNodeEffects; //scaled evolution of the probability of association between particular microbial clades and particular host clades
    scales
        = sqrt(3 * stDProps)
          * aveStDPriorExpect
          * aveStD;
    metaScales
        = sqrt(3 * metaVarProps)
          * aveStDMetaPriorExpect
          * aveStDMeta;
    microbeVarRaw
        = microbeEdges
          .* exp((phyloLogVarMultPrev
                  * metaScales[1])
                 * microbeAncestorsT); //time between divergences multiplied by the rate of evolution
    microbeScales
        = sqrt(microbeVarRaw
               / mean(microbeVarRaw * microbeTipAncestorsT[2:,])); //normalize by mean root-to-tip distance
    hostVarRaw
        = hostEdges
          .* exp(hostAncestors
              * (phyloLogVarMultADiv
                 * metaScales[2])); //time between divergences multiplied by the rate of evolution
    hostScales
        = scales[1]
          * sqrt(hostVarRaw
                 / mean(hostTipAncestors * hostVarRaw)); //normalize by mean root-to-tip distance
    phyloVarRaw
        = exp(hostAncestors
              * (phyloLogVarMultRaw
                 * metaScales[3])
              * microbeAncestorsT)
          .* (hostVarRaw * microbeVarRaw); //time between divergences multiplied by the rate of evolution
    phyloScales
        = scales[2]
          * sqrt(phyloVarRaw
                 / mean(hostTipAncestors
                        * (phyloVarRaw
                           * microbeTipAncestorsT[2:,]))); //normalize by mean root-to-tip distance
    scaledMicrobeNodeEffects
        = append_col(
                append_row(1.0,
                           hostScales),
                append_row(scales[3] * microbeScales,
                           phyloScales))
          .* rawMicrobeNodeEffects; //scale evolutionary effects
}
model {
    matrix[NSamples, NMicrobeTips] sampleTipEffects;
    vector[NObs] logit_ratios;
    target += exponential_lpdf(aveStD | 1.0); //mean rate of evolution is shrunk toward 0
    target += dirichlet_lpdf(stDProps | rep_vector(1, 3)); //no shrinkage for variance partitioning
    target += exponential_lpdf(aveStDMeta | 1.0); //mean scale of rate evolution is shrunk toward 0
    target += dirichlet_lpdf(metaVarProps | rep_vector(1, 3)); //no shrinkage for meta-variance partitioning
    target += std_normal_lpdf(phyloLogVarMultPrev); //rate evolution is shrunk toward 0
    target += std_normal_lpdf(phyloLogVarMultADiv); //rate evolution is shrunk toward 0
    target += std_normal_lpdf(to_vector(phyloLogVarMultRaw)); //rate evolution is shrunk toward 0
    target += std_normal_lpdf(to_vector(rawMicrobeNodeEffects)[2:]); //evolution of associations is shrunk toward 0
    target += logistic_lpdf(rawMicrobeNodeEffects[1,1] | 0,1); //no shrinkage for overall mean probability of association
    sampleTipEffects = modelMat * (scaledMicrobeNodeEffects * microbeTipAncestorsT); //sum effects of all higher nodes to get expectation for each combination of sample and microbe
    for (n in 1:NObs)
        logit_ratios[n] = sampleTipEffects[sampleNames[n], microbeTipNames[n]]; //map expectations to observations
    target += bernoulli_logit_lpmf(present | logit_ratios); //likelihood of data given parameters
}
