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
    matrix[NEffects, sum(NSubPerFactor)] subfactLevelMat;
    matrix[NSamples, NEffects + NHostNodes + 1] modelMat;
    int NSumTo0;
    matrix[NSumTo0, NEffects] baseLevelMat;
    matrix[NMicrobeNodes, NMicrobeNodes] microbeAncestorsT;
    matrix[NMicrobeNodes + 1, NMicrobeTips] microbeTipAncestorsT;
    matrix[NHostNodes, NHostNodes] hostAncestors;
    matrix[NHostTips, NHostNodes] hostTipAncestors;
    vector[NHostNodes] hostEdges;
    row_vector[NMicrobeNodes] microbeEdges;
}
transformed data {
    row_vector[NHostNodes] PDescHost
        = rep_row_vector(1.0 / NHostTips, NHostTips)
          * hostTipAncestors;
    vector[NHostNodes] logPDescHost
        = log(PDescHost)';
    vector[NMicrobeNodes] PDescMicrobe
        = microbeTipAncestorsT[2:,]
          * rep_vector(1.0 / NMicrobeTips, NMicrobeTips);
    row_vector[NMicrobeNodes] logPDescMicrobe
        = log(PDescMicrobe)';
    matrix[NHostNodes, NMicrobeNodes] logPDescBoth
        = rep_matrix(logPDescMicrobe, NHostNodes) + rep_matrix(logPDescHost, NMicrobeNodes);
    matrix[NMicrobeNodes, NMicrobeNodes] microbeAncestorsTCont
        = microbeAncestorsT;
    matrix[NHostNodes, NHostNodes] hostAncestorsCont
        = hostAncestors;
    vector[NHostNodes] logHostEdges = log(hostEdges);
    row_vector[NMicrobeNodes] logMicrobeEdges = log(microbeEdges);
    vector[NHostNodes] sqrtHostEdges = sqrt(hostEdges);
    row_vector[NMicrobeNodes] sqrtMicrobeEdges = sqrt(microbeEdges);
    matrix[NHostNodes, NMicrobeNodes] outerEdges
        = sqrt(hostEdges
               * microbeEdges);
    int NSubfactorGammas = 0;
    int NSubfactors = sum(NSubPerFactor);
    for(i in 1:NMicrobeNodes) {
        microbeAncestorsTCont[i,i] = 1.0 - exp(-1);
    }
    for(i in 1:NHostNodes) {
        hostAncestorsCont[i,i] = 1.0 - exp(-1);
    }
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
    real<lower=0> aveStDMeta;
    simplex[NFactors + 3] metaVarProps;
    vector<lower=0>[NSubfactorGammas] subfactMetaPropsRaw;
    row_vector[NMicrobeNodes] phyloLogVarMultPrev;
    vector[NHostNodes] phyloLogVarMultADiv;
    matrix[NHostNodes, NMicrobeNodes] phyloLogVarMultRaw;
    matrix[NSubfactors, NMicrobeNodes] phyloLogVarMultFacts;
    matrix[NHostNodes, NMicrobeNodes] rawHostMicrobeSpecificity;
    matrix[NEffects, NMicrobeNodes] rawMicrobeEffectSpecificity;
    vector[NHostNodes] rawLatentMicrobeHostSpecificity;
    vector[NEffects] rawLatentMicrobeEffectSpecificity;
    row_vector[NMicrobeNodes] rawMicrobeLatentHostSpecificity;
    row_vector[NMicrobeNodes] rawMicrobeLatentIndividualSpecificity
    real rawLatentHostADiv;
    vector[NHostNodes] rawHostADiv;
    real rawLatentMicrobePrevalence;
    row_vector[NMicrobeNodes] rawMicrobePrevalence;
    vector[NEffects] rawEffectsADiv;
    real rawLatentPhyEffect;
    real rawLatentIndivMicrobeEffect;
    vector<lower=-1, upper=1>[3] varEffectCor;
    real<lower=0> hostLatent1;
    vector[NHostNodes - 1] hostLatent;
    real<lower=0> individualLatent1;
    vector[NIndividuals - 1] individualLatent;
    real<lower=0> microbeLatent1;
    row_vector[NMicrobeNodes - 1] microbeLatent;
}
transformed parameters {
    simplex[2 * NSubfactors + 3] subfactProps;
    vector<lower=0>[2 * NSubfactors + 3] scales;
    simplex[NSubfactors + 3] subfactMetaProps;
    vector<lower=0>[NSubfactors + 3] metaScales;
    row_vector<lower=0>[NMicrobeNodes] microbeScales;
    vector<lower=0>[NHostNodes] hostScales;
    matrix<lower=0>[NSubfactors, NMicrobeNodes] factScales;
    matrix[NEffects + NHostNodes + 3, NMicrobeNodes + 2] scaledMicrobeNodeEffects;
    vector<lower=0, upper=1>[3] varEffectChol2 = sqrt(1 - square(varEffectCor));
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
                dirichSubFact_lpdf += lmultiply(-NSubPerFactor[i], sum_gamma)
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
                dirichSubFact_lpdf += lmultiply(-NSubPerFactor[i], sum_gamma)
                                      + dirichlet_lpdf(subfactProps[(NSubfactors + normStart):(NSubfactors + normStart - 1 + NSubPerFactor[i])] | rep_vector(1, NSubPerFactor[i]));
                subfactProps[(NSubfactors + normStart):(NSubfactors + normStart - 1 + NSubPerFactor[i])]
                    = subfactProps[(NSubfactors + normStart):(NSubfactors + normStart - 1 + NSubPerFactor[i])]
                      * stDProps[NFactors + i];
                sum_gamma = 1 + sum(segment(subfactMetaPropsRaw,
                                            rawStart,
                                            NSubPerFactor[i] - 1));
                subfactMetaProps[normStart:(normStart - 1 + NSubPerFactor[i])]
                    = append_row(1, segment(subfactMetaPropsRaw, rawStart, NSubPerFactor[i] - 1))
                      / sum_gamma;
                dirichSubFact_lpdf += lmultiply(-NSubPerFactor[i], sum_gamma)
                                      + dirichlet_lpdf(subfactMetaProps[normStart:(normStart - 1 + NSubPerFactor[i])] | rep_vector(1, NSubPerFactor[i]));
                subfactMetaProps[normStart:(normStart - 1 + NSubPerFactor[i])]
                    = subfactMetaProps[normStart:(normStart - 1 + NSubPerFactor[i])]
                      * metaVarProps[i];
                rawStart += NSubPerFactor[i] - 1;
            } else {
                subfactProps[normStart]
                    = stDProps[i];
                subfactProps[NSubfactors + normStart]
                    = stDProps[NFactors + i];
                subfactMetaProps[normStart]
                    = metaVarProps[i];
            }
            normStart += NSubPerFactor[i];
        }
    }
    subfactProps[(2 * NSubfactors + 1):(2 * NSubfactors + 3)]
        = stDProps[(2 * NFactors + 1):(2 * NFactors + 3)];
    subfactMetaProps[(NSubfactors + 1):(NSubfactors + 3)]
        = metaVarProps[(NFactors + 1):(NFactors + 3)];
    scales
        = sqrt((NHostNodes + NSubfactors + 1) * (NMicrobeNodes + 1) * subfactProps)
          * aveStD;
    metaScales
        = sqrt(((NHostNodes + 1) * (NMicrobeNodes + 1) + NSubfactors * NMicrobeNodes) * subfactMetaProps)
          * aveStDMeta;
    {
        row_vector[NMicrobeNodes] logMicrobeVarRaw;
        vector[NHostNodes] logHostVarRaw;
        matrix[NHostNodes, NMicrobeNodes] logPhyloVarRaw;
        matrix[NSubfactors, NMicrobeNodes] logFactVarRaw;
        matrix[NHostNodes, NMicrobeNodes] phyloScales;
        logMicrobeVarRaw
            = logMicrobeEdges
              + (sqrtMicrobeEdges
                 .* phyloLogVarMultPrev
                 * metaScales[NSubfactors + 1])
                * microbeAncestorsTCont;
        logHostVarRaw
            = logHostEdges
              + hostAncestorsCont
                * (sqrtHostEdges
                   .* phyloLogVarMultADiv
                   * metaScales[NSubfactors + 2]);
        logPhyloVarRaw
            = hostAncestorsCont
                  * (outerEdges
                     .* phyloLogVarMultRaw
                     * metaScales[NSubfactors + 3])
                  * microbeAncestorsTCont
              + rep_matrix(logHostVarRaw, NMicrobeNodes)
              + rep_matrix(logMicrobeVarRaw, NHostNodes);
        logFactVarRaw
            = metaScales[1:NSubfactors]
              * sqrtMicrobeEdges
              .* phyloLogVarMultFacts
              * microbeAncestorsTCont
              + rep_matrix(logMicrobeVarRaw, NSubfactors);
        microbeScales
            = scales[2 * NSubfactors + 3]
              * exp((logMicrobeVarRaw
                     - log_sum_exp(logPDescMicrobe + logMicrobeVarRaw))
                    * 0.5);
        hostScales
            = scales[2 * NSubfactors + 1]
              * exp((logHostVarRaw
                     - log_sum_exp(logPDescHost + logHostVarRaw))
                    * 0.5);
        phyloScales
            = scales[2 * NSubfactors + 2]
              * exp((logPhyloVarRaw
                     - log_sum_exp(logPDescBoth + logPhyloVarRaw))
                    * 0.5);
        for(f in 1:NSubfactors) {
            factScales[f,]
                = scales[NSubfactors + f]
                  * exp((logFactVarRaw[f,]
                        - log_sum_exp(logPDescMicrobe + logFactVarRaw[f,]))
                        * 0.5);
        }
        scaledMicrobeNodeEffects[1,1]
            = rawMicrobeNodeEffects[1,1]; // intercept
        scaledMicrobeNodeEffects[1, 2:(NMicrobeNodes + 1)]
            = microbeScales .* (varEffectChol2[1] * rawMicrobePrevalence
                                + varEffectCor[1] * phyloLogVarMultPrev); //microbe prevalence
        scaledMicrobeNodeEffects[1, NMicrobeNodes + 2]
            = rawLatentMicrobePrevalence; // effect on microbe prevalence of the latent microbe factor

        scaledMicrobeNodeEffects[2:(NEffects + 1), 1]
            = subfactLevelMat
              * segment(scales, 1, NSubfactors)
              .* rawEffectsADiv; // samplewise factor alpha diversity
        scaledMicrobeNodeEffects[2:(NEffects + 1), 2:(NMicrobeNodes + 1)]
            = subfactLevelMat
              * factScales
              .* rawMicrobeEffectSpecificity; // samplewise factor microbe interactions
        scaledMicrobeNodeEffects[2:(NEffects + 1), NMicrobeNodes + 2]
            = rawLatentMicrobeEffectSpecificity; // interaction of samplewise factors with latent microbe factor

        scaledMicrobeNodeEffects[(NEffects + 2):(NEffects + NHostNodes + 1), 1]
            = hostScales .* (varEffectChol2[2] * rawHostADiv
                             + varEffectCor[2] * phyloLogVarMultADiv); // host alpha diversity
        scaledMicrobeNodeEffects[(NEffects + 2):(NEffects + NHostNodes + 1), 2:(NMicrobeNodes + 1)]
            = phyloScales .* (varEffectChol2[3] * rawHostMicrobeSpecificity
                              + varEffectCor[3] * phyloLogVarMultRaw); // host microbe interactions
        scaledMicrobeNodeEffects[(NEffects + 2):(NEffects + NHostNodes + 1), NMicrobeNodes + 2]
            = rawLatentMicrobeHostSpecificity; // interaction of hosts with latent microbe factor

        scaledMicrobeNodeEffects[NEffects + NHostNodes + 2, 1]
            = rawLatentHostADiv; // effect on alpha diversity of the latent host factor
        scaledMicrobeNodeEffects[NEffects + NHostNodes + 2, 2:(NMicrobeNodes + 1)]
            = rawMicrobeLatentHostSpecificity; // interaction of microbes with latent host factor
        scaledMicrobeNodeEffects[NEffects + NHostNodes + 2, NMicrobeNodes + 2]
            = rawLatentPhyEffect; // interaction of latent host and microbe factors

        scaledMicrobeNodeEffects[NEffects + NHostNodes + 3, 1]
            = rawLatentIndividualADiv; // effect on alpha diversity of the latent colony factor
        scaledMicrobeNodeEffects[NEffects + NHostNodes + 3, 2:(NMicrobeNodes + 1)]
            = rawMicrobeLatentIndividualSpecificity; // interaction of microbes with the latent colony factor
        scaledMicrobeNodeEffects[NEffects + NHostNodes + 3, NMicrobeNodes + 2]
            = rawLatentIndivMicrobeEffect; // interaction of latent individual and microbe factors
    }
}
model {
    matrix[NSamples, NMicrobeTips] sampleTipEffects;
    vector[NObs] logit_ratios;
    target += exponential_lpdf(aveStD | 1.0 / aveStDPriorExpect);
    target += dirichSubFact_lpdf;
    target += dirichlet_lpdf(stDProps | rep_vector(1, 2 * NFactors + 3));
    target += exponential_lpdf(aveStDMeta | 1.0 / aveStDMetaPriorExpect);
    target += dirichlet_lpdf(metaVarProps | rep_vector(1, NFactors + 3));
    target += std_normal_lpdf(phyloLogVarMultPrev);
    target += std_normal_lpdf(phyloLogVarMultADiv);
    target += std_normal_lpdf(to_vector(phyloLogVarMultRaw));
    target += std_normal_lpdf(to_vector(phyloLogVarMultFacts));
    target += std_normal_lpdf(to_vector(rawMicrobeNodeEffects)[2:]);
    target += logistic_lpdf(rawMicrobeNodeEffects[1,1] | 0,1);
    target += std_normal_lpdf(to_vector(baseLevelMat * rawMicrobeNodeEffects[2:(NEffects + 1),]));
    sampleTipEffects = append_col(append_col(modelMat, append_row(hostLatent1, hostLatent)), append_row(individualLatent1, individualLatent)) * (scaledMicrobeNodeEffects * append_row(microbeTipAncestorsT, append_col(microbeLatent1, microbeLatent)));
    for (n in 1:NObs)
        logit_ratios[n] = sampleTipEffects[sampleNames[n], microbeTipNames[n]];
    target += bernoulli_logit_lpmf(present | logit_ratios);
}
generated quantities {
    row_vector[NMicrobeNodes] microbeNewEdges
        = square(microbeScales / scales[2 * NSubfactors + 3]);
    vector[NHostNodes] hostNewEdges
        = square(hostScales / scales[2 * NSubfactors + 1]);
    matrix[NSumTo0, NMicrobeNodes + 1] baseLevelEffects
        = baseLevelMat * scaledMicrobeNodeEffects[2:(NEffects + 1),];
}
