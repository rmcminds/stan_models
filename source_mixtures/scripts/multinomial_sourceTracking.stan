data {
    int NSamples;
    int NMicrobeNodes;
    int NMicrobeTips;
    int NEnvs;
    int NFactors;
    int envs[NSamples];
    int NSepFacts;
    matrix[NEnvs, NSepFacts] envFact;
    matrix[NSamples, NFactors] modelMat;
    int y[NMicrobeTips, NSamples];
    real taxAveStDPriorExpect;
    real envPropAveStDPriorExpect;
    matrix[NMicrobeNodes, NMicrobeTips] microbeAncestors;
    vector[NMicrobeNodes] microbeNodeScales;
}
transformed data {
    real<lower=2.0> taxEffectsDF = 5.0;
    real<lower=2.0> envPropsDF = 5.0;
    int NEnvSamps = NEnvs * NSamples;
    int NTotalFactors = NSepFacts * NFactors;
    matrix[NEnvs + 1, NEnvs + 1] envTaxMat;
    matrix[NEnvs + 1, NSepFacts] envFactWithUnk;
    envTaxMat[1:NEnvs, 2:(NEnvs + 1)] = diag_matrix(rep_vector(1, NEnvs));
    envTaxMat[NEnvs + 1, 2:(NEnvs + 1)] = rep_row_vector(-1, NEnvs);
    envTaxMat[,1] = rep_vector(1, NEnvs + 1);
    envFactWithUnk = append_row(envFact, rep_row_vector(0, NSepFacts));
}
parameters {
    real<lower=0> taxAveStD;
    real<lower=0> envPropAveStD;
    simplex[NEnvs] envPropResidualProps;
    vector[NEnvs] basePropsRaw;
    matrix[NEnvs, NEnvs] envPropMeansRaw;
    matrix[NEnvs, NSamples] envPropResiduals;
    row_vector[NMicrobeNodes] taxIntercepts;
    simplex[NTotalFactors + 2] taxVarianceProps;
    simplex[NMicrobeTips] taxResidualVarianceProps;
    vector[(NEnvs + NTotalFactors + NSamples) * NMicrobeNodes] taxEffects;
    vector[NSamples] multinomialNuisance;
}
transformed parameters {
    vector<lower=0>[NEnvs] envPropResidualScales;
    matrix<lower=0>[NFactors, NSepFacts] taxFactScales;
    real<lower=0> taxEnvScale;
    vector<lower=0>[NMicrobeTips] taxResidualScales;
    matrix[NEnvs + 1, NEnvs] envPropMeans;
    vector[NEnvs + 1] log_samp_props[NSamples];
    matrix[NEnvs + 1, NMicrobeNodes] taxNode_normEnvEffects;
    matrix[NEnvs + 1, NMicrobeTips] tax_normEnvEffects;
    matrix[NFactors, NMicrobeNodes] taxNode_normFacts[NSepFacts];
    matrix[NFactors, NMicrobeTips] tax_normFacts[NSepFacts];
    vector[NEnvs + 1] baseProps;
    envPropResidualScales
        = sqrt((envPropsDF - 2)
               / envPropsDF
               * NEnvs
               * envPropResidualProps)
          * envPropAveStD;
    taxFactScales
        = to_matrix(
            sqrt((taxEffectsDF - 2)
                 / taxEffectsDF
                 * (NTotalFactors + 2)
                 * taxVarianceProps[1:NTotalFactors])
            * taxAveStD,
            NFactors, NSepFacts);
    taxEnvScale
        = sqrt((taxEffectsDF - 2)
               / taxEffectsDF
               * (NTotalFactors + 2)
               * taxVarianceProps[NTotalFactors + 1])
          * taxAveStD;
    taxResidualScales
        = sqrt((taxEffectsDF - 2)
               / taxEffectsDF
               * (NTotalFactors + 2)
               * taxVarianceProps[NTotalFactors + 2]
               * NMicrobeTips
               * taxResidualVarianceProps)
          * taxAveStD;
    baseProps = append_row(basePropsRaw, 0);
    for (i in 1:NEnvs) {
        envPropMeans[(1:(i - 1)]           = baseProps[(1:(i - 1)]
                                             - NEnvs * exp(envPropMeansRaw[, i]);
        envPropMeans[i, i]                 = baseProps[i];
        envPropMeans[(i + 1):(NEnvs + 1))] = baseProps[(i + 1):(NEnvs + 1))]
                                             - NEnvs * exp(envPropMeansRaw[, i]);
    }
    for (k in 1:NSamples) {
        vector[NEnvs + 1] samp_prop_normal;
        samp_prop_normal[(1:(envs[k] - 1)]           = envPropMeans[(1:(envs[k] - 1), envs[k]]
                                                       + envPropResidualScales[envs[k]] * envPropResiduals[, k];
        samp_prop_normal[envs[k]]                    = envPropMeans[envs[k], envs[k]];
        samp_prop_normal[(envs[k] + 1):(NEnvs + 1))] = envPropMeans[(envs[k] + 1):(NEnvs + 1)), envs[k]]
                                                       + envPropResidualScales[envs[k]] * envPropResiduals[, k];
        log_samp_props[k]                            = log_softmax(samp_prop_normal);
    }
    {
        taxEffectMat
            = to_matrix(taxEffects[1:(NEnvs * NMicrobeNodes)],
                        NEnvs, NMicrobeNodes);
        taxNode_normEnvEffects
            = diag_post_multiply(
                append_row(taxIntercepts,
                           taxEnvScale * taxEffectMat),
                microbeNodeScales);
    }
    tax_normEnvEffects = taxNode_normEnvEffects * microbeAncestors;
    for (h in 1:NSepFacts) {
        int startIndex = NMicrobeNodes * (NEnvs + NFactors * (h - 1)) + 1;
        matrix[NFactors, NMicrobeNodes] taxEffectMat;
        taxEffectMat
            = to_matrix(segment(taxEffects,
                                startIndex,
                                NFactors * NMicrobeNodes),
                        NFactors, NMicrobeNodes);
        taxNode_normFacts[h]
            = diag_post_multiply(
                diag_pre_multiply(
                    taxFactScales[,h],
                    taxEffectMat),
                microbeNodeScales);
        tax_normFacts[h] = taxNode_normFacts[h] * microbeAncestors;
    }
}
model {
    matrix[NSamples, NMicrobeTips] sampleTaxEffects[NSepFacts];
    matrix[NMicrobeTips, NSamples] logProps;
    matrix[NEnvs + 1, NMicrobeTips] tax_normEnvs;
    taxAveStD ~ exponential(1.0 / taxAveStDPriorExpect);
    envPropAveStD ~ exponential(1.0 / envPropAveStDPriorExpect);
    basePropsRaw ~ student_t(envPropsDF, 0, 10);
    to_vector(envPropMeansRaw) ~ student_t(envPropsDF, 0, 10);
    to_vector(envPropResiduals) ~ student_t(envPropsDF, 0, 1);
    taxIntercepts ~ cauchy(0, 2.5);
    taxEffects ~ student_t(taxEffectsDF, 0, 1);
    tax_normEnvs = envTaxMat * tax_normEnvEffects;
    for (h in 1:NSepFacts)
        sampleTaxEffects[h] = modelMat * tax_normFacts[h];
    for (k in 1:NSamples) {
        matrix[NSepFacts, NMicrobeTips] sampleTaxEffectsTransposed;
        matrix[NMicrobeTips, NEnvs + 1] allEnvSampleTaxEffects;
        matrix[NMicrobeTips, NEnvs + 1] allEnvSampleTaxLogProps;
        for (h in 1:NSepFacts)
            sampleTaxEffectsTransposed[h,] = sampleTaxEffects[h,k,];
        allEnvSampleTaxEffects = (tax_normEnvs
                                  + envFactWithUnk * sampleTaxEffectsTransposed)';
        for (i in 1:(NEnvs + 1)) {
            int startIndex = NMicrobeNodes * (NEnvs + NTotalFactors)
                             + (k - 1) * NMicrobeTips
                             + 1;
            allEnvSampleTaxLogProps[,i]
                = log_samp_props[k,i]
                  + log_softmax(allEnvSampleTaxEffects[,i]
                                + segment(taxEffects, startIndex, NMicrobeTips)
                                  .* taxResidualScales);
        }
        for (o in 1:NMicrobeTips)
            logProps[o,k] = log_sum_exp(allEnvSampleTaxLogProps[o,]) + multinomialNuisance[k];
    }
    to_array_1d(y) ~ poisson_log(to_array_1d(logProps));
}
generated quantities {
    simplex[NEnvs + 1] env_props[NEnvs];
    simplex[NEnvs + 1] samp_props[NSamples];
    for (j in 1:NEnvs)
        env_props[j] = softmax(envPropMeans[, j]);
    for (k in 1:NSamples)
        samp_props[k] = exp(log_samp_props[k]);
}
