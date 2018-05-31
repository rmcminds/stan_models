data {
    int NEstSamples;
	int NEstTaxNodes;
	int NEstTips;
	int NEnvs;
	int NFactors;
	int envs[NEstSamples];
	int NSepFacts;
	matrix[NEnvs, NSepFacts] envFact;
	matrix[NEstSamples, NFactors] modelMat;
	int y[NEstTips, NEstSamples];
	real taxAveStDPriorExpect;
	real envPropAveStDPriorExpect;
	matrix[NEstTaxNodes, NEstTips] ancestors;
	vector[NEstTaxNodes] taxNodeScales;
}
transformed data {
	real<lower=2.0> taxEffectsDF = 5.0;
	real<lower=2.0> envPropsDF = 5.0;
	int NEnvSamps = NEnvs * NEstSamples;
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
	matrix<upper=0>[NEnvs, NEnvs] envPropMeansRaw;
	matrix[NEnvs, NEstSamples] envPropResiduals;
	row_vector[NEstTaxNodes] taxIntercepts;
	simplex[NSepFacts * NFactors + 2] taxVarianceProps;
	simplex[NEstTips] taxResidualVarianceProps;
	vector[(NEnvs + NSepFacts * NFactors + NEstSamples) * NEstTaxNodes] taxEffects;
	vector[NEstSamples] multinomialNuisance;
}
transformed parameters {
	vector<lower=0>[NEnvs] envPropResidualScales = sqrt((envPropsDF-2) / envPropsDF * NEnvs * envPropResidualProps) * envPropAveStD;
	matrix<lower=0>[NFactors, NSepFacts] taxFactScales = to_matrix(sqrt((taxEffectsDF-2) / taxEffectsDF * (NSepFacts * NFactors + 2) * segment(taxVarianceProps, 1, NSepFacts * NFactors)) * taxAveStD, NFactors, NSepFacts);
	real<lower=0> taxEnvScale = sqrt((taxEffectsDF-2) / taxEffectsDF * (NSepFacts * NFactors + 2) * taxVarianceProps[NSepFacts * NFactors + 1]) * taxAveStD;
	vector<lower=0>[NEstTips] taxResidualScales = sqrt((taxEffectsDF-2) / taxEffectsDF * (NSepFacts * NFactors + 2) * taxVarianceProps[NSepFacts * NFactors + 2] * NEstTips * taxResidualVarianceProps) * taxAveStD;
	matrix[NEnvs + 1, NEnvs] envPropMeans;
	vector[NEnvs + 1] log_samp_props[NEstSamples];
	matrix[NEnvs + 1, NEstTaxNodes] taxNode_normEnvEffects;
	matrix[NEnvs + 1, NEstTips] tax_normEnvEffects;
	matrix[NFactors, NEstTaxNodes] taxNode_normFacts[NSepFacts];
	matrix[NFactors, NEstTips] tax_normFacts[NSepFacts];
	vector[NEnvs + 1] baseProps = append_row(basePropsRaw, 0);
    for (i in 1:NEnvs) {
        envPropMeans[1:(i-1), i] = baseProps[1:(i-1)] + envPropMeansRaw[1:(i-1), i];
        envPropMeans[i, i] = baseProps[i];
        envPropMeans[(i+1):(NEnvs+1), i] = baseProps[(i+1):(NEnvs+1)] + envPropMeansRaw[i:NEnvs, i];
    }
    for (k in 1:NEstSamples) {
        vector[NEnvs + 1] samp_prop_normal;
        samp_prop_normal[1:(envs[k]-1)] = envPropResidualScales[envs[k]] * envPropResiduals[1:(envs[k]-1), k] + envPropMeans[1:(envs[k]-1), envs[k]];
        samp_prop_normal[envs[k]] = envPropMeans[envs[k], envs[k]];
        samp_prop_normal[(envs[k]+1):(NEnvs+1)] = envPropResidualScales[envs[k]] * envPropResiduals[envs[k]:NEnvs, k] + envPropMeans[(envs[k]+1):(NEnvs+1), envs[k]];
        log_samp_props[k] = log_softmax(samp_prop_normal);
    }
	taxNode_normEnvEffects = diag_post_multiply(append_row(taxIntercepts, taxEnvScale * to_matrix(segment(taxEffects, 1, NEnvs * NEstTaxNodes), NEnvs, NEstTaxNodes)), taxNodeScales);
	tax_normEnvEffects = taxNode_normEnvEffects * ancestors;
	for (h in 1:NSepFacts) {
		taxNode_normFacts[h] = diag_post_multiply(diag_pre_multiply(taxFactScales[,h], to_matrix(segment(taxEffects, NEnvs * NEstTaxNodes + NFactors * NEstTaxNodes * (h-1) + 1, NFactors * NEstTaxNodes), NFactors, NEstTaxNodes)), taxNodeScales);
		tax_normFacts[h] = taxNode_normFacts[h] * ancestors;
	}
}
model {
	matrix[NEstSamples, NEstTips] sampleTaxEffects[NSepFacts];
	matrix[NEstTips, NEstSamples] logProps;
	matrix[NEnvs + 1, NEstTips] tax_normEnvs;
	taxAveStD ~ exponential(1.0 / taxAveStDPriorExpect);
	envPropAveStD ~ exponential(1.0 / envPropAveStDPriorExpect);
	basePropsRaw ~ student_t(envPropsDF, 0, 10);
	to_vector(envPropMeansRaw) ~ student_t(envPropsDF, 0, NEnvs * 10);
	to_vector(envPropResiduals) ~ student_t(envPropsDF, 0, 1);
	taxIntercepts ~ cauchy(0, 2.5);
	taxEffects ~ student_t(taxEffectsDF, 0, 1);
	tax_normEnvs = envTaxMat * tax_normEnvEffects;
	for (h in 1:NSepFacts)
		sampleTaxEffects[h] = modelMat * tax_normFacts[h];
	for (k in 1:NEstSamples) {
		matrix[NSepFacts, NEstTips] sampleTaxEffectsTransposed;
		matrix[NEstTips, NEnvs + 1] allEnvSampleTaxEffects;
		matrix[NEstTips, NEnvs + 1] allEnvSampleTaxLogProps;
		for (h in 1:NSepFacts)
			sampleTaxEffectsTransposed[h,] = sampleTaxEffects[h,k,];
		allEnvSampleTaxEffects = (tax_normEnvs + envFactWithUnk * sampleTaxEffectsTransposed)';
		for (i in 1:(NEnvs + 1))
			allEnvSampleTaxLogProps[,i] = log_samp_props[k,i] + log_softmax(allEnvSampleTaxEffects[,i] + segment(taxEffects, NEnvs * NEstTaxNodes + NSepFacts * NFactors * NEstTaxNodes + (k - 1) * NEstTips + 1, NEstTips) .* taxResidualScales);
		for (o in 1:NEstTips)
			logProps[o,k] = log_sum_exp(allEnvSampleTaxLogProps[o,]) + multinomialNuisance[k];
	}
	to_array_1d(y) ~ poisson_log(to_array_1d(logProps));
}
generated quantities {
	simplex[NEnvs + 1] env_props[NEnvs];
	simplex[NEnvs + 1] samp_props[NEstSamples];
	for (j in 1:NEnvs)
		env_props[j] = softmax(envPropMeans[,j]);
	for (k in 1:NEstSamples)
		samp_props[k] = exp(log_samp_props[k]);
}
