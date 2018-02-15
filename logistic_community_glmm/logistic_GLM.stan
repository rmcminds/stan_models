data {
    int NEstSamples;
	int NObs;
	int NEstTaxNodes;
	int NFactors;
	int NFactsPerTax[NEstTaxNodes];
	int fiidx[NEstTaxNodes];
	int factIndices[sum(NFactsPerTax)];
	int NEffects;
    int present[NObs];
	int samplenames[NObs];
	int nodenames[NObs];
	matrix[NEstSamples, NEffects] modelMat;
	matrix[NEffects, NFactors] factLevelMat;
}
transformed data {
	int NFactTax = sum(NFactsPerTax);
	matrix[NEffects, NEstTaxNodes] tax_normFactsZeros = rep_matrix(0, NEffects, NEstTaxNodes);
}
parameters {
	vector<lower=0, upper=pi()/2>[NFactors * 2] scalesUnif;
	vector[NFactTax + NEffects] normals;
}
transformed parameters {
	vector<lower=0>[NFactors * 2] scales = 2.5 * tan(scalesUnif);
	vector[NEffects] global = factLevelMat * segment(scales, NFactors + 1, NFactors) .* segment(normals, NFactTax + 1, NEffects);
	matrix[NEffects, NEstTaxNodes] tax_normFactsRaw = tax_normFactsZeros;
	matrix[NEffects, NEstTaxNodes] tax_normFacts;
	for (o in 1:NEstTaxNodes) {
		int factidx[NFactsPerTax[o]] = segment(factIndices, fiidx[o], NFactsPerTax[o]);
		tax_normFactsRaw[factidx, o] = segment(normals, fiidx[o], NFactsPerTax[o]);
	}
	tax_normFacts = diag_pre_multiply(factLevelMat * scales[1:NFactors], tax_normFactsRaw);
}
model {
	vector[NEstSamples] sampleEffects;
	matrix[NEstSamples, NEstTaxNodes] sampleTaxEffects;
	vector[NObs] logit_ratios;
	normals ~ normal(0,1);
	sampleEffects = modelMat * global;
	sampleTaxEffects = modelMat * tax_normFacts;
	for (n in 1:NObs)
        logit_ratios[n] = sampleEffects[samplenames[n]] + sampleTaxEffects[samplenames[n],nodenames[n]];
	present ~ bernoulli_logit(logit_ratios);
}
