data {
    int NEstSamples;
	int NObs;
	int NEstTaxNodes;
	int NFactors;
	int NFactorsDisp;
	int NFactsPerTax[NEstTaxNodes];
	int fiidx[NEstTaxNodes];
	int factIndices[sum(NFactsPerTax)];
	int sampledcounts[NObs];
	int datacounts[NObs];
	int samplenames[NObs];
	int nodenames[NObs];
	matrix[NEstSamples, NFactors] modelMat;
	matrix[NEstSamples, NFactorsDisp] dispMat;
}
transformed data {
	int NFactTax = sum(NFactsPerTax);
	matrix[NFactors, NEstTaxNodes] tax_normFactsZeros = rep_matrix(0, NFactors, NEstTaxNodes);
}
parameters {
	vector<lower=0, upper=pi()/2>[NFactors + 1] scalesUnif;
	vector[NFactTax + NFactorsDisp + NEstTaxNodes + NObs] normals;
}
transformed parameters {
	vector<lower=0>[NFactors + 1] scales = 2.5 * tan(scalesUnif);
	matrix[NFactors, NEstTaxNodes] tax_normFactsRaw = tax_normFactsZeros;
	matrix[NFactors, NEstTaxNodes] tax_normFacts;
	for (o in 1:NEstTaxNodes) {
		tax_normFactsRaw[segment(factIndices, fiidx[o], NFactsPerTax[o]), o] = segment(normals, fiidx[o], NFactsPerTax[o]);
	}
	tax_normFacts = diag_pre_multiply(scales[1:NFactors], tax_normFactsRaw);
}
model {
	matrix[NEstSamples, NEstTaxNodes] sampleTaxEffects;
	matrix[NEstSamples, NEstTaxNodes] sampleTaxDispersion;
	vector[NObs] logit_ratios;
	normals ~ normal(0,1);
	sampleTaxEffects = modelMat * tax_normFacts;
	sampleTaxDispersion = scales[NFactors + 1] * exp(rep_matrix(dispMat * segment(normals, NFactTax + 1, NFactorsDisp), NEstTaxNodes) + rep_matrix(segment(normals, NFactTax + NFactorsDisp + 1, NEstTaxNodes)', NEstSamples));
	for (n in 1:NObs)
		logit_ratios[n] = sampleTaxEffects[samplenames[n],nodenames[n]] + normals[NFactTax + NFactorsDisp + NEstTaxNodes + n] * sampleTaxDispersion[samplenames[n],nodenames[n]];
	datacounts ~ binomial_logit(sampledcounts, logit_ratios);
}
generated quantities {
	vector[NFactorsDisp] factorDispersionMultipliers = exp(segment(normals, NFactTax + 1, NFactorsDisp));
}
