data {
    int NEstSamples;
	int NObs;
	int NEstTaxNodes;
	int NFactors;
	int NFactorsDisp;
	int NIndiv;
	int NIndivPerTax[NEstTaxNodes];
	int iiidx[NEstTaxNodes]; //start of indices contained in indivIndices
	int indivIndices[sum(NIndivPerTax)]; // indices for each taxon that specify the individuals to estimate
	int NFactsPerTax[NEstTaxNodes];
	int fiidx[NEstTaxNodes];
	int factIndices[sum(NFactsPerTax)];
	int sampledcounts[NObs];
	int datacounts[NObs];
	int samplenames[NObs];
	int nodenames[NObs];
	matrix[NEstSamples, NFactors + NIndiv] modelMat;
	matrix[NEstSamples, NFactorsDisp] dispMat;
	matrix[NIndiv, NIndiv] distIndiv;
}
transformed data {
	int NIndivTax = sum(NIndivPerTax);
	int NFactTax = sum(NFactsPerTax);
	int NFIT = NFactTax + NIndivTax;
	matrix[NFactors, NEstTaxNodes] tax_normFactsZeros = rep_matrix(0, NFactors, NEstTaxNodes);
	matrix[NIndiv, NEstTaxNodes] tax_normIndivZeros = rep_matrix(0, NIndiv, NEstTaxNodes);
}
parameters {
	vector<lower=0, upper=pi()/2>[NFactors + 3] scalesUnif;
	vector[NFIT + NFactorsDisp + NEstTaxNodes + NObs] normals;
}
transformed parameters {
	vector<lower=0>[NFactors + 3] scales = 2.5 * tan(scalesUnif);
	cov_matrix[NIndiv] covIndivR = exp(-distIndiv/scales[NFactors+1]) * scales[NFactors+2]^2;
	matrix[NFactors, NEstTaxNodes] tax_normFactsRaw = tax_normFactsZeros;
	matrix[NFactors, NEstTaxNodes] tax_normFacts;
	matrix[NIndiv, NEstTaxNodes] tax_normIndiv = tax_normIndivZeros;
	for (o in 1:NEstTaxNodes) {
		int curridx[NIndivPerTax[o]] = segment(indivIndices, iiidx[o], NIndivPerTax[o]);
		tax_normFactsRaw[segment(factIndices, fiidx[o], NFactsPerTax[o]), o] = segment(normals, fiidx[o], NFactsPerTax[o]);
		tax_normIndiv[curridx, o] = cholesky_decompose(covIndivR[curridx, curridx]) * segment(normals, NFactTax + iiidx[o], NIndivPerTax[o]); //replace zeros in tax_normIndiv for taxon o only for individuals specified with the appropriate indivIndices
	}
	tax_normFacts = diag_pre_multiply(scales[1:NFactors], tax_normFactsRaw);
}
model {
	matrix[NEstSamples, NEstTaxNodes] sampleTaxEffects;
	matrix[NEstSamples, NEstTaxNodes] sampleTaxDispersion;
	vector[NObs] logit_ratios;
	normals ~ normal(0,1);
	sampleTaxEffects = modelMat * append_row(tax_normFacts, tax_normIndiv[1:NIndiv,]);
	sampleTaxDispersion = scales[NFactors + 3] * exp(rep_matrix(dispMat * segment(normals, NFIT + 1, NFactorsDisp), NEstTaxNodes) + rep_matrix(segment(normals, NFIT + NFactorsDisp + 1, NEstTaxNodes)', NEstSamples));
	for (n in 1:NObs)
		logit_ratios[n] = sampleTaxEffects[samplenames[n],nodenames[n]] + normals[NFIT + NFactorsDisp + NEstTaxNodes + n] * sampleTaxDispersion[samplenames[n],nodenames[n]];
	datacounts ~ binomial_logit(sampledcounts, logit_ratios);
}
generated quantities {
	vector[NFactorsDisp] factorDispersionMultipliers = exp(segment(normals, NFIT + 1, NFactorsDisp));
}
