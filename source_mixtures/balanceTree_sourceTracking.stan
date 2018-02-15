data {
    int NSamples;
	int NObs;
	int NEstTaxNodes;
	int NEnvs;
	int NFactors;
	int NFactsPerTax[NEstTaxNodes];
	int fiidx[NEstTaxNodes];
	int factIndices[sum(NFactsPerTax)];
	int envs[NSamples];
	matrix[NEstSamples, NFactors] modelMat;
	int sampledcounts[NObs];
	int datacounts[NObs];
	int samplenames[NObs];
	int nodenames[NObs];
}
transformed data {
	int NEnvSamps = NEnvs * NSamples;
	int NEnvTax = NEnvs * NEstTaxNodes;
	int NFactTax = sum(NFactsPerTax);
	matrix[NFactors, NEstTaxNodes] tax_normFactsZeros = rep_matrix(0, NFactors, NEstTaxNodes);
}
parameters {
	vector<lower=0, upper=pi()/2>[3 + NEnvs + NEstTaxNodes + NFactors] scalesUnif;
	matrix<upper=0>[NEnvs,NEnvs] env_prop_normal_raw;
	vector[NEnvs * (1 + NSamples + NEstTaxNodes) + NObs + NFactTax] normals;
}
transformed parameters {
	vector<lower=0>[3 + NEnvs + NEstTaxNodes + NFactors] scales = 2.5 * tan(scalesUnif);
	matrix[NEnvs,NEnvs] env_prop_normal;
	vector[NEnvs] log_samp_props[NSamples];
	matrix[NEnvs, NEstTaxNodes] log_tax_props_perEnv;
	matrix[NFactors, NEstTaxNodes] tax_normFactsRaw = tax_normFactsZeros;
	matrix[NFactors, NEstTaxNodes] tax_normFacts;
	for (i in 1:NEnvs) {
		vector[NEnvs] env_prop_normal_diffs = env_prop_normal_raw[,i];
		env_prop_normal_diffs[i] = -env_prop_normal_raw[i,i];
		env_prop_normal[,i] = scales[1] * env_prop_normal_diffs + scales[2] * normals[i];
	}
	for (k in 1:NSamples)
		log_samp_props[k] = log_softmax(env_prop_normal[envs[k],]' + segment(normals, k * NEnvs + 1, NEnvs) * scales[2 + envs[k]]);
	log_tax_props_perEnv = to_matrix(log_inv_logit(segment(normals, NEnvs + NEnvSamps + 1, NEnvTax) * scales[3 + NEnvs]), NEnvs, NEstTaxNodes);
	for (o in 1:NEstTaxNodes)
		tax_normFactsRaw[segment(factIndices, fiidx[o], NFactsPerTax[o]), o] = segment(normals, NEnvs + NEnvSamps + NEnvTax + fiidx[o], NFactsPerTax[o]);
	tax_normFacts = diag_pre_multiply(segment(scales, 3 + NEnvs + NEstTaxNodes + 1, NFactors), tax_normFactsRaw);
}
model {
	matrix[NEstSamples, NEstTaxNodes] sampleTaxEffects;
	vector[NObs] logit_ratios;
    to_vector(env_prop_normal_raw) ~ normal(0,1);
	normals ~ normal(0,1);
	sampleTaxEffects = modelMat * tax_normFacts;
	for (n in 1:NObs) {
		real log_compInt = log_sum_exp(log_tax_props_perEnv[,nodenames[n]] + log_samp_props[samplenames[n]]);
		logit_ratios[n] = log_compInt - log1m_exp(log_compInt) + sampleTaxEffects[samplenames[n],nodenames[n]] + normals[NEnvs + NEnvSamps + NEnvTax + n] * scales[3 + NEnvs + nodenames[n]];
	}
	datacounts ~ binomial_logit(sampledcounts, logit_ratios);
}
generated quantities {
	simplex[NEnvs] env_props[NEnvs];
	simplex[NEnvs] samp_props[NSamples];
	for (j in 1:NEnvs)
		env_props[j] = softmax(env_prop_normal[j,]');
	for (k in 1:NSamples)
		samp_props[k] = exp(log_samp_props[k]);
}
