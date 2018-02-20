data {
    int NEstSamples;
	int NObs;
	int NEstTaxNodes;
	int NEnvs;
	int NFactors;
	int envs[NEstSamples];
	matrix[NEstSamples, NFactors] modelMat;
	int sampledcounts[NObs];
	int datacounts[NObs];
	int samplenames[NObs];
	int nodenames[NObs];
}
transformed data {
	int NEnvSamps = NEnvs * NEstSamples;
	int NEnvTax = NEnvs * NEstTaxNodes;
	int NFactTax = NFactors * NEstTaxNodes;
}
parameters {
	vector<lower=0>[3 + NEnvs + NFactors + NEstTaxNodes] scales;
	matrix<upper=0>[NEnvs,NEnvs] env_prop_normal_raw;
	vector[NEnvs + NEnvSamps + NEnvTax + NFactTax + NObs] normals;
}
transformed parameters {
	matrix[NEnvs,NEnvs] env_prop_normal;
	vector[NEnvs] log_samp_props[NEstSamples];
	matrix[NEnvs, NEstTaxNodes] tax_normEnvs;
	matrix[NFactors, NEstTaxNodes] tax_normFacts;
	for (i in 1:NEnvs) {
		vector[NEnvs] env_prop_normal_diffs = env_prop_normal_raw[,i];
		env_prop_normal_diffs[i] = -env_prop_normal_raw[i,i];
		env_prop_normal[,i] = scales[1] * env_prop_normal_diffs + scales[2] * normals[i];
	}
	for (k in 1:NEstSamples)
		log_samp_props[k] = log_softmax(env_prop_normal[envs[k],]' + segment(normals, k * NEnvs + 1, NEnvs) * scales[2 + envs[k]]);
	tax_normEnvs = scales[2 + NEnvs + 1] * to_matrix(segment(normals, NEnvs + NEnvSamps + 1, NEnvTax), NEnvs, NEstTaxNodes);
	tax_normFacts = diag_pre_multiply(segment(scales, 3 + NEnvs + 1, NFactors), to_matrix(segment(normals, NEnvs + NEnvSamps + NEnvTax + 1, NFactTax), NFactors, NEstTaxNodes));
}
model {
	matrix[NEstSamples, NEstTaxNodes] sampleTaxEffects;
	vector[NObs] logit_ratios;
	scales ~ student_t(5,0,2.5);
    to_vector(env_prop_normal_raw) ~ normal(0,1);
	normals ~ normal(0,1);
	sampleTaxEffects = modelMat * tax_normFacts;
	for (n in 1:NObs) {
		vector[NEnvs] sampleTaxEnvErr = log_inv_logit(tax_normEnvs[,nodenames[n]] + sampleTaxEffects[samplenames[n], nodenames[n]] + normals[NEnvs + NEnvSamps + NEnvTax + NFactTax + n] * scales[3 + NEnvs + NFactors + nodenames[n]]);
		real log_ratios = log_sum_exp(sampleTaxEnvErr + log_samp_props[samplenames[n]]);
		logit_ratios[n] = log_ratios - log1m_exp(log_ratios);
	}
	datacounts ~ binomial_logit(sampledcounts, logit_ratios);
}
generated quantities {
	simplex[NEnvs] env_props[NEnvs];
	simplex[NEnvs] samp_props[NEstSamples];
	for (j in 1:NEnvs)
		env_props[j] = softmax(env_prop_normal[j,]');
	for (k in 1:NEstSamples)
		samp_props[k] = exp(log_samp_props[k]);
}
