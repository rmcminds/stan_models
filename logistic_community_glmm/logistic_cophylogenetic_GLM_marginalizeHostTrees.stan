data {
    int NEstSamples;
	int NObs;
	int NEstTaxNodes;
	int NEstTips;
	int NFactors;
	int NEffects;
	int NHostNodes;
	int NTimeBins;
    int present[NObs];
	int samplenames[NObs];
	int tipnames[NObs];
	int NTrees;
	real<lower=0> taxAveStDPriorExpect;
	vector[NTimeBins] timeBinSizes[NTrees];
	matrix[NEstTaxNodes, NEstTips] ancestors;
	matrix[NEffects, NFactors] factLevelMat;
	matrix[NEstSamples, NEffects] modelMat;
	matrix[NHostNodes, NTimeBins] edgetobin[NTrees];
	matrix[NEstSamples, NHostNodes] hostAncestors[NTrees];
	vector[NEstTaxNodes] taxNodeScales;
}
parameters {
	real<lower=0> aveStD[NTrees];
	simplex[2 * NFactors + 1] stDPropList[NTrees];
	simplex[NTimeBins] timeBinPropsADiv[NTrees];
	simplex[NTimeBins] timeBinPropsSpec[NTrees];
	simplex[NTrees] treeProbs[NTrees];
	real globalIntercept[NTrees];
	vector[NEffects - 1 + NHostNodes] rawAlphaDivEffects[NTrees];
	matrix[NEffects + NHostNodes, NEstTaxNodes] rawTaxNodeEffects[NTrees];
}
transformed parameters {
	vector<lower=0>[2 * NFactors + 1] scales[NTrees];
	vector<lower=0>[NHostNodes] newEdgesADiv[NTrees];
	vector<lower=0>[NHostNodes] hostNodeScalesADiv[NTrees];
	vector[NEffects + NHostNodes] scaledAlphaDivEffects[NTrees];
	vector<lower=0>[NHostNodes] newEdgesSpec[NTrees];
	vector<lower=0>[NHostNodes] hostNodeScalesSpec[NTrees];
	matrix[NEffects + NHostNodes, NEstTaxNodes] scaledTaxNodeEffects[NTrees];
	matrix[NEffects + NHostNodes, NEstTips] scaledTipEffects[NTrees];
	for (t in 1:NTrees) {
		scales[t] = sqrt((2 * NFactors + 1) * stDPropList[t]) * aveStD[t];
        newEdgesADiv[t] = edgetobin[t] * timeBinPropsADiv[t];
        hostNodeScalesADiv[t] = sqrt(newEdgesADiv[t]) * scales[t][NFactors];
        scaledAlphaDivEffects[t] = append_row(globalIntercept[t], append_row(factLevelMat[2:NEffects, 2:NFactors] * segment(scales[t], 1, NFactors - 1), hostNodeScalesADiv[t]) .* rawAlphaDivEffects[t]);
        newEdgesSpec[t] = edgetobin[t] * timeBinPropsSpec[t];
        hostNodeScalesSpec[t] = sqrt(newEdgesSpec[t]) * scales[t][2 * NFactors + 1];
        scaledTaxNodeEffects[t] = diag_post_multiply(diag_pre_multiply(append_row(factLevelMat * segment(scales[t], NFactors + 1, NFactors), hostNodeScalesSpec[t]), rawTaxNodeEffects[t]),taxNodeScales);
        scaledTipEffects[t] = scaledTaxNodeEffects[t] * ancestors;
	}
}
model {
	matrix[NEstSamples, NEstTips] logit_ratios[NTrees];
	real logTreeProbs;
	treeProbs ~ dirichlet(rep_vector(1,NTrees));
	logTreeProbs = log(treeProbs)
	for (t in 1:NTrees) {
		aveStD[t] ~ exponential(1.0/taxAveStDPriorExpect);
		globalIntercept[t] ~ normal(0, 50);
		stDPropList[t] ~ dirichlet(rep_vector(1, 2 * NFactors + 1));
    	timeBinPropsADiv[t] ~ dirichlet(NTimeBins * timeBinSizes[t]);
    	timeBinPropsSpec[t] ~ dirichlet(NTimeBins * timeBinSizes[t]);
		rawAlphaDivEffects[t] ~ normal(0,1);
		to_vector(rawTaxNodeEffects[t]) ~ normal(0,1);
	}
	for (t in 1:NTrees)
		logit_ratios[t] = append_col(modelMat,hostAncestors[t]) * (rep_matrix(scaledAlphaDivEffects[t], NEstTips) + scaledTipEffects[t]);
	for (n in 1:NObs) {
		vector[NTrees] lps;
		for(t in 1:NTrees) {
			lps[t] = logTreeProbs[t] + bernoulli_logit_lpmf(present[n] | logit_ratios[t][samplenames[n], tipnames[n]]);
		}
		target += log_sum_exp(lps);
	}
}
generated quantities {
	int treeSample = categorical_rng(treeProbs);
	vector[NTimeBins] relativeEvolRatesADiv = timeBinPropsADiv[treeSample] ./ timeBinSizes[treeSample];
	vector[NTimeBins] relativeEvolRatesSpec = timeBinPropsSpec[treeSample] ./ timeBinSizes[treeSample];
	vector[2 * NFactors + 1] stDProps = stDPropList[treeSample];
}
