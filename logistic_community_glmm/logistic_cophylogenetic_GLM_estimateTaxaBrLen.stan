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
	real<lower=0> taxAveStDPriorExpect;
	vector[NTimeBins] timeBinSizes;
	matrix[NEstTaxNodes, NEstTips] ancestors;
	matrix[NEffects, NFactors] factLevelMat;
	matrix[NEstSamples, NEffects] modelMat;
	matrix[NHostNodes, NTimeBins] edgetobin;
	matrix[NEstSamples, NHostNodes] hostAncestors;
}
parameters {
    row_vector<lower=0>[NEstTaxNodes] taxNodeScalesRaw;
	real<lower=0> aveStD;
	simplex[2 * NFactors + 1] stDProps;
	simplex[NTimeBins] timeBinPropsADiv;
	simplex[NTimeBins] timeBinPropsSpec;
	real globalIntercept;
	vector[NEffects - 1 + NHostNodes] rawAlphaDivEffects;
	matrix[NEffects + NHostNodes, NEstTaxNodes] rawTaxNodeEffects;
}
transformed parameters {
    row_vector[NEstTaxNodes] taxNodeScales;
	vector<lower=0>[2 * NFactors + 1] scales;
	vector<lower=0>[NHostNodes] newEdgesADiv;
	vector<lower=0>[NHostNodes] hostNodeScalesADiv;
	vector[NEffects + NHostNodes] scaledAlphaDivEffects;
	vector<lower=0>[NHostNodes] newEdgesSpec;
	vector<lower=0>[NHostNodes] hostNodeScalesSpec;
	matrix[NEffects + NHostNodes, NEstTaxNodes] scaledTaxNodeEffects;

    taxNodeScales
        = sqrt(square(taxNodeScalesRaw)
               / mean(square(taxNodeScalesRaw) * ancestors));
	scales
        = sqrt((2 * NFactors + 1)
               * stDProps)
          * aveStD;
	newEdgesADiv
        = edgetobin * timeBinPropsADiv;
	hostNodeScalesADiv
        = sqrt(newEdgesADiv) * scales[NFactors];
	scaledAlphaDivEffects
        = append_row(globalIntercept,
            append_row(factLevelMat[2:NEffects, 2:NFactors]
                       * segment(scales, 1, NFactors - 1),
                       hostNodeScalesADiv)
            .* rawAlphaDivEffects);
	newEdgesSpec
        = edgetobin * timeBinPropsSpec;
	hostNodeScalesSpec
        = sqrt(newEdgesSpec) * scales[2 * NFactors + 1];
	scaledTaxNodeEffects
        = diag_post_multiply(
            diag_pre_multiply(
                append_row(factLevelMat * segment(scales, NFactors + 1, NFactors),
                    hostNodeScalesSpec),
                rawTaxNodeEffects),
            taxNodeScales);
}
model {
	matrix[NEstSamples, NEstTips] sampleTipEffects;
	vector[NObs] logit_ratios;

    taxNodeScalesRaw ~ exponential(1);
	aveStD ~ exponential(1.0 / taxAveStDPriorExpect);
	stDProps ~ dirichlet(rep_vector(1, 2 * NFactors + 1));
	timeBinPropsADiv ~ dirichlet(NTimeBins * timeBinSizes);
	timeBinPropsSpec ~ dirichlet(NTimeBins * timeBinSizes);
	globalIntercept ~ normal(0, 50);
	rawAlphaDivEffects ~ normal(0,1);
	to_vector(rawTaxNodeEffects) ~ normal(0,1);
	sampleTipEffects = append_col(modelMat,hostAncestors)
                       * (rep_matrix(scaledAlphaDivEffects, NEstTips)
                          + scaledTaxNodeEffects
                            * ancestors);
	for (n in 1:NObs)
		logit_ratios[n] = sampleTipEffects[samplenames[n], tipnames[n]];
	present ~ bernoulli_logit(logit_ratios);
}
generated quantities {
	vector[NTimeBins] relativeEvolRatesADiv = timeBinPropsADiv ./ timeBinSizes;
	vector[NTimeBins] relativeEvolRatesSpec = timeBinPropsSpec ./ timeBinSizes;
}
