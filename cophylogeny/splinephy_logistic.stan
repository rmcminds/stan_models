functions {
	vector build_b_spline(real[] t, real[] extKnots, int i, int order);
    vector build_b_spline(real[] t, real[] extKnots, int i, int order) {
        vector[size(t)] b_spline;
        vector[size(t)] w1 = rep_vector(0, size(t));
        vector[size(t)] w2 = rep_vector(0, size(t));
        if (order==1)
            for (j in 1:size(t))
                b_spline[j] = (extKnots[i] <= t[j]) && (t[j] < extKnots[i+1]);
        else {
        	if (extKnots[i] != extKnots[i+order-1])
            	w1 = (to_vector(t) - rep_vector(extKnots[i], size(t))) / (extKnots[i+order-1] - extKnots[i]);
        	if (extKnots[i+1] != extKnots[i+order])
            	w2 = 1 - (to_vector(t) - rep_vector(extKnots[i+1], size(t))) / (extKnots[i+order] - extKnots[i+1]);
        	// Calculating the B-spline recursively as linear interpolation of two lower-order splines
        	b_spline = w1 .* build_b_spline(t, extKnots, i, order-1) + w2 .* build_b_spline(t, extKnots, i+1, order-1);
    	}
    	return b_spline;
	} // http://mc-stan.org/users/documentation/case-studies/splines_in_stan.html
	matrix build_B(int splineDegree, vector knots, real[] cophLT, int NBasis) {
		int NKnots = rows(knots);
		int NLT = size(cophLT);
		vector[splineDegree + NKnots] extKnots_temp = append_row(rep_vector(knots[1], splineDegree), knots);
		vector[2 * splineDegree + NKnots] extKnots = append_row(extKnots_temp, rep_vector(knots[NKnots], splineDegree));
		matrix[NBasis, NLT] B;
		for (i in 1:NBasis)
			B[i,] = to_row_vector(build_b_spline(cophLT, to_array_1d(extKnots), i, splineDegree + 1));
		B[NKnots + splineDegree - 1, NLT] = 1;
		return B;
	}
	matrix phyCorr(real intercept, matrix cophenetic, row_vector a, matrix B) {
		int NCol = cols(cophenetic);
		matrix[NCol, NCol] rawCov;
		matrix[NCol, NCol] cov;
		matrix[NCol, NCol] corr;
		row_vector[cols(B)] aB = a * B;
		matrix[NCol, NCol] raw = rep_matrix(0, NCol, NCol);
		{
			int ind = 1;
            for (i in 1:NCol) {
                raw[i, 1:i] = segment(aB, ind, i);
				raw[1:(i-1), i] = segment(aB, ind, i-1)';
                ind = ind + i;
            }
        }
		rawCov = intercept * cophenetic + raw;
		cov = rawCov - diag_matrix(rep_vector(1.1 * min(eigenvalues_sym(rawCov)), NCol));
		corr = quad_form_diag(cov, inv_sqrt(diagonal(cov)));
		return corr;
	}
}
data {
	int NSamples;
	int NHosts;
	int NSym;
	int y[NSym, NSamples];
	matrix[NSamples, NHosts] modelMat;
	matrix[NHosts, NHosts] hostCophenetic;
	matrix[NSym, NSym] symCophenetic;
	int NHostCophLT;
	real hostCophLT[NHostCophLT];
	int NSymCophLT;
	real symCophLT[NSymCophLT];
	int NHostKnots;
	vector[NHostKnots] hostKnots;
	int NSymKnots;
	vector[NSymKnots] symKnots;
	int splineDegree;
	vector[NSamples] sampleCount;
}
transformed data {
	int NBasisHost = NHostKnots + splineDegree - 1;
	int NBasisSym = NSymKnots + splineDegree - 1;
	matrix[NBasisHost, NHostCophLT] BHost = build_B(splineDegree, hostKnots, hostCophLT, NBasisHost);
	matrix[NBasisSym, NSymCophLT] BSym = build_B(splineDegree, symKnots, symCophLT, NBasisSym);
	int NHostSym = NHosts * NSym;
	int NHostSymBH = NHostSym + NBasisHost;
	int NHostSymBHS = NHostSymBH + NBasisSym;
	vector[NSamples] log_sampleCount = log(sampleCount);
}
parameters {
	vector<lower=0, upper=pi()/2>[NHosts + NSym + 2] scalesUnif;
	vector<lower=0>[2] splineScales;
	vector[NHostSymBHS + 4] effects;
}
transformed parameters {
	vector<lower=0>[NHosts + NSym + 2] scales = 2.5 * tan(scalesUnif);
	row_vector[NBasisHost] aHost = to_row_vector(splineScales[1] * segment(effects, NHostSym + 1, NBasisHost));
	row_vector[NBasisSym] asym = to_row_vector(splineScales[2] * segment(effects, NHostSymBH + 1, NBasisSym));
	matrix[NHosts,NHosts] hostCorr = phyCorr(effects[NHostSymBHS + 1], hostCophenetic, aHost, BHost);
	matrix[NSym,NSym] symCorr = phyCorr(effects[NHostSymBHS + 2], symCophenetic, asym, BSym);
	cholesky_factor_cov[NHosts] cholesky_hostVCV = diag_pre_multiply(segment(scales, 1, NHosts), cholesky_decompose(hostCorr));
	cholesky_factor_cov[NSym] cholesky_symVCV = diag_pre_multiply(segment(scales, NHosts + 1, NSym), cholesky_decompose(symCorr));
	matrix[NHosts, NSym] hostBySymEfects = cholesky_hostVCV * to_matrix(segment(effects, 1, NHostSym), NHosts, NSym) * cholesky_symVCV'; //kronecker(A,B)*vec(X) = B*X*t(A)
	real intercept = scales[NHosts + NSym + 1] * effects[NHostSymBHS + 3];
	real sampleCountEffect = scales[NHosts + NSym + 2] * effects[NHostSymBHS + 4];
}
model {
	matrix[NSamples, NSym] sampleEstimates;
	effects ~ normal(0,1);
	splineScales ~ normal(0,1);
    sampleEstimates = intercept + modelMat * hostBySymEfects + rep_matrix(sampleCountEffect * log_sampleCount, NSym);
	to_array_1d(y) ~ bernoulli_logit(to_vector(sampleEstimates));
}
