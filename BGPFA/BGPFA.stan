//strong inspiration from https://www.cs.helsinki.fi/u/sakaya/tutorial/code/cca.R
#include functions.stan
data {
    int<lower=0> N;                                     // Number of samples
    int<lower=0> D;                                     // Number of compositional count datasets
    int<lower=0> R;                                     // Number of continuous datasets
    int<lower=0> C;                                     // Number of categorical datasets
    int<lower=0> nVarGroups;                            // Number of groups of samples that have a shared single observation
    int<lower=0> M[D+R+C];                              // Number of observed variables in each dataset
    int<lower=0> Mplus[D+R+C];                          // Number of parameters for each dataset
    int<lower=0> K_linear;                              // Number of residual latent linear dimensions
    int<lower=0> K_gp;                                  // Number of residual latent GP dimensions per kernel
    int<lower=0> KG;                                    // Number of unique GP kernels
    int<lower=0> I[D, N + nVarGroups];                  // Sample indices for count datasets
    int<lower=0> O;                                     // Total number of observed continuous variables
    int<lower=0> IR[O, N + nVarGroups];                 // Sample indices for continuous variables
    int<lower=0> C_vars;                                // Total number of observed categorical variables
    int<lower=0> IC[C_vars, N + nVarGroups];            // Sample indices for categorical datasets
    int<lower=0> Mc[C_vars];                            // Number of levels for each categorical variable
    int<lower=0> F;                                     // Total number of compositional count observations
    int<lower=0> X[F];                                  // compositional count data
    int<lower=0> G;                                     // Total number of categorical observations
    vector<lower=0>[G] Y;                               // categorical data
    int<lower=0> H;                                     // Total number of continuous observations
    vector[H] P;                                        // continuous data
    real<lower=0> global_scale_prior;
    real<lower=0> ortho_scale_prior;                    // prior degree of orthogonalization regularization
    vector<lower=0>[sum(Mplus)] prior_scales;
    vector<lower=0>[sum(M)] prior_intercept_scales;
    vector[sum(M)] prior_intercept_centers;
    vector[sum(M[1:D])] binary_count_intercept_centers;
    int<lower=0> sizeMM;                                // size of model matrix for higher level variables
    vector[sizeMM] mm;                                  // model matrix for higher level variables
    int<lower=0> nMP;                                   // number of P measured as -inf
    int<lower=0> indMP[nMP];                            // indices for missing Ps
    real P_max[nMP];                                    // lowest measured value for the variable corresponding to a missing P
    real rate_gamma_fact;                               // scale of inverse gamma prior on nu for W
    real shape_gamma_fact;                              // shape of inverse gamma prior on nu for W
    vector[choose(M[size(M)],2)] dist_sites;            // distance between sites
    real rho_sites_prior;                               // prior expectation for length scale of gaussian process on sites
    matrix[N,nVarGroups] samp2group;                    // Matrix to calculate means of grouped samples
    int site_smoothness;                                // Matern smoothness parameter for sites
    real nu_residuals;                                  // residual robustness
    vector[D] inv_log_max_contam;                       // prior expectation of contamination rate
    real<lower=0> gnorm_shape;                          // strength of prior pulling contamination toward zero
}
transformed data {
    int K = K_linear + KG * K_gp;                     // Total number of latent axes
    int DRC = D+R+C;                                  // Total number of datasets
    int V = 0;                                        // Total number of observed compositional count variables
    int Vplus = 0;                                    // Total number of observed compositional count parameters
    int VOB = sum(M);                                 // Total number of observed variables
    int VOBplus = sum(Mplus);                         // Total number of parameters per dim
    int sumM[DRC+1];                                  // Cumulative number of base level variables
    int sumMplus[DRC+1];                              // Cumulative number of base level and higher variables
    int MMplus[DRC];                                  // size of model matrix for higher level variables
    int sumMMplus[DRC+1];                             // Cumulative size of model matrices for higher level variables
    int sumMc[C_vars+1];                              // Cumulative number of categorical levels
    int nmulti = 0;                                   // Total number of multinomial samples
    int sumID[D];                                     // Number of samples in each count dataset
    int IDInds[D,N+nVarGroups];                       // Indices mapping each sample to each count dataset
    int sumMhigherD[D];
    int F_higher;
    int sumIR[R,N+nVarGroups];                        // Number of samples in each continuous dataset
    int IRInds[R,N+nVarGroups,max(M[(D+1):(D+R)])];   // Indices mapping each sample to each continuous dataset
    int segInds1[R,N+nVarGroups,max(M[(D+1):(D+R)])]; // Various indices dealing with independent sets of continuous variables
    int segInds2[R,N+nVarGroups,max(M[(D+1):(D+R)])];
    int nIsolate[R,N+nVarGroups]
        = rep_array(0,R,N+nVarGroups);
    int nMat[R,N+nVarGroups]
        = rep_array(0,R,N+nVarGroups);
    int matchInds[R,N+nVarGroups]
        = rep_array(0,R,N+nVarGroups);
    int matchIndsInv[R,N+nVarGroups,N+nVarGroups]
        = rep_array(0,R,N+nVarGroups,N+nVarGroups);
    int nMatches[R,N+nVarGroups]
        = rep_array(0,R,N+nVarGroups);
    int nUniqueR[R] = rep_array(1,R);
    int uniqueRInds[R,N+nVarGroups,max(M[(D+1):(D+R)])]
        = rep_array(0,R,N+nVarGroups,max(M[(D+1):(D+R)]));
    int sumIRUnique[R,N+nVarGroups];                         //
    real rho_Z_shape = 1.0 / K_linear;                       // More 'independent' latent variables means more need for sparsity in feature selection
    real rho_Z_scale = 6.0 / K_gp;                           // More 'dependent' latent variables per group means, in order to encourage orthogonality, length scale must be smaller
    matrix[site_smoothness-1, 2 * site_smoothness] chooseRJ; // precomputed part of Matern covariance
    matrix[site_smoothness, 2 * site_smoothness] ffKJ;       // precomputed part of Matern covariance
    // create indices for all datasets
    sumM[1] = 0;
    sumMplus[1] = 0;
    sumMMplus[1] = 0;
    for(drc in 1:DRC) {
        sumM[drc+1] = sum(M[1:drc]);
        sumMplus[drc+1] = sum(Mplus[1:drc]);
        MMplus[drc] = M[drc] * (Mplus[drc] - M[drc]);
        sumMMplus[drc+1] = sum(MMplus[1:drc]);
    }
    //
    // create indices for count datasets
    V = sumM[D+1];
    Vplus = sumMplus[D+1];
    for(d in 1:D) {
        int nplace = 1;
        sumID[d] = sum(I[d,]);
        sumMhigherD[d] = sumID[d] * (Mplus[d] - M[d]);
        for(n in 1:(N+nVarGroups)) {
            if(I[d,n]) {
                IDInds[d,nplace] = n;
                nplace += 1;
            }
        }
    }
    nmulti = sum(sumID);
    F_higher = sum(sumMhigherD);
    //
    // create indices for sets of independent continuous variables. would love simplification.
    for(r in 1:R) {
        matrix[M[D+r],M[D+r]] cov
            = tcrossprod(to_matrix(segment(mm, sumMMplus[D+r] + 1, MMplus[D+r]),
                                   M[D+r],
                                   Mplus[D+r] - M[D+r]));
        for(n in 1:(N+nVarGroups)) {
            int segplace = 1;
            sumIR[r,n] = sum(IR[(sumM[D+r]-V+1):(sumM[D+r+1]-V),n]);
            for(o in 1:M[D+r]) {
                if(IR[sumM[D+r]-V+o, n]) {
                    IRInds[r,n,segplace] = o;
                    segplace += 1;
                }
            }
            for(o in 1:sumIR[r,n]) {
                if(sum(cov[IRInds[r,n,1:sumIR[r,n]],IRInds[r,n,o]]) == cov[IRInds[r,n,o],IRInds[r,n,o]]) {
                    nIsolate[r,n] += 1;
                    segInds1[r,n,nIsolate[r,n]] = o;
                }
            }
        }
        uniqueRInds[r,1,] = IRInds[r,1,];
        nMatches[r,1] = 1;
        matchInds[r,1] = 1;
        matchIndsInv[r,1,1] = 1;
        for(n in 2:(N+nVarGroups)) {
            for(m in 1:nUniqueR[r]) {
                if(sumIR[r,n] == sumIR[r,matchIndsInv[r,m,1]]) {
                    for(o in 1:sumIR[r,n]) {
                        if(IRInds[r,n,o] != uniqueRInds[r,m,o]) {
                            break;
                        }
                        if(o == sumIR[r,n]) {
                            nMatches[r,m] += 1;
                            matchInds[r,n] = m;
                            matchIndsInv[r,m,nMatches[r,m]] = n;
                        }
                    }
                }
                if((matchInds[r,n] == 0) && (m == nUniqueR[r])) {
                    nUniqueR[r] += 1;
                    uniqueRInds[r,nUniqueR[r],] = IRInds[r,n,];
                    nMatches[r,nUniqueR[r]] += 1;
                    matchInds[r,n] = nUniqueR[r];
                    matchIndsInv[r,nUniqueR[r],nMatches[r,nUniqueR[r]]] = n;
                }
            }
        }
        for(m in 1:nUniqueR[r]) {
            sumIRUnique[r,m] = sumIR[r,matchIndsInv[r,m,1]];
            for(o in 1:sumIRUnique[r,m]) {
                if(sum(cov[uniqueRInds[r,m,1:sumIRUnique[r,m]],uniqueRInds[r,m,o]]) != cov[uniqueRInds[r,m,o],uniqueRInds[r,m,o]]) {
                    nMat[r,m] += 1;
                    segInds2[r,m,nMat[r,m]] = o;
                }
            }
        }
    }
    // end creating continuous variable indices
    // create indices for categorical variables
    sumMc[1] = 0;
    for(cv in 1:C_vars) sumMc[cv+1] = sum(Mc[1:cv]);
    //
    // pre-compute some parts of the circular matern covariance function
    for(r in 0:(site_smoothness-2)) {
        for(j in 0:(2*r+1)) {
            chooseRJ[r+1,j+1] = choose(2*r+1,j);
        }
    }
    for(k in 0:(site_smoothness-1)) {
        for(j in 0:k) {
            ffKJ[k+1,j+1] = ff(k,j);
        }
    }
    //
}
parameters {
    matrix[K_linear + KG * K_gp,N] Z;              // PCA axis scores
    matrix[VOBplus+Vplus+D,K] W_norm;              // PCA variable loadings
    vector<lower=0>[VOBplus+Vplus+D] sds;          // variable scales
    real<lower=0> global_effect_scale;             // overall scale of variable loadings
    real<lower=0> ortho_scale;                     // inverse strength of orthogonality shrinkage
    row_vector<lower=0>[K] latent_scales;          // overall scale of each axis
    vector<lower=0>[DRC+D] dataset_scales;         // overall scale of each dataset
    matrix<lower=0>[DRC+D,K] weight_scales;        // per-dataset-and-axis scales
    vector[VOB] intercepts;
    vector[V] binary_count_intercepts;
    vector[F] abundance_true_vector;               // count dataset latent log abundance
    vector[F_higher] abundance_higher_vector;
    vector[D] binary_count_dataset_intercepts;     // each count dataset could have biased prior intercepts
    vector[nmulti] multinomial_nuisance;           // per-sample parameter to convert independent poisson distributions to a multinomial one
    vector<upper=0>[nMP] P_missing;                // latent estimates of continuous variables with partial information (truncated observations)
    matrix<lower=0>[DRC+D,K] nu_factors_raw;       // per-dataset-and-axis sparsity of variable loadings
    vector<lower=0>[K] rho_sites;                  // length scale for site gaussian process
    vector<lower=0, upper=1>[K] site_prop;         // importance of site covariance compared to nugget effect
    matrix<lower=0>[K_linear,KG] rho_Z;            // length scale for latent axis gaussian process
    vector<upper=0>[D] inv_log_less_contamination; // smaller = less average contamination
    vector<lower=0>[D] contaminant_overdisp;       // dispersion parameter for amount of contamination in true negative count observations
}
transformed parameters {
    matrix[K,nVarGroups] Z_higher = Z * samp2group;
    matrix[K_linear + KG * K_gp, N] Z_ortho = diag_pre_multiply(sqrt(rows_dot_self(Z)), svd_V(Z')) * svd_U(Z')';
    matrix[VOBplus+Vplus+D,K] W_ortho = svd_U(W_norm) * diag_post_multiply(svd_V(W_norm)', sqrt(columns_dot_self(W_norm)));
    matrix[VOB,K] W;
    matrix[V,K] W_binary_counts;
    vector<lower=0>[VOBplus+Vplus+D] var_scales
        = sds
          .* append_row(prior_scales,
                        ones_vector(Vplus+D));
    vector<upper=0>[D] log_less_contamination = inv(inv_log_less_contamination);
    vector[H] P_filled = P;
    cov_matrix[M[DRC]] cov_sites[K];
    matrix<lower=2>[DRC,K] nu_factors = nu_factors_raw + 2;
    for(miss in 1:nMP) P_filled[indMP[miss]] = P_missing[miss] + P_max[miss];
    for(k in 1:K) {
        cov_sites[k]
            = fill_sym(site_prop[k]
                       * square(weight_scales[DRC,k])
                   //  * circular_matern(dist_sites, site_smoothness, inv(rho_sites[k]), ffKJ, chooseRJ),
                       * exp(-dist_sites / rho_sites[k]),
                       M[DRC],
                       square(weight_scales[DRC,k]) + 1e-9);
    } // determine covariance of sites
    for(drc in 1:DRC) {
        W[(sumM[drc] + 1):(sumM[drc] + M[drc]),]
            = diag_pre_multiply(segment(var_scales, sumMplus[drc] + 1, M[drc]),
                                W_norm[(sumMplus[drc] + 1):(sumMplus[drc] + M[drc]),]); // effects for each base variable
        if(drc == DRC) {
            for(k in 1:K) {
                W[(sumM[drc] + 1):(sumM[drc] + M[drc]),k]
                    = cholesky_decompose(cov_sites[k])
                      * sqrt(nu_factors_raw[DRC,k] / nu_factors[DRC,k])
                      * W[(sumM[drc] + 1):(sumM[drc] + M[drc]),k];
            }
        } // induce correlation among sites
        if(Mplus[drc] > M[drc]) {
            W[(sumM[drc] + 1):(sumM[drc] + M[drc]),]
                += to_matrix(segment(mm, sumMMplus[drc] + 1, MMplus[drc]), M[drc], Mplus[drc] - M[drc])
                   * diag_pre_multiply(segment(var_scales, sumMplus[drc] + M[drc] + 1, Mplus[drc] - M[drc]),
                                       W_norm[(sumMplus[drc] + M[drc] + 1):(sumMplus[drc] + Mplus[drc]),]); // add higher level effects
        } // add higher level effects
    } // scale PCA factor loadings and combine linear effects
    for(d in 1:D) {
        W_binary_counts[(sumM[d] + 1):(sumM[d] + M[d]),]
            = diag_pre_multiply(segment(var_scales, VOBplus + sumMplus[d] + 1, M[d]),
                                W_norm[(VOBplus + sumMplus[d] + 1):(VOBplus + sumMplus[d] + M[d]),]);
        if(Mplus[d] > M[d]) {
            W_binary_counts[(sumM[d] + 1):(sumM[d] + M[d]),]
                += to_matrix(segment(mm, sumMMplus[d] + 1, MMplus[d]), M[d], Mplus[d] - M[d])
                   * diag_pre_multiply(segment(var_scales, VOBplus + sumMplus[d] + M[d] + 1, Mplus[d] - M[d]),
                                       W_norm[(VOBplus + sumMplus[d] + M[d] + 1):(VOBplus + sumMplus[d] + Mplus[d]),]);
        }
        W_binary_counts[(sumM[d] + 1):(sumM[d] + M[d]),]
            += rep_matrix(var_scales[VOBplus+Vplus+d] * W_norm[VOBplus+Vplus+d,], M[d]);
    } // same as above, for prevalence estimates in count data
}
model {
    // data wrangling (should replace some with transformed data indices)
    matrix[K, N + nVarGroups] Z_Z_higher = append_col(Z,Z_higher);
    vector[VOBplus+Vplus+D] dsv;
    matrix[VOBplus+Vplus+D,K] num;
    matrix[VOBplus+Vplus+D,K] wsn;
    int Xplace = 1;
    int X_higher_place = 1;
    int Pplace = 1;
    int Yplace = 1;
    int multinomPlace = 1;
    for(drc in 1:DRC) {
        dsv[(sumMplus[drc] + 1):(sumMplus[drc] + Mplus[drc])] = rep_vector(dataset_scales[drc], Mplus[drc]);
        for(k in 1:K) {
            num[(sumMplus[drc] + 1):(sumMplus[drc] + Mplus[drc]),k]
                = rep_vector(nu_factors[drc,k],
                             Mplus[drc]);
            wsn[(sumMplus[drc] + 1):(sumMplus[drc] + Mplus[drc]),k]
                = rep_vector(weight_scales[drc,k]
                             * sqrt(nu_factors_raw[drc,k] / nu_factors[drc,k]),
                             Mplus[drc]);
        }
        if(drc <= D) {
            dsv[(VOBplus + sumMplus[drc] + 1):(VOBplus + sumMplus[drc] + Mplus[drc])] = rep_vector(dataset_scales[DRC+drc], Mplus[drc]);
            dsv[VOBplus + Vplus + drc] = dataset_scales[drc];
            for(k in 1:K) {
                num[(VOBplus + sumMplus[drc] + 1):(VOBplus + sumMplus[drc] + Mplus[drc]),k]
                    = rep_vector(nu_factors[drc,k],
                                 Mplus[drc]);
                num[VOBplus+Vplus+drc,k]
                    = nu_factors[DRC+drc,k];
                wsn[(VOBplus + sumMplus[drc] + 1):(VOBplus + sumMplus[drc] + Mplus[drc]),k]
                    = rep_vector(weight_scales[drc,k]
                                 * sqrt(nu_factors_raw[drc,k] / nu_factors[drc,k]),
                                 Mplus[drc]);
                wsn[VOBplus+Vplus+drc,k]
                    = weight_scales[DRC+drc,k]
                      * sqrt(nu_factors_raw[DRC+drc,k] / nu_factors[DRC+drc,k]);
            }
        }
    }
    // end data wrangling
    // priors
    target += inv_gamma_lupdf(to_vector(nu_factors_raw) | shape_gamma_fact, rate_gamma_fact);                // PCA variable loadings have dataset and axis specific sparsity
    target += std_normal_lupdf(dataset_scales);                                                              // entire datasets may have uniformly biased difference in scale compared to priors
    target += cauchy_lupdf(intercepts | prior_intercept_centers, 2.5 * prior_intercept_scales);              // overall abundance and center of variables may differ from priors
    target += cauchy_lupdf(binary_count_intercepts | binary_count_intercept_centers, 2.5);                   // overall prevalence may differ from priors
    target += cauchy_lupdf(binary_count_dataset_intercepts | 0, 2.5);                                        // allow entire count datasets to have biased difference in prevalence compared to priors
    target += normal_lupdf(global_effect_scale | 0, global_scale_prior);                                     // shrink global scale of effects toward zero
    target += normal_lupdf(ortho_scale | 0, ortho_scale_prior);                                              // estimate necessary strength of orthonogonalization
    target += student_t_lupdf(latent_scales | 2, 0, global_effect_scale);                                    // sparse selection of number of axes
    target += student_t_lupdf(to_vector(weight_scales) | 2, 0, to_vector(rep_matrix(latent_scales, DRC+D))); // sparse selection of datasets per axis
    target += generalized_normal_lpdf(inv_log_less_contamination | 0, inv_log_max_contam, gnorm_shape);      // shrink amount of contamination in 'true zeros' toward zero
    target += std_normal_lupdf(contaminant_overdisp);                                                        // shrink overdispersion of contaminant counts in 'true zeros' toward zero
    target += normal_lupdf(to_vector(W_norm) | to_vector(W_ortho), global_effect_scale * ortho_scale);       // shrink PCA variable loadings toward closes orthogonal matrix
    target += normal_lupdf(to_vector(Z) | to_vector(Z_ortho), ortho_scale);                                  // shrink PCA axis scores toward closes orthogonal matrix
    target += std_normal_lupdf(to_vector(Z[1:K_linear,]));                                                   // first PCA axis scores are independent of one another
    target += inv_gamma_lupdf(rho_sites | 5, 5 * rho_sites_prior);                                           // length scale for gaussian process on sites
    target += inv_gamma_lupdf(to_vector(rho_Z) | rho_Z_shape, rho_Z_scale);                                  // length scale for gaussian process on PCA axis scores
    for(g in 1:KG) {
        target += multi_gp_cholesky_lupdf(Z[(K_linear + (K_gp * (g-1)) + 1):(K_linear + K_gp * g),] |
                                          L_cov_exp_quad_ARD(Z[1:K_linear,], rho_Z[,g], 1e-9),
                                          ones_vector(K_gp));
    }                                                                                                        // final KG * K_gp PCA axis scores are functions of first K_linear ones
    target += normal_lupdf(sds | 0, dsv);                                                                    // per-variable sigmas shrink toward dataset scales
    target += student_t_lupdf(to_vector(W_norm) | to_vector(num), 0, to_vector(wsn));                        // PCA loadings shrink to zero with axis-and-dataset-specific nu and variance
    // end priors
    // likelihoods
    for(d in 1:D) {
        matrix[M[d],sumID[d]] predicted
            = rep_matrix(intercepts[(sumM[d] + 1):(sumM[d] + M[d])], sumID[d])
              + W[(sumM[d] + 1):(sumM[d] + M[d]),]
                * Z_Z_higher[,IDInds[d,1:sumID[d]]];
        matrix[M[d],sumID[d]] prevalence
            = rep_matrix(binary_count_dataset_intercepts[d]
                         + segment(binary_count_intercepts, sumM[d] + 1, M[d]),
                         sumID[d])
              + W_binary_counts[(sumM[d] + 1):(sumM[d] + M[d]),]
                * Z_Z_higher[,IDInds[d,1:sumID[d]]];
        vector[M[d]] abundance_contam
            = segment(intercepts, sumM[d] + 1, M[d])
              + log_inv_logit(binary_count_dataset_intercepts[d]
                              + segment(binary_count_intercepts, sumM[d] + 1, M[d]))
              + log_less_contamination[d];
        matrix[M[d],sumID[d]] abundance_true
            = to_matrix(segment(abundance_true_vector, Xplace, M[d] * sumID[d]),
                        M[d], sumID[d]);
        vector[M[d]] phi;
        if(Mplus[d] > M[d]) {
            matrix[M[d],Mplus[d] - M[d]] MM
                = to_matrix(segment(mm, sumMMplus[d] + 1, MMplus[d]),
                            M[d],
                            Mplus[d] - M[d]);
            matrix[M[d],sumID[d]] resid_higher
                = MM
                  * to_matrix(segment(abundance_higher_vector, X_higher_place, (Mplus[d] - M[d]) * sumID[d]),
                              Mplus[d] - M[d], sumID[d]);
            target += student_t_lpdf(segment(abundance_higher_vector, X_higher_place, (Mplus[d] - M[d]) * sumID[d]) |
                                     nu_residuals,
                                     0,
                                     to_vector(rep_matrix(segment(var_scales, sumMplus[d] + M[d] + 1, Mplus[d] - M[d]), sumID[d])));
            X_higher_place += (Mplus[d] - M[d]) * sumID[d];
            target += student_t_lpdf(to_vector(abundance_true) |
                                     nu_residuals,
                                     to_vector(predicted + resid_higher),
                                     to_vector(rep_matrix(segment(var_scales, sumMplus[d] + 1, M[d]), sumID[d])));
            phi
                = inv_square(contaminant_overdisp[d])
                  * inv(rows_dot_self(diag_post_multiply(MM,
                                                         segment(var_scales, sumMplus[d] + M[d] + 1, Mplus[d] - M[d])))
                        + square(segment(var_scales, sumMplus[d] + 1, M[d])));
        } else {
            target += student_t_lupdf(to_vector(abundance_true) |
                                      nu_residuals,
                                      to_vector(predicted),
                                      to_vector(rep_matrix(segment(var_scales, sumMplus[d] + 1, M[d]), sumID[d])));
            phi = inv_square(contaminant_overdisp[d] * segment(var_scales, sumMplus[d] + 1, M[d]));
        }
        for(n in 1:sumID[d]) {
            for(m in 1:M[d]) {
                target += log_sum_exp(log1m_inv_logit(prevalence[m,n])
                                      + neg_binomial_2_log_lpmf(X[Xplace + m - 1] |
                                                                abundance_contam[m] + multinomial_nuisance[multinomPlace],
                                                                phi[m]), //estimated abundance if true negative
                                        log_inv_logit(prevalence[m,n])
                                        + poisson_log_lpmf(X[Xplace + m - 1] |
                                                           log_sum_exp(abundance_contam[m], abundance_true[m,n])
                                                           + multinomial_nuisance[multinomPlace])); //estimated abundance if true positive
            }
            multinomPlace += 1;
            Xplace += M[d];
        }
    } // count likelihood
    for(r in 1:R) {
        if(Mplus[D+r] > M[D+r]) {
            matrix[max(nMat[r,]), N+nVarGroups] observed;
            matrix[M[D+r],M[D+r]] cov
                = add_diag(tcrossprod(diag_post_multiply(
                                   to_matrix(segment(mm, sumMMplus[D+r] + 1, MMplus[D+r]),
                                             M[D+r],
                                             Mplus[D+r] - M[D+r]),
                                   segment(var_scales, sumMplus[D+r] + M[D+r] + 1, Mplus[D+r] - M[D+r]))),
                           square(segment(var_scales, sumMplus[D+r] + 1, M[D+r])) + 1e-9);
            for(n in 1:(N+nVarGroups)) {
                if(nIsolate[r,n] > 0) {
                    int inds[nIsolate[r,n]] = IRInds[r,n,1:sumIR[r,n]][segInds1[r,n,1:nIsolate[r,n]]];
                    target += student_t_lupdf(segment(P_filled, Pplace, sumIR[r,n])[segInds1[r,n,1:nIsolate[r,n]]] |
                                              nu_residuals,
                                              intercepts[(sumM[D+r] + 1):sumM[D+r+1]][inds]
                                              + W[(sumM[D+r] + 1):sumM[D+r+1],][inds,]
                                                * Z_Z_higher[,n],
                                              sqrt(diagonal(cov)[inds]));
                } // likelihood for variables with no higher level covariance, within datasets that do have such effects
                if(nMat[r,matchInds[r,n]] > 0) {
                    int inds[nMat[r,matchInds[r,n]]] = segInds2[r,matchInds[r,n],1:nMat[r,matchInds[r,n]]];
                    observed[1:nMat[r,matchInds[r,n]],n] = segment(P_filled, Pplace, sumIR[r,n])[inds];
                }
                Pplace += sumIR[r,n];
            }
            for(m in 1:nUniqueR[r]) {
                if(nMat[r,m] > 0) {
                    int inds[nMat[r,m]] = uniqueRInds[r,m,1:sumIRUnique[r,m]][segInds2[r,m,1:nMat[r,m]]];
                    matrix[nMat[r,m], nMatches[r,m]] predicted
                        = rep_matrix(intercepts[(sumM[D+r] + 1):(sumM[D+r+1])][inds], nMatches[r,m])
                          + W[(sumM[D+r] + 1):(sumM[D+r] + M[D+r]),][inds,]
                            * Z_Z_higher[,matchIndsInv[r,m,1:nMatches[r,m]]];
                    target += multi_student_t_cholesky_lpdf(observed[1:nMat[r,m], matchIndsInv[r,m,1:nMatches[r,m]]] |
                                                            nu_residuals,
                                                            predicted,
                                                           cholesky_decompose(cov[inds,inds]));
                } // likelihood for sets of variables that share a subsettable covariance matrix
            }
        } else {
            for(n in 1:(N+nVarGroups)) {
                if(sumIR[r,n] > 0) {
                    target += student_t_lupdf(segment(P_filled, Pplace, sumIR[r,n]) |
                                              nu_residuals,
                                              intercepts[(sumM[D+r] + 1):sumM[D+r+1]][IRInds[r,n,1:sumIR[r,n]]]
                                              + W[(sumM[D+r] + 1):sumM[D+r+1],][IRInds[r,n,1:sumIR[r,n]],]
                                                * Z_Z_higher[,n],
                                              segment(var_scales, sumM[D+r] + 1, M[D+r])[IRInds[r,n,1:sumIR[r,n]]]);
                    Pplace += sumIR[r,n];
                }
            }
        } // likelihood for entire datasets with no higher level covariance
    } // continuous likelihood
    for(cv in 1:C_vars) {
        for(n in 1:(N+nVarGroups)) {
            if(IC[cv,n]) {
                if(Mc[cv] > 1) {
                    vector[Mc[cv]] predicted
                        = log_softmax(segment(intercepts, V + O + sumMc[cv] + 1, Mc[cv])
                                      + W[(V + O + sumMc[cv] + 1):(V + O + sumMc[cv] + Mc[cv]),]
                                        * Z_Z_higher[,n]);
                    vector[Mc[cv]] observed = segment(Y, Yplace, Mc[cv]);
                    real terms = 0;
                    for(c in 1:Mc[cv]) {
                        if(observed[c] > 0) {
                            if(terms == 0) {
                                terms = log(observed[c]) + predicted[c];
                            } else {
                                terms = log_sum_exp(terms, log(observed[c]) + predicted[c]);
                            }
                        }
                    }
                    target += terms; // likelihood for categorical variables
                } else {
                    real predicted
                        = intercepts[V + O + sumMc[cv] + 1]
                          + W[V + O + sumMc[cv] + 1,]
                            * Z_Z_higher[,n];
                    if(Y[Yplace] == 1) {
                        target += log_inv_logit(predicted);
                    } else if(Y[Yplace] == 0) {
                        target += log1m_inv_logit(predicted);
                    } else {
                        target += log_mix(Y[Yplace],
                                          log_inv_logit(predicted),
                                          log1m_inv_logit(predicted));
                    }
                } // likelihood for binary variables
                Yplace += Mc[cv];
            }
        }
    } // categorical/binary likelihood allowing uncertain (real-valued) data
}
