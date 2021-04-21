//strong inspiration from https://www.cs.helsinki.fi/u/sakaya/tutorial/code/cca.R
#include functions.stan
data {
    int<lower=0> N; // Number of samples
    int<lower=0> D; // Number of compositional count datasets
    int<lower=0> R; // Number of continuous datasets
    int<lower=0> C; // Number of categorical datasets
    int<lower=0> nVarGroups;
    int<lower=0> M[D+R+C]; // Number of observed variables in each dataset
    int<lower=0> Mplus[D+R+C]; // Number of parameters for each dataset
    int<lower=0> K_linear; // Number of residual latent linear dimensions
    int<lower=0> K_gp; // Number of residual latent GP dimensions per kernel
    int<lower=0> KG; // Number of unique GP kernels
    int<lower=0> I[D, N + nVarGroups]; // Sample indices for count datasets
    int<lower=0> O; // Total number of observed continuous variables
    int<lower=0> IR[O, N + nVarGroups]; // Sample indices for continuous variables
    int<lower=0> C_vars; // Total number of observed categorical variables
    int<lower=0> IC[C_vars, N + nVarGroups]; // Sample indices for categorical datasets
    int<lower=0> Mc[C_vars]; // Number of levels for each categorical variable
    int<lower=0> F; // Total number of compositional count observations
    int<lower=0> X[F]; // compositional count data
    int<lower=0> G; // Total number of categorical observations
    vector<lower=0>[G] Y; // categorical data
    int<lower=0> H; // Total number of continuous observations
    vector[H] P; // continuous data
    real<lower=0> global_scale_prior;
    real<lower=0> ortho_scale_prior;
    vector<lower=0>[sum(Mplus)] prior_scales;
    vector<lower=0>[sum(M)] prior_intercept_scales;
    vector[sum(M)] prior_intercept_centers;
    vector[sum(M[1:D])] binary_count_intercept_centers;
    int<lower=0> sizeMM;
    vector[sizeMM] mm; // model matrix for higher level variables
    int<lower=0> nMP; // number of P measured as -inf
    int<lower=0> indMP[nMP]; // indices for missing Ps
    real P_max[nMP]; // lowest measured value for the variable corresponding to a missing P
    real rate_gamma_fact;
    real shape_gamma_fact;
    vector[choose(M[size(M)],2)] dist_sites;
    real rho_sites_prior;
    int nVarsWGroups;
    matrix[N,nVarGroups] samp2group;
    int varsWGroupsInds[nVarsWGroups];
    int ms; // Matern smoothness parameter for sites
    real nu_residuals;
    vector[D] inv_log_max_contam; // prior expectation of contamination rate
    real<lower=0> gnorm_shape;
}
transformed data {
    int K = K_linear + KG * K_gp;
    int DRC = D+R+C; // Total number of datasets
    int V = 0; // Total number of observed compositional count variables
    int Vplus = 0; // Total number of observed compositional count parameters
    int VOB = sum(M); // Total number of observed variables
    int VOBplus = sum(Mplus); // Total number of parameters per dim
    int sumM[DRC+1]; // Cumulative number of base level variables
    int sumMplus[DRC+1]; // Cumulative number of base level and higher variables
    int MMplus[DRC]; // size of model matrix for higher level variables
    int sumMMplus[DRC+1]; // Cumulative size of model matrices for higher level variables
    int sumMc[C_vars+1]; // Cumulative number of categorical levels
    real rho_Z_shape = 1.0 / K_linear;
    real rho_Z_scale = 6.0 / K_gp; // * 2 * tgamma((K_linear+1)/2.0) / tgamma(K_linear/2.0);
    int nmulti = 0;
    int sumID[D];
    int IDInds[D,N+nVarGroups];
    int X0Inds[D,N+nVarGroups];
    int X1Inds[D,N+nVarGroups];
    int NX1[D,N+nVarGroups] = rep_array(0,D,N+nVarGroups);
    int NX0[D,N+nVarGroups] = rep_array(0,D,N+nVarGroups);
    int xplace = 0;
    int sumIR[R,N+nVarGroups];
    int IRInds[R,N+nVarGroups,max(M[(D+1):(D+R)])];
    int segInds1[R,N+nVarGroups,max(M[(D+1):(D+R)])];
    int segInds2[R,N+nVarGroups,max(M[(D+1):(D+R)])];
    int nIsolate[R,N+nVarGroups] = rep_array(0,R,N+nVarGroups);
    int nMat[R,N+nVarGroups] = rep_array(0,R,N+nVarGroups);
    int matchInds[R,N+nVarGroups] = rep_array(0,R,N+nVarGroups);
    int matchIndsInv[R,N+nVarGroups,N+nVarGroups] = rep_array(0,R,N+nVarGroups,N+nVarGroups);
    int nMatches[R,N+nVarGroups] = rep_array(0,R,N+nVarGroups);
    int nUniqueR[R] = rep_array(1,R);
    int uniqueRInds[R,N+nVarGroups,max(M[(D+1):(D+R)])] = rep_array(0,R,N+nVarGroups,max(M[(D+1):(D+R)]));
    int sumIRUnique[R,N+nVarGroups];
    matrix[ms-1, 2 * ms] chooseRJ;
    matrix[ms, 2 * ms] ffKJ;
    sumM[1] = 0;
    sumMplus[1] = 0;
    sumMMplus[1] = 0;
    for(drc in 1:DRC) {
        sumM[drc+1] = sum(M[1:drc]);
        sumMplus[drc+1] = sum(Mplus[1:drc]);
        MMplus[drc] = M[drc] * (Mplus[drc] - M[drc]);
        sumMMplus[drc+1] = sum(MMplus[1:drc]);
    }
    V = sumM[D+1];
    Vplus = sumMplus[D+1];
    for(d in 1:D) {
        int nplace = 1;
        sumID[d] = sum(I[d,]);
        for(n in 1:(N+nVarGroups)) {
            if(I[d,n]) {
                IDInds[d,nplace] = n;
                nplace += 1;
                for(m in 1:M[d]) {
                    if(X[xplace + m]) {
                        NX1[d,n] += 1;
                        X1Inds[d,NX1[d,n]] = m;
                    } else {
                        NX0[d,n] += 1;
                        X0Inds[d,NX0[d,n]] = m;
                    }
                }
                xplace += M[d];
            }
        }
    }
    nmulti = sum(sumID);
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
    sumMc[1] = 0;
    for(cv in 1:C_vars) sumMc[cv+1] = sum(Mc[1:cv]);
    for(r in 0:(ms-2)) {
        for(j in 0:(2*r+1)) {
            chooseRJ[r+1,j+1] = choose(2*r+1,j);
        }
    }
    for(k in 0:(ms-1)) {
        for(j in 0:k) {
            ffKJ[k+1,j+1] = ff(k,j);
        }
    }
}
parameters {
    matrix[K_linear + KG * K_gp,N] Z; // PCA axis scores
    matrix[VOBplus+Vplus+D,K] W_norm; // PCA variable loadings
    vector<lower=0>[VOBplus+Vplus+D] sds; // variable scales
    real<lower=0> global_effect_scale;
    real<lower=0> ortho_scale;
    row_vector<lower=0>[K] latent_scales;
    vector<lower=0>[DRC+D] dataset_scales;
    matrix<lower=0>[DRC+D,K] weight_scales;
    vector[VOB] intercepts;
    vector[V] binary_count_intercepts;
    vector[F] abundance_true_vector;
    vector[D] binary_count_dataset_intercepts;
    vector[nmulti] multinomial_nuisance;
    vector<upper=0>[nMP] P_missing;
    matrix<lower=0>[DRC+D,K] nu_factors_raw;
    vector<lower=0>[K] rho_sites;
    vector<lower=0, upper=1>[K] site_prop;
    matrix<lower=0>[K_linear,KG] rho_Z;
    vector<upper=0>[D] inv_log_less_contamination;
    vector<lower=0>[D] contaminant_overdisp;
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
//                     * circular_matern(dist_sites, ms, inv(rho_sites[k]), ffKJ, chooseRJ),
                       * exp(-dist_sites / rho_sites[k]),
                       M[DRC],
                       square(weight_scales[DRC,k]) + 1e-9);
    }
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
    matrix[K, N + nVarGroups] Z_Z_higher = append_col(Z,Z_higher);
    vector[VOBplus+Vplus+D] dsv;
    matrix[VOBplus+Vplus+D,K] num;
    matrix[VOBplus+Vplus+D,K] wsn;
    int Xplace = 1;
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
            dsv[(VOBplus + sumMplus[d] + 1):(VOBplus + sumMplus[d] + Mplus[d])] = rep_vector(dataset_scales[DRC+d], Mplus[DRC+d]);
            dsv[VOBplus + Vplus + d] = dataset_scales[d];
            for(k in 1:K) {
                num[(VOBplus + sumMplus[d] + 1):(VOBplus + sumMplus[d] + Mplus[d]),k]
                    = rep_vector(nu_factors[d,k],
                                 Mplus[d]);
                num[VOBplus+Vplus+d,k]
                    = rep_vector(nu_factors[DRC+d,k],
                                 Mplus[d]);
                wsn[(VOBplus + sumMplus[d] + 1):(VOBplus + sumMplus[d] + Mplus[d]),k]
                    = rep_vector(weight_scales[d,k]
                                 * sqrt(nu_factors_raw[d,k] / nu_factors[d,k]),
                                 Mplus[d]);
                wsn[VOBplus+Vplus+d,k]
                    = rep_vector(weight_scales[DRC+d,k]
                                 * sqrt(nu_factors_raw[DRC+d,k] / nu_factors[DRC+d,k]),
                                 Mplus[d]);
            }
        }
    } // data wrangling (should replace with transformed data indices)
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
            matrix[M[d],M[d]] cov
                = add_diag(tcrossprod(diag_post_multiply(
                                   to_matrix(segment(mm, sumMMplus[d] + 1, MMplus[d]),
                                             M[d],
                                             Mplus[d] - M[d]),
                                   segment(var_scales, sumMplus[d] + M[d] + 1, Mplus[d] - M[d]))),
                           square(segment(var_scales, sumMplus[d] + 1, M[d])) + 1e-9);
            target += multi_student_t_cholesky_lpdf(abundance_true |
                                                    nu_residuals,
                                                    predicted,
                                                    cholesky_decompose(cov));
            phi = inv_square(contaminant_overdisp[d]) * inv(diagonal(cov));
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
                }
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
                }
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
        }
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
                    target += terms;
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
                }
                Yplace += Mc[cv];
            }
        }
    } // categorical/binary likelihood allowing uncertain (real-valued) data
}
