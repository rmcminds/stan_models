functions {
    real generalized_normal_lpdf(vector y, real mu, vector alpha, real beta) {
        return sum(log(beta) - log(2) - log(alpha) - lgamma(inv(beta)) - exp(beta * log(fabs(y-mu)./alpha)));
    }
    real ff(int k, int j) {
        if(j)
          return(falling_factorial(k,j));
        else
          return(1);
    }
    vector circular_matern(vector d, int n, real alpha, matrix ffKJ, matrix chooseRJ) {
        real ap = alpha * pi();
        row_vector[2] csap = [cosh(ap), sinh(ap)];
        matrix[n-1,n] H;
        real annm1 = inv((-2 * square(alpha))^(n-1) * tgamma(n));
        vector[n] a;
        vector[rows(d)] adp = alpha * (d-pi());
        matrix[rows(d),2] csadp = append_col(cosh(adp), sinh(adp));
        vector[rows(d)] cov = zeros_vector(rows(d));
        for(k in 0:(n-1)) {
          for(r in 0:(n-2)) {
            H[r+1,k+1] = 0;
            for(j in 0:(2*r+1)) {
              if(j <= k) {
                H[r+1,k+1]
                 += chooseRJ[r+1,j+1]
                    * ffKJ[k+1,j+1]
                    * ap^(k-j)
                    * csap[((k-j+1) % 2) + 1];
              }
            }
          }
        }
        a = append_row(-annm1 * (H[,1:(n-1)] \ H[,n-1]), annm1);
        for(k in 0:(n-1)) {
          cov += a[k+1] * adp^k .* csadp[,(k % 2) + 1];
        }
        return(cov / (2*alpha*csap[2]));
    }
    matrix fill_sym(vector lt, int N, real c) {
        matrix[N,N] s_mat;
        int iter = 1;
        for(j in 1:(N-1)) {
            s_mat[j,j] = c;
            for(i in (j+1):N) {
                s_mat[i,j] = lt[iter];
                s_mat[j,i] = lt[iter];
                iter += 1;
            }
        }
        s_mat[N,N] = c;
        return(s_mat);
    }
    matrix L_cov_exp_quad_ARD(matrix x, vector rho, real delta) {
        int N = cols(x);
        matrix[N,N] cov;
        for (i in 1:(N-1)) {
            cov[i,i] = 1 + delta;
            for (j in (i+1):N) {
                cov[j,i] = exp(-0.5 * dot_self((x[,i] - x[,j]) ./ rho));
            }
        }
        cov[N,N] = 1 + delta;
        return cholesky_decompose(symmetrize_from_lower_tri(cov));
    }
}
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
    vector<lower=0>[sum(Mplus)] priorScales;
    vector<lower=0>[sum(M)] priorInterceptScales;
    vector[sum(M)] priorInterceptCenters;
    vector[sum(M[1:D])] binary_count_intercept_centers;
    int<lower=0> sizeMM;
    vector[sizeMM] mm; // model matrix for higher level variables
    int<lower=0> nMP; // number of P measured as -inf
    int<lower=0> indMP[nMP]; // indices for missing Ps
    real pMax[nMP]; // lowest measured value for the variable corresponding to a missing P
    real gammaRateFact;
    real gammaShapeFact;
    vector[choose(M[size(M)],2)] distSites;
    real rhoSitesPrior;
    int nVarsWGroups;
    matrix[N,nVarGroups] samp2group;
    int varsWGroupsInds[nVarsWGroups];
    int ms; // Matern smoothness parameter for sites
    real nu_residuals;
    vector[D] log_inv_max_contam; // prior expectation of contamination rate
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
    real meanDistZ = 2 * tgamma((K_linear+1)/2.0) / tgamma(K_linear/2.0);
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
    matrix[K_linear + KG * K_gp,N] Z_raw; // PCA axis scores
    matrix[VOBplus+Vplus+D,K] W_raw; // PCA variable loadings
    vector<lower=0>[VOBplus+Vplus+D] sds;
    real<lower=0> global_effect_scale;
    vector<lower=0>[K] latent_scales;
    vector<lower=0>[DRC+D] dataset_scales;
    matrix<lower=0>[DRC+D,K] weight_scales;
    vector[VOB] intercepts;
    vector[V] binary_count_intercepts;
    vector[F] latent_props_raw;
    vector[D] binary_count_dataset_intercepts;
    vector[nmulti] multinomial_nuisance;
    vector<upper=0>[nMP] missingP;
    matrix<lower=0>[DRC+D,K] nu_factors_raw;
    vector<lower=0>[K] rhoSites;
    vector<lower=0, upper=1>[K] siteProp;
    matrix<lower=0>[K_linear,KG] rhoZ;
    vector<upper=0>[D] inv_log_less_contamination;
    vector<lower=0>[D] contaminant_overDisp;
}
transformed parameters {
    vector[H] PFilled = P;
    vector<upper=0>[D] log_less_contamination = inv(inv_log_less_contamination);
    matrix[VOB,K] W;
    matrix[V,K] W_binary_counts;
    matrix[K,nVarGroups] Z_higher;
    vector<lower=0>[VOBplus+Vplus+D] var_scales
        = sds
          .* append_row(priorScales,
                        ones_vector(Vplus+D));
    matrix<lower=2>[DRC,K] nu_factors = nu_factors_raw + 2;
    cov_matrix[M[size(M)]] covSites[K];
    matrix[VOBplus+Vplus+D,K] W_norm = svd_U(W_raw) * diag_post_multiply(svd_V(W_raw)', sqrt(columns_dot_self(W_raw)));
    matrix[K_linear + KG * K_gp, N] Z = diag_pre_multiply(sqrt(rows_dot_self(Z_raw)), svd_V(Z_raw')) * svd_U(Z_raw')';
    for(miss in 1:nMP) PFilled[indMP[miss]] = missingP[miss] + pMax[miss];
    for(drc in 1:DRC) {
        W[(sumM[drc] + 1):(sumM[drc] + M[drc]),]
            = diag_pre_multiply(segment(var_scales, sumMplus[drc] + 1, M[drc]),
                                W_norm[(sumMplus[drc] + 1):(sumMplus[drc] + M[drc]),]); // effects for each base variable
        if(Mplus[drc] > M[drc]) {
            W[(sumM[drc] + 1):(sumM[drc] + M[drc]),]
                += to_matrix(segment(mm, sumMMplus[drc] + 1, MMplus[drc]), M[drc], Mplus[drc] - M[drc])
                   * diag_pre_multiply(segment(var_scales, sumMplus[drc] + M[drc] + 1, Mplus[drc] - M[drc]),
                                       W_norm[(sumMplus[drc] + M[drc] + 1):(sumMplus[drc] + Mplus[drc]),]); // add higher level effects
        }
    }
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
    }
    Z_higher = Z * samp2group;
    for(k in 1:K) {
        //covSites[k]
        //    = fill_sym(siteProp[k]
        //               * square(weight_scales[DRC,k])
        //               * circular_matern(distSites, ms, inv(rhoSites[k]), ffKJ, chooseRJ),
        //               M[size(M)],
        //               square(weight_scales[DRC,k]) + 1e-10);
        covSites[k]
            = fill_sym(siteProp[k]
                       * square(weight_scales[DRC,k])
                       * exp(-distSites / rhoSites[k]),
                       M[size(M)],
                       square(weight_scales[DRC,k]) + 1e-10);
    }
}
model {
    matrix[K, N + nVarGroups] Z_Z_higher = append_col(Z,Z_higher);
    int Xplace = 1;
    int Pplace = 1;
    int Yplace = 1;
    int multinomPlace = 1;
    target += inv_gamma_lupdf(to_vector(nu_factors_raw) | gammaShapeFact,gammaRateFact);
    target += std_normal_lupdf(dataset_scales);
    target += inv_gamma_lupdf(rhoSites | 5, 5 * rhoSitesPrior);
    target += cauchy_lupdf(intercepts | priorInterceptCenters, 2.5 * priorInterceptScales);
    target += cauchy_lupdf(binary_count_intercepts | binary_count_intercept_centers, 2.5);
    target += cauchy_lupdf(binary_count_dataset_intercepts | 0, 2.5);
    target += student_t_lupdf(global_effect_scale | 2, 0, global_scale_prior);
    target += student_t_lupdf(latent_scales | 2, 0, global_effect_scale);
    target += generalized_normal_lpdf(inv_log_less_contamination | 0, log_inv_max_contam, 15);
    target += std_normal_lupdf(contaminant_overDisp);
    target += std_normal_lupdf(to_vector(Z[1:K_linear,]));
    target += inv_gamma_lupdf(to_vector(rhoZ) | 1.0 / K_linear, 6.0 / K_gp);
    for(g in 1:KG) {
        target += multi_gp_cholesky_lupdf(Z[(K_linear + (K_gp * (g-1)) + 1):(K_linear + K_gp * g),] |
                                          L_cov_exp_quad_ARD(Z[1:K_linear,], rhoZ[,g], 1e-10),
                                          ones_vector(K_gp));
    }
    for(drc in 1:DRC) {
        target += normal_lupdf(sds[(sumMplus[drc] + 1):(sumMplus[drc] + Mplus[drc])] |
                               0,
                               dataset_scales[drc]);
    }
    for(d in 1:D) {
        target += normal_lupdf(sds[(VOBplus + sumMplus[d] + 1):(VOBplus + sumMplus[d] + Mplus[d])] |
                               0,
                               dataset_scales[DRC+d]);
        target += normal_lupdf(sds[VOBplus + Vplus + d] |
                               0,
                               dataset_scales[DRC+d]);
    }
    for(k in 1:K) {
        target += student_t_lupdf(weight_scales[,k] | 2, 0, latent_scales[k]);
        for(drc in 1:(DRC-1)) {
            target += student_t_lupdf(W_norm[(sumMplus[drc] + 1):(sumMplus[drc] + Mplus[drc]),k] |
                                      nu_factors[drc,k],
                                      0,
                                      weight_scales[drc,k]
                                      * sqrt(nu_factors_raw[drc,k] / nu_factors[drc,k]));
        }
        target += multi_student_t_lupdf(W_norm[(sumMplus[DRC] + 1):(sumMplus[DRC] + M[DRC]),k] |
                                        nu_factors[DRC,k],
                                        rep_vector(0,M[DRC]),
                                        covSites[k]
                                        * sqrt(nu_factors_raw[DRC,k] / nu_factors[DRC,k]));
        target += student_t_lupdf(W_norm[(sumMplus[DRC] + M[DRC] + 1):(sumMplus[DRC] + Mplus[DRC]),k] |
                                  nu_factors[DRC,k],
                                  0,
                                  weight_scales[DRC,k]
                                  * sqrt(nu_factors_raw[DRC,k] / nu_factors[DRC,k]));
        for(d in 1:D) {
            target += student_t_lupdf(W_norm[(VOBplus + sumMplus[d] + 1):(VOBplus + sumMplus[d] + Mplus[d]),k] |
                                      nu_factors[DRC+d,k],
                                      0,
                                      weight_scales[DRC+d,k]
                                      * sqrt(nu_factors_raw[DRC+d,k] / nu_factors[DRC+d,k]));
            target += student_t_lupdf(W_norm[VOBplus+Vplus+d,k] |
                                      nu_factors[DRC+d,k],
                                      0,
                                      weight_scales[DRC+d,k]
                                      * sqrt(nu_factors_raw[DRC+d,k] / nu_factors[DRC+d,k]));
        }
    }
    for(d in 1:D) {
        matrix[M[d],sumID[d]] predictTemp
            = rep_matrix(intercepts[(sumM[d] + 1):(sumM[d] + M[d])], sumID[d])
              + W[(sumM[d] + 1):(sumM[d] + M[d]),] * Z_Z_higher[,IDInds[d,1:sumID[d]]];
        matrix[M[d],sumID[d]] logit_prob_present
            = rep_matrix(binary_count_dataset_intercepts[d]
                         + segment(binary_count_intercepts, sumM[d] + 1, M[d]),
                         sumID[d])
              + W_binary_counts[(sumM[d] + 1):(sumM[d] + M[d]),]
                * Z_Z_higher[,IDInds[d,1:sumID[d]]];
        vector[M[d]] contamination
            = segment(intercepts, sumM[d] + 1, M[d])
              + log_inv_logit(binary_count_dataset_intercepts[d]
                              + segment(binary_count_intercepts, sumM[d] + 1, M[d]))
              + log_less_contamination[d];
        if(Mplus[d] > M[d]) {
            vector[M[d]] noise[sumID[d]];
            vector[M[d]] predicted[sumID[d]];
            int nObs = 1;
            matrix[M[d],M[d]] cov
                = add_diag(tcrossprod(diag_post_multiply(
                                   to_matrix(segment(mm, sumMMplus[d] + 1, MMplus[d]),
                                             M[d],
                                             Mplus[d] - M[d]),
                                   segment(var_scales, sumMplus[d] + M[d] + 1, Mplus[d] - M[d]))),
                           square(segment(var_scales, sumMplus[d] + 1, M[d])) + 1e-10);
            vector[M[d]] phi = inv_square(contaminant_overDisp[d]) * inv(diagonal(cov));
            for(n in 1:sumID[d]) {
                noise[n] = segment(latent_props_raw, Xplace, M[d]);
                predicted[n] = predictTemp[,n];
                for(m in 1:M[d]) {
                    target += log_sum_exp(log1m_inv_logit(logit_prob_present[m,n])
                                          + neg_binomial_2_log_lpmf(X[Xplace + m - 1] |
                                                                    contamination[m] + multinomial_nuisance[multinomPlace],
                                                                    phi[m]), //estimated abundance if true negative
                                            log_inv_logit(logit_prob_present[m,n])
                                            + poisson_log_lpmf(X[Xplace + m - 1] |
                                                               log_sum_exp(contamination[m], noise[n,m]) + multinomial_nuisance[multinomPlace])); //estimated abundance if true positive
                }
                multinomPlace += 1;
                Xplace += M[d];
            }
            target += multi_student_t_lupdf(noise |
                                            nu_residuals,
                                            predicted,
                                            cov);
        } else {
            vector[M[d]] phi = inv_square(contaminant_overDisp[d] * segment(var_scales, sumMplus[d] + 1, M[d]));
            for(n in 1:sumID[d]) {
                vector[M[d]] seg = segment(latent_props_raw, Xplace, M[d]);
                for(m in 1:M[d]) {
                    target += log_sum_exp(log1m_inv_logit(logit_prob_present[m,n])
                                          + neg_binomial_2_log_lpmf(X[Xplace + m - 1] |
                                                                    contamination[m] + multinomial_nuisance[multinomPlace],
                                                                    phi[m]), //estimated abundance if true negative
                                            log_inv_logit(logit_prob_present[m,n])
                                            + poisson_log_lpmf(X[Xplace + m - 1] |
                                                               log_sum_exp(contamination[m], seg[m]) + multinomial_nuisance[multinomPlace])); //estimated abundance if true positive
                }
                target += student_t_lupdf(seg |
                                          nu_residuals,
                                          predictTemp[,n],
                                          segment(var_scales, sumMplus[d] + 1, M[d]));
                multinomPlace += 1;
                Xplace += M[d];
            }
        }
    } // count likelihood
    for(r in 1:R) {
        if(Mplus[D+r] > M[D+r]) {
            vector[max(nMat[r,])] noise[N+nVarGroups];
            matrix[M[D+r],M[D+r]] cov
                = add_diag(tcrossprod(diag_post_multiply(
                                   to_matrix(segment(mm, sumMMplus[D+r] + 1, MMplus[D+r]),
                                             M[D+r],
                                             Mplus[D+r] - M[D+r]),
                                   segment(var_scales, sumMplus[D+r] + M[D+r] + 1, Mplus[D+r] - M[D+r]))),
                           square(segment(var_scales, sumMplus[D+r] + 1, M[D+r])) + 1e-10);
            for(n in 1:(N+nVarGroups)) {
                if(nIsolate[r,n] > 0) {
                    int inds[nIsolate[r,n]] = IRInds[r,n,1:sumIR[r,n]][segInds1[r,n,1:nIsolate[r,n]]];
                    target += student_t_lupdf(segment(PFilled, Pplace, sumIR[r,n])[segInds1[r,n,1:nIsolate[r,n]]] |
                                              nu_residuals,
                                              intercepts[(sumM[D+r] + 1):sumM[D+r+1]][inds]
                                              + W[(sumM[D+r] + 1):sumM[D+r+1],][inds,]
                                                * Z_Z_higher[,n],
                                              sqrt(diagonal(cov)[inds]));
                }
                if(nMat[r,matchInds[r,n]] > 0) {
                    int inds[nMat[r,matchInds[r,n]]] = segInds2[r,matchInds[r,n],1:nMat[r,matchInds[r,n]]];
                    noise[n,1:nMat[r,matchInds[r,n]]] = segment(PFilled, Pplace, sumIR[r,n])[inds];
                }
                Pplace += sumIR[r,n];
            }
            for(m in 1:nUniqueR[r]) {
                if(nMat[r,m] > 0) {
                    int inds[nMat[r,m]] = uniqueRInds[r,m,1:sumIRUnique[r,m]][segInds2[r,m,1:nMat[r,m]]];
                    matrix[nMat[r,m],nMatches[r,m]] predictTemp
                        = rep_matrix(intercepts[(sumM[D+r] + 1):(sumM[D+r+1])][inds], nMatches[r,m])
                          + W[(sumM[D+r] + 1):(sumM[D+r] + M[D+r]),][inds,]
                            * Z_Z_higher[,matchIndsInv[r,m,1:nMatches[r,m]]];
                    vector[nMat[r,m]] predicted[nMatches[r,m]];
                    for(n in 1:nMatches[r,m]) {
                        predicted[n] = predictTemp[,n];
                    }
                    target += multi_student_t_lupdf(noise[matchIndsInv[r,m,1:nMatches[r,m]],1:nMat[r,m]] |
                                                    nu_residuals,
                                                    predicted,
                                                    cov[inds,inds]);
                }
            }
        } else {
            for(n in 1:(N+nVarGroups)) {
                if(sumIR[r,n] > 0) {
                    target += student_t_lupdf(segment(PFilled, Pplace, sumIR[r,n]) |
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
                    vector[Mc[cv]] seg
                        = log_softmax(segment(intercepts, V + O + sumMc[cv] + 1, Mc[cv])
                                      + W[(V + O + sumMc[cv] + 1):(V + O + sumMc[cv] + Mc[cv]),]
                                        * Z_Z_higher[,n]);
                    vector[Mc[cv]] segY = segment(Y, Yplace, Mc[cv]);
                    real terms = 0;
                    for(c in 1:Mc[cv]) {
                        if(segY[c] > 0) {
                            if(terms == 0) {
                                terms = log(segY[c]) + seg[c];
                            } else {
                                terms = log_sum_exp(terms, log(segY[c]) + seg[c]);
                            }
                        }
                    }
                    target += terms;
                } else {
                    real seg
                        = intercepts[V + O + sumMc[cv] + 1]
                          + W[V + O + sumMc[cv] + 1,]
                            * Z_Z_higher[,n];
                    if(Y[Yplace] == 1) {
                        target += log_inv_logit(seg);
                    } else if(Y[Yplace] == 0) {
                        target += log1m_inv_logit(seg);
                    } else {
                        target += log_mix(Y[Yplace],
                                          log_inv_logit(seg),
                                          log1m_inv_logit(seg));
                    }
                }
                Yplace += Mc[cv];
            }
        }
    } // categorical/binary likelihood allowing uncertain (real-valued) data
}
//strong inspiration from https://www.cs.helsinki.fi/u/sakaya/tutorial/code/cca.R
