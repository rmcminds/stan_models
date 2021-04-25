//strong inspiration from https://www.cs.helsinki.fi/u/sakaya/tutorial/code/cca.R
#include functions.stan
data {
    int<lower=0> N;                                     // Number of samples
    int<lower=0> N_var_groups;                          // Number of groups of samples that have a shared single observation
    int<lower=0> N_all;                                 // Number of samples plus number of groups of samples
    int<lower=0> D;                                     // Number of compositional count datasets
    int<lower=0> R;                                     // Number of continuous datasets
    int<lower=0> C;                                     // Number of categorical datasets
    int<lower=0> M[D+R+C];                              // Number of observed variables in each dataset
    int<lower=0> M_higher[D+R+C];                       // Number of higher parameters for each dataset
    int<lower=0> M_all[D+R+C];                          // Number of parameters for each dataset
    int<lower=0> K_linear;                              // Number of residual latent linear dimensions
    int<lower=0> K_gp;                                  // Number of residual latent GP dimensions per kernel
    int<lower=0> KG;                                    // Number of unique GP kernels
    int<lower=0> I[D,N_all];                            // Sample indices for count datasets
    int<lower=0> O;                                     // Total number of observed continuous variables
    int<lower=0> O_higher;                              // Total number of higher effects for observed continuous variables
    int<lower=0> IR[O, N_all];                          // Sample indices for continuous variables
    int<lower=0> IR_higher[O_higher,N_all];             // Sample indices for higher effects of continuous variables
    int<lower=0> B;                                     // Total number of observed categorical variables
    int<lower=0> IC[B,N_all];                           // Sample indices for categorical datasets
    int<lower=0> B_higher;                              // Total number of observed categorical variables
    int<lower=0> IC_higher[B_higher,N_all];
    int<lower=0> C_vars;                                // Total number of observed categorical variables
    int<lower=0> ICv[C_vars,N_all];                     // Sample indices for categorical datasets
    int<lower=0> Mc[C_vars];                            // Number of levels for each categorical variable
    int<lower=0> F;                                     // Total number of compositional count observations
    int<lower=0> F_higher;                              // Total number of higher effects for compositional count observations
    int<lower=0> X[F];                                  // compositional count data
    int<lower=0> G;                                     // Total number of categorical observations
    int<lower=0> G_higher;                              // Total number of higher effects for categorical observations
    vector<lower=0>[G] Y;                               // categorical data
    int<lower=0> H;                                     // Total number of continuous observations
    int<lower=0> H_higher;                              // Total number of higher effects for continuous observations
    vector[H] P;                                        // continuous data
    real<lower=0> global_scale_prior;
    real<lower=0> ortho_scale_prior;                    // prior degree of orthogonalization regularization
    vector<lower=0>[sum(M_all)] prior_scales;
    vector<lower=0>[sum(M)] prior_intercept_scales;
    vector[sum(M)] prior_intercept_centers;
    vector[sum(M[1:D])] binary_count_intercept_centers;
    int<lower=0> size_mm;                               // size of model matrix for higher level variables
    vector[size_mm] mm;                                 // model matrix for higher level variables
    int<lower=0> N_Pm;                                  // number of P measured as -inf
    int<lower=0> idx_Pm[N_Pm];                          // indices for missing Ps
    real P_max[N_Pm];                                   // lowest measured value for the variable corresponding to a missing P
    real rate_gamma_fact;                               // scale of inverse gamma prior on nu for W
    real shape_gamma_fact;                              // shape of inverse gamma prior on nu for W
    vector[choose(M[size(M)],2)] dist_sites;            // distance between sites
    real rho_sites_prior;                               // prior expectation for length scale of gaussian process on sites
    matrix[N,N_var_groups] samp2group;                  // Matrix to calculate means of grouped samples
    int site_smoothness;                                // Matern smoothness parameter for sites
    real nu_residuals;                                  // residual robustness
    vector[D] inv_log_max_contam;                       // prior expectation of contamination rate
    real<lower=0> shape_gnorm;                          // strength of prior pulling contamination toward zero
    real<lower=0> skew_Z_prior;
}
transformed data {
    int K = K_linear + KG * K_gp;                      // Total number of latent axes
    int DR = D+R;                                      // Total number of count and continuous datasets
    int DRC = DR+C;                                    // Total number of datasets
    int V = 0;                                         // Total number of observed compositional count variables
    int V_all = 0;                                     // Total number of observed compositional count parameters
    int VOB = sum(M);                                  // Total number of observed variables
    int VOB_all = sum(M_all);                          // Total number of parameters per dim
    int sum_M[DRC+1];                                  // Cumulative number of base level variables
    int sum_M_higher[DRC+1];                           // Cumulative number of base level variables
    int sum_M_all[DRC+1];                              // Cumulative number of base level and higher variables
    int MxM_all[DRC];                                  // size of model matrix for higher level variables
    int sum_MxM_all[DRC+1];                            // Cumulative size of model matrices for higher level variables
    int sum_Mc[C_vars+1];                              // Cumulative number of categorical levels
    int N_multinom = 0;                                // Total number of multinomial samples
    int sum_ID[D];                                     // Number of samples in each count dataset
    int idx_ID[D,N_all];                               // Indices mapping each sample to each count dataset
    int sum_IR[R,N_all];                               // Number of continuous observations in each sample, per dataset
    int idx_IR[R,N_all,max(M[(D+1):(DR)])];            // Indices mapping each sample to each continuous dataset
    int sum_IR_higher[R,N_all];                        // Number of continuous higher effects in each sample, per dataset
    int idx_IR_higher[R,N_all,max(M_higher[(D+1):(DR)])];    // Indices mapping each sample to each continuous dataset
    int sum_IC[C,N_all];                                     // Number of categorical observations in each sample, per dataset
    int idx_IC[C,N_all,max(M[(DR+1):(DRC)])];                // Indices mapping each sample to each categorical dataset
    int sum_IC_higher[C,N_all];                              // Number of categorical higher effects in each sample, per dataset
    int idx_IC_higher[C,N_all,max(M_higher[(DR+1):(DRC)])];  // Indices mapping each sample to each categorical dataset
    real rho_Z_shape = 1.0 / K_linear;                       // More 'independent' latent variables means more need for sparsity in feature selection
    real rho_Z_scale = 6.0 / K_gp;                           // More 'dependent' latent variables per group means, in order to encourage orthogonality, length scale must be smaller
    matrix[site_smoothness-1, 2 * site_smoothness] chooseRJ; // precomputed part of Matern covariance
    matrix[site_smoothness, 2 * site_smoothness] ffKJ;       // precomputed part of Matern covariance
    // create indices for all datasets
    sum_M[1] = 0;
    sum_M_higher[1] = 0;
    sum_M_all[1] = 0;
    sum_MxM_all[1] = 0;
    for(drc in 1:DRC) {
        sum_M[drc+1] = sum(M[1:drc]);
        sum_M_higher[drc+1] = sum(M_higher[1:drc]);
        sum_M_all[drc+1] = sum(M_all[1:drc]);
        MxM_all[drc] = M[drc] * M_higher[drc];
        sum_MxM_all[drc+1] = sum(MxM_all[1:drc]);
    }
    //
    // create indices for count datasets
    V = sum_M[D+1];
    V_all = sum_M_all[D+1];
    for(d in 1:D) {
        int nplace = 1;
        sum_ID[d] = sum(I[d,]);
        for(n in 1:N_all) {
            if(I[d,n]) {
                idx_ID[d,nplace] = n;
                nplace += 1;
            }
        }
    }
    N_multinom = sum(sum_ID);
    //
    // create indices for continuous variables
    for(r in 1:R) {
        for(n in 1:N_all) {
            int segplace = 1;
            int segplace_higher = 1;
            sum_IR[r,n] = sum(IR[(sum_M[D+r]-V+1):(sum_M[D+r+1]-V),n]);
            sum_IR_higher[r,n] = sum(IR_higher[(sum_M_higher[D+r]-sum_M_higher[D+1]+1):(sum_M_higher[D+r+1]-sum_M_higher[D+1]),n]);
            for(o in 1:M[D+r]) {
                if(IR[sum_M[D+r]-V+o, n]) {
                    idx_IR[r,n,segplace] = o;
                    segplace += 1;
                }
            }
            for(o in 1:M_higher[D+r]) {
                if(IR_higher[sum_M_higher[D+r]-sum_M_higher[D+1]+o, n]) {
                    idx_IR_higher[r,n,segplace_higher] = o;
                    segplace_higher += 1;
                }
            }
        }
    }
    // end creating continuous variable indices
    // create indices for categorical variables
    sum_Mc[1] = 0;
    for(cv in 1:C_vars) {
        sum_Mc[cv+1] = sum(Mc[1:cv]);
    }
    for(c in 1:C) {
        for(n in 1:N_all) {
            int segplace = 1;
            int segplace_higher = 1;
            sum_IC[c,n] = sum(IC[(sum_M[DR+c]-V-O+1):(sum_M[DR+1]-V-O),n]);
            sum_IC_higher[c,n] = sum(IC_higher[(sum_M_higher[DR+c]-sum_M_higher[DR+1]+1):(sum_M_higher[DR+c+1]-sum_M_higher[DR+1]),n]);
            for(b in 1:M[DR+c]) {
                if(IC[sum_M[DR+c]-V-O+b, n]) {
                    idx_IC[c,n,segplace] = b;
                    segplace += 1;
                }
            }
            for(b in 1:M_higher[DR+c]) {
                if(IC_higher[sum_M_higher[DR+c]-sum_M_higher[DR+1]+b, n]) {
                    idx_IC_higher[c,n,segplace_higher] = b;
                    segplace_higher += 1;
                }
            }
        }
    }
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
    matrix[K,N] Z1;             // PCA axis scores, normal part
    matrix<lower=0>[K,N] Z2;    // PCA axis scores, skew part
    vector<lower=0>[K] skew_Z;
    matrix[VOB_all+V_all+D,K] W_norm;              // PCA variable loadings
    vector<lower=0>[VOB_all+V_all+D] sds;          // variable scales
    real<lower=0> global_effect_scale;             // overall scale of variable loadings
    real<lower=0> ortho_scale;                     // inverse strength of orthogonality shrinkage
    real<lower=0,upper=1> order_prior_scales;      //
    row_vector<lower=0>[K] latent_scales;      // overall scale of each axis
    vector<lower=0>[DRC+D] dataset_scales;         // overall scale of each dataset
    matrix<lower=0>[DRC+D,K] weight_scales;        // per-dataset-and-axis scales
    vector[VOB] intercepts;
    vector[V] binary_count_intercepts;
    vector[F] abundance_true_vector;               // count dataset latent log abundance
    vector[F_higher] abundance_higher_vector;
    vector[F_higher] prevalence_higher_vector;
    vector[H_higher] P_higher;
    vector[G_higher] Y_higher;
    vector[D] binary_count_dataset_intercepts;     // each count dataset could have biased prior intercepts
    vector[N_multinom] multinomial_nuisance;       // per-sample parameter to convert independent poisson distributions to a multinomial one
    vector<upper=0>[N_Pm] P_missing;               // latent estimates of continuous variables with partial information (truncated observations)
    matrix<lower=0>[DRC+D,K] nu_factors_raw;       // per-dataset-and-axis sparsity of variable loadings
    vector<lower=0>[K] rho_sites;                  // length scale for site gaussian process
    vector<lower=0, upper=1>[K] site_prop;         // importance of site covariance compared to nugget effect
    matrix<lower=0>[K_linear,KG] rho_Z;            // length scale for latent axis gaussian process
    vector<upper=0>[D] inv_log_less_contamination; // smaller = less average contamination
    vector<lower=0>[D] contaminant_overdisp;       // dispersion parameter for amount of contamination in true negative count observations
}
transformed parameters {
    vector[K] delta = skew_Z / sqrt(1 - skew_Z^2);
    matrix[K,N] Z
        = diag_pre_multiply(inv_sqrt(1 - 2 * square(delta) / pi()),
                            diag_pre_multiply(delta ./ skew_Z,
                                              Z1 + diag_pre_multiply(skew_Z, Z2))
                            - rep_matrix(delta * sqrt(2 / pi()), N));             // skew-normal PCA axis scores with mean 0 and sd 1, assuming input Z1 and Z2 have those location and scales
    matrix[K,N_var_groups] Z_higher = Z * samp2group;
    matrix[K_linear + KG * K_gp, N] Z_ortho = diag_pre_multiply(sqrt(rows_dot_self(Z)), svd_V(Z')) * svd_U(Z')';
    matrix[VOB_all+V_all+D,K] W_ortho = svd_U(W_norm) * diag_post_multiply(svd_V(W_norm)', sqrt(columns_dot_self(W_norm)));
    matrix[VOB,K] W;
    matrix[V,K] W_binary_counts;
    vector<lower=0>[VOB_all+V_all+D] var_scales
        = sds
          .* append_row(prior_scales,
                        ones_vector(V_all+D));
    vector<upper=0>[D] log_less_contamination = inv(inv_log_less_contamination);
    vector[H] P_filled = P;
    cov_matrix[M[DRC]] cov_sites[K];
    matrix<lower=2>[DRC,K] nu_factors = nu_factors_raw + 2;
    for(miss in 1:N_Pm) P_filled[idx_Pm[miss]] = P_missing[miss] + P_max[miss];
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
        W[(sum_M[drc] + 1):(sum_M[drc] + M[drc]),]
            = diag_pre_multiply(segment(var_scales, sum_M_all[drc] + 1, M[drc]),
                                W_norm[(sum_M_all[drc] + 1):(sum_M_all[drc] + M[drc]),]); // effects for each base variable
        if(drc == DRC) {
            for(k in 1:K) {
                W[(sum_M[drc] + 1):(sum_M[drc] + M[drc]),k]
                    = cholesky_decompose(cov_sites[k])
                      * sqrt(nu_factors_raw[DRC,k] / nu_factors[DRC,k])
                      * W[(sum_M[drc] + 1):(sum_M[drc] + M[drc]),k];
            }
        } // induce correlation among sites
        if(M_all[drc] > M[drc]) {
            W[(sum_M[drc] + 1):(sum_M[drc] + M[drc]),]
                += to_matrix(segment(mm, sum_MxM_all[drc] + 1, MxM_all[drc]), M[drc], M_higher[drc])
                   * diag_pre_multiply(segment(var_scales, sum_M_all[drc] + M[drc] + 1, M_higher[drc]),
                                       W_norm[(sum_M_all[drc] + M[drc] + 1):(sum_M_all[drc] + M_all[drc]),]); // add higher level effects
        } // add higher level effects
    } // scale PCA factor loadings and combine linear effects
    for(d in 1:D) {
        W_binary_counts[(sum_M[d] + 1):(sum_M[d] + M[d]),]
            = diag_pre_multiply(segment(var_scales, VOB_all + sum_M_all[d] + 1, M[d]),
                                W_norm[(VOB_all + sum_M_all[d] + 1):(VOB_all + sum_M_all[d] + M[d]),]);
        if(M_all[d] > M[d]) {
            W_binary_counts[(sum_M[d] + 1):(sum_M[d] + M[d]),]
                += to_matrix(segment(mm, sum_MxM_all[d] + 1, MxM_all[d]), M[d], M_higher[d])
                   * diag_pre_multiply(segment(var_scales, VOB_all + sum_M_all[d] + M[d] + 1, M_higher[d]),
                                       W_norm[(VOB_all + sum_M_all[d] + M[d] + 1):(VOB_all + sum_M_all[d] + M_all[d]),]);
        }
        W_binary_counts[(sum_M[d] + 1):(sum_M[d] + M[d]),]
            += rep_matrix(var_scales[VOB_all+V_all+d] * W_norm[VOB_all+V_all+d,], M[d]);
    } // same as above, for prevalence estimates in count data
}
model {
    // data wrangling (should replace some with transformed data indices)
    matrix[K,N_all] Z_Z_higher = append_col(Z,Z_higher);
    vector[VOB_all+V_all+D] dsv;
    matrix[VOB_all+V_all+D,K] num;
    matrix[VOB_all+V_all+D,K] wsn;
    vector[H] P_predicted;
    vector[H] var_P;
    vector[H_higher] var_P_higher;
    int i_X = 1;
    int i_X_higher = 1;
    int i_P = 1;
    int i_P_higher = 1;
    int i_Y = 1;
    int i_Y_higher = 1;
    int i_multinom = 1;
    for(drc in 1:DRC) {
        dsv[(sum_M_all[drc] + 1):(sum_M_all[drc] + M_all[drc])] = rep_vector(dataset_scales[drc], M_all[drc]);
        for(k in 1:K) {
            num[(sum_M_all[drc] + 1):(sum_M_all[drc] + M_all[drc]),k]
                = rep_vector(nu_factors[drc,k],
                             M_all[drc]);
            wsn[(sum_M_all[drc] + 1):(sum_M_all[drc] + M_all[drc]),k]
                = rep_vector(weight_scales[drc,k]
                             * sqrt(nu_factors_raw[drc,k] / nu_factors[drc,k]),
                             M_all[drc]);
        }
        if(drc <= D) {
            dsv[(VOB_all + sum_M_all[drc] + 1):(VOB_all + sum_M_all[drc] + M_all[drc])] = rep_vector(dataset_scales[DRC+drc], M_all[drc]);
            dsv[VOB_all + V_all + drc] = dataset_scales[drc];
            for(k in 1:K) {
                num[(VOB_all + sum_M_all[drc] + 1):(VOB_all + sum_M_all[drc] + M_all[drc]),k]
                    = rep_vector(nu_factors[drc,k],
                                 M_all[drc]);
                num[VOB_all+V_all+drc,k]
                    = nu_factors[DRC+drc,k];
                wsn[(VOB_all + sum_M_all[drc] + 1):(VOB_all + sum_M_all[drc] + M_all[drc]),k]
                    = rep_vector(weight_scales[drc,k]
                                 * sqrt(nu_factors_raw[drc,k] / nu_factors[drc,k]),
                                 M_all[drc]);
                wsn[VOB_all+V_all+drc,k]
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
    target += student_t_lupdf(latent_scales[K_linear] | 2, 0, global_effect_scale * order_prior_scales^(0.5*K_linear));// final axis scale centered on global scale diminished by the distance between scales K/2 times
    target += student_t_lupdf(latent_scales[1:(K_linear-1)] | 2, 0, latent_scales[2:K_linear] / order_prior_scales);   // each axis scale has prior expectation to be larger than the next by a fit distance
    target += student_t_lupdf(to_vector(weight_scales) | 2, 0, to_vector(rep_matrix(latent_scales, DRC+D))); // sparse selection of datasets per axis
    target += generalized_normal_lpdf(inv_log_less_contamination | 0, inv_log_max_contam, shape_gnorm);      // shrink amount of contamination in 'true zeros' toward zero
    target += std_normal_lupdf(contaminant_overdisp);                                                        // shrink overdispersion of contaminant counts in 'true zeros' toward zero
    target += normal_lupdf(to_vector(W_norm) | to_vector(W_ortho), global_effect_scale * ortho_scale);       // shrink PCA variable loadings toward closes orthogonal matrix
    target += normal_lupdf(to_vector(Z) | to_vector(Z_ortho), ortho_scale);                                  // shrink PCA axis scores toward closes orthogonal matrix
    target += std_normal_lupdf(to_vector(Z1[1:K_linear,]));                                                  // first PCA axis scores are independent of one another
    target += std_normal_lupdf(to_vector(Z2));                                                                // PCA scores have idependent positive skew to help identify
    target += inv_gamma_lupdf(skew_Z | 5, 5 * skew_Z_prior);                                                // all Z are positively skewed, but each axis varies
    target += inv_gamma_lupdf(rho_sites | 5, 5 * rho_sites_prior);                                           // length scale for gaussian process on sites
    target += inv_gamma_lupdf(to_vector(rho_Z) | rho_Z_shape, rho_Z_scale);                                  // length scale for gaussian process on PCA axis scores
    for(g in 1:KG) {
        target += student_t_lupdf(latent_scales[K_linear + K_gp * g] | 2, 0, global_effect_scale * order_prior_scales^(0.5*K_gp));
        target += student_t_lupdf(latent_scales[(K_linear + (K_gp * (g-1)) + 1):(K_linear + K_gp * g - 1)] |
                                  2,
                                  0,
                                  latent_scales[(K_linear + (K_gp * (g-1)) + 1):(K_linear + K_gp * g - 1)] / order_prior_scales);
        matrix[N_all,N_all] L = L_cov_exp_quad_ARD(Z[1:K_linear,], rho_Z[,g], 1e-9);
        target += multi_gp_cholesky_lupdf(Z1[(K_linear + (K_gp * (g-1)) + 1):(K_linear + K_gp * g),] |
                                          L,
                                          ones_vector(K_gp));
        target += multi_gp_cholesky_lupdf(Z2[(K_linear + (K_gp * (g-1)) + 1):(K_linear + K_gp * g),] |
                                          L,
                                          ones_vector(K_gp));
    }                                                                                                        // final KG * K_gp PCA axis scores are functions of first K_linear ones
    target += normal_lupdf(sds | 0, dsv);                                                                    // per-variable sigmas shrink toward dataset scales
    target += student_t_lupdf(to_vector(W_norm) | to_vector(num), 0, to_vector(wsn));                        // PCA loadings shrink to zero with axis-and-dataset-specific nu and variance
    // end priors
    // likelihoods
    for(d in 1:D) {
        matrix[M[d],sum_ID[d]] predicted
            = rep_matrix(intercepts[(sum_M[d] + 1):(sum_M[d] + M[d])], sum_ID[d])
              + W[(sum_M[d] + 1):(sum_M[d] + M[d]),]
                * Z_Z_higher[,idx_ID[d,1:sum_ID[d]]];
        matrix[M[d],sum_ID[d]] prevalence
            = rep_matrix(binary_count_dataset_intercepts[d]
                         + segment(binary_count_intercepts, sum_M[d] + 1, M[d]),
                         sum_ID[d])
              + W_binary_counts[(sum_M[d] + 1):(sum_M[d] + M[d]),]
                * Z_Z_higher[,idx_ID[d,1:sum_ID[d]]];
        vector[M[d]] abundance_contam
            = segment(intercepts, sum_M[d] + 1, M[d])
              + log_inv_logit(binary_count_dataset_intercepts[d]
                              + segment(binary_count_intercepts, sum_M[d] + 1, M[d]))
              + log_less_contamination[d];
        matrix[M[d],sum_ID[d]] abundance_true
            = to_matrix(segment(abundance_true_vector, i_X, M[d] * sum_ID[d]),
                        M[d], sum_ID[d]);
        vector[M[d]] phi;
        if(M_all[d] > M[d]) {
            matrix[M[d],M_higher[d]] MM
                = to_matrix(segment(mm, sum_MxM_all[d] + 1, MxM_all[d]),
                            M[d],
                            M_higher[d]);
            phi = inv_square(contaminant_overdisp[d])
                  * inv(rows_dot_self(diag_post_multiply(MM,
                                                         segment(var_scales, sum_M_all[d] + M[d] + 1, M_higher[d])))
                        + square(segment(var_scales, sum_M_all[d] + 1, M[d])));
            predicted
                += MM * to_matrix(segment(abundance_higher_vector, i_X_higher, M_higher[d] * sum_ID[d]),
                                  M_higher[d], sum_ID[d]);
            prevalence
                += MM * to_matrix(segment(prevalence_higher_vector, i_X_higher, M_higher[d] * sum_ID[d]),
                                  M_higher[d], sum_ID[d]);
            target += student_t_lpdf(segment(abundance_higher_vector, i_X_higher, M_higher[d] * sum_ID[d]) |
                                     nu_residuals,
                                     0,
                                     to_vector(rep_matrix(segment(var_scales, sum_M_all[d] + M[d] + 1, M_higher[d]), sum_ID[d])));
            target += student_t_lpdf(segment(prevalence_higher_vector, i_X_higher, M_higher[d] * sum_ID[d]) |
                                     nu_residuals,
                                     0,
                                     to_vector(rep_matrix(segment(var_scales, VOB_all + sum_M_all[d] + M[d] + 1, M_higher[d]), sum_ID[d])));
            i_X_higher += M_higher[d] * sum_ID[d];
        } else {
            phi = inv_square(contaminant_overdisp[d] * segment(var_scales, sum_M_all[d] + 1, M[d]));
        }
        target += student_t_lpdf(to_vector(abundance_true) |
                                 nu_residuals,
                                 to_vector(predicted),
                                 to_vector(rep_matrix(segment(var_scales, sum_M_all[d] + 1, M[d]), sum_ID[d])));
        for(n in 1:sum_ID[d]) {
            for(m in 1:M[d]) {
                target += log_sum_exp(log1m_inv_logit(prevalence[m,n])
                                      + neg_binomial_2_log_lpmf(X[i_X + m - 1] |
                                                                abundance_contam[m] + multinomial_nuisance[i_multinom],
                                                                phi[m]), //estimated abundance if true negative
                                        log_inv_logit(prevalence[m,n])
                                        + poisson_log_lpmf(X[i_X + m - 1] |
                                                           log_sum_exp(abundance_contam[m], abundance_true[m,n])
                                                           + multinomial_nuisance[i_multinom])); //estimated abundance if true positive
            }
            i_multinom += 1;
            i_X += M[d];
        }
    } // count likelihood
    for(r in 1:R) {
        matrix[M[D+r],M_higher[D+r]] MM;
        if(M_all[D+r] > M[D+r]) {
            MM = to_matrix(segment(mm, sum_MxM_all[D+r] + 1, MxM_all[D+r]),
                           M[D+r],
                           M_higher[D+r]);
        }
        for(n in 1:N_all) {
            if(sum_IR[r,n] > 0) {
                P_predicted[i_P:(i_P + sum_IR[r,n] - 1)]
                    = intercepts[(sum_M[D+r] + 1):sum_M[D+r+1]][idx_IR[r,n,1:sum_IR[r,n]]]
                      + W[(sum_M[D+r] + 1):sum_M[D+r+1],][idx_IR[r,n,1:sum_IR[r,n]],]
                        * Z_Z_higher[,n];
                var_P[i_P:(i_P + sum_IR[r,n] - 1)] = segment(var_scales, sum_M[D+r] + 1, M[D+r])[idx_IR[r,n,1:sum_IR[r,n]]];
                if(sum_IR_higher[r,n] > 0) {
                        P_predicted[i_P:(i_P + sum_IR[r,n] - 1)]
                            += MM[idx_IR[r,n,1:sum_IR[r,n]], idx_IR_higher[r,n,1:sum_IR_higher[r,n]]]
                               * segment(P_higher, i_P_higher, sum_IR_higher[r,n]);
                        var_P_higher[i_P_higher:(i_P_higher + sum_IR_higher[r,n] - 1)]
                            = segment(var_scales, sum_M_all[D+r] + M[D+r] + 1, M_higher[D+r])[idx_IR_higher[r,n,1:sum_IR_higher[r,n]]];
                        i_P_higher += sum_IR_higher[r,n];
                }
                i_P += sum_IR[r,n];
            }
        }
    }
    target += student_t_lupdf(P_higher |
                              nu_residuals,
                              0,
                              var_P_higher);
    target += student_t_lupdf(P_filled |
                              nu_residuals,
                              P_predicted,
                              var_P);  // continuous likelihood
    {
        matrix[sum(M[(DR+1):DRC]), N_all] resids = rep_matrix(0,sum(M[(DR+1):DRC]), N_all);
        for(c in 1:C) {
            matrix[M[DR+c], M_higher[DR+c]] MM
                = to_matrix(segment(mm, sum_MxM_all[DR+c] + 1, MxM_all[DR+c]),
                            M[DR+c],
                            M_higher[DR+c]);
            for(n in 1:N_all) {
                resids[(sum_M[DR+c] - sum_M[DR+1] + 1):(sum_M[DR+c] - sum_M[DR+1] + M[DR+c]),n]
                    = MM[, idx_IC_higher[c,n,1:sum_IC_higher[c,n]]]
                      * segment(Y_higher, i_Y_higher, sum_IC_higher[c,n]);
                target += student_t_lupdf(segment(Y_higher, i_Y_higher, sum_IC_higher[c,n]) |
                                          nu_residuals,
                                          0,
                                          segment(var_scales, sum_M_all[DR+c] + M[DR+c] + 1, M_higher[DR+c])[idx_IC_higher[c,n,1:sum_IC_higher[c,n]]]);
                i_Y_higher += sum_IC_higher[c,n];
            }
        }
        for(cv in 1:C_vars) {
            for(n in 1:N_all) {
                if(ICv[cv,n]) {
                    if(Mc[cv] > 1) {
                        vector[Mc[cv]] predicted
                            = log_softmax(segment(intercepts, V + O + sum_Mc[cv] + 1, Mc[cv])
                                          + W[(V + O + sum_Mc[cv] + 1):(V + O + sum_Mc[cv] + Mc[cv]),]
                                            * Z_Z_higher[,n]
                                          + resids[(sum_Mc[cv] + 1):(sum_Mc[cv] + Mc[cv]),n]);
                        vector[Mc[cv]] observed = segment(Y, i_Y, Mc[cv]);
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
                            = intercepts[V + O + sum_Mc[cv] + 1]
                              + W[V + O + sum_Mc[cv] + 1,]
                                * Z_Z_higher[,n]
                              + resids[sum_Mc[cv] + 1, n];
                        if(Y[i_Y] == 1) {
                            target += log_inv_logit(predicted);
                        } else if(Y[i_Y] == 0) {
                            target += log1m_inv_logit(predicted);
                        } else {
                            target += log_mix(Y[i_Y],
                                              log_inv_logit(predicted),
                                              log1m_inv_logit(predicted));
                        }
                    } // likelihood for binary variables
                    i_Y += Mc[cv];
                }
            }
        }
    } // categorical/binary likelihood allowing uncertain (real-valued) data
}
