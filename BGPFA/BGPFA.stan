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
    vector[sum(M[1:D])] scales_observed;
    vector<lower=0>[sum(M_all) + sum(M_all[1:D]) + D] prior_scales;
    vector[sum(M) + sum(M[1:D]) + D] prior_intercept_centers;
    int<lower=0> size_mm;                               // size of model matrix for higher level variables
    vector[size_mm] mm;                                 // model matrix for higher level variables
    int<lower=0> N_Pm;                                  // number of P measured as -inf
    int<lower=0> idx_Pm[N_Pm];                          // indices for missing Ps
    int<lower=0> idx_Pnm[H-N_Pm];                       // indices for nonmissing Ps
    vector[choose(M[size(M)],2)] dist_sites;            // distance between sites
    real rho_sites_prior;                               // prior expectation for length scale of gaussian process on sites
    matrix[N,N_var_groups] samp2group;                  // Matrix to calculate means of grouped samples
    int site_smoothness;                                // Matern smoothness parameter for sites
    real nu_residuals;                                  // residual robustness
    vector[D] inv_log_max_contam;                       // prior expectation of contamination rate
    real<lower=0> shape_gnorm;                          // strength of prior pulling contamination toward zero
    real mass_slow;
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
    int sum_M_higher[DRC+1];                           // Cumulative number of higher variables
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
    real rho_Z_scale = 2.0 * K_linear / K_gp;                // More 'dependent' latent variables per group means, in order to encourage orthogonality, length scale must be smaller. Also compensate for shape
    matrix[site_smoothness-1, 2 * site_smoothness] chooseRJ; // precomputed part of Matern covariance
    matrix[site_smoothness, 2 * site_smoothness] ffKJ;       // precomputed part of Matern covariance
    int idx_scales[KG+1];
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
    //
    idx_scales[1] = 1;
    for(g in 1:KG) idx_scales[g+1] = K_linear + K_gp * (g-1) + 1;
    //
}
parameters {
    matrix[K_linear,N] Z_linear_raw;               // PCA axis scores raw
    matrix[KG*K_gp,N] Z_gp_raw;                    // PCA axis scores for gaussian process raw
    vector[K] alpha_Z_raw;                         // controls concentration for each axis
    vector[K] beta_Z_prop_raw;                     // controls degree of skew for each axis
    matrix[VOB_all+V_all+D,K] W_raw;               // PCA variable loadings
    vector[VOB_all+V_all+D] sds_raw;               // variable scales
    vector[K] latent_var_raw;                   // overall scale of each axis. consider constraining to be positive? no change in priors needed except first entry
    vector[DRC+D] dataset_scales_raw;              // overall scale of each dataset
    matrix[DRC+D,K] weight_scales_raw;             // per-dataset-and-axis scales
    vector[VOB+V+D] intercepts;
    vector[F] abundance_observed_vector;           // count dataset latent log abundance
    vector[F_higher] abundance_higher_vector;      // count dataset residual higher effects on log abundance
    vector[F_higher] prevalence_higher_vector;     // count dataset residual higher effects on prevalence
    vector[H_higher] P_higher;                     // continuous dataset residual higher effects
    vector[G_higher] Y_higher_vector;              // categorical dataset residual higher effects
    row_vector[N_multinom] multinomial_nuisance;   // per-sample parameter to convert independent poisson distributions to a multinomial one
    vector<lower=0>[K] rho_sites;                  // length scale for site gaussian process
    vector<lower=0, upper=1>[K] site_prop;         // importance of site covariance compared to nugget effect
    matrix<lower=0>[K_linear,KG] rho_Z;            // length scale for latent axis gaussian process
    vector<upper=0>[D] inv_log_less_contamination; // smaller = less average contamination
    vector[D] contaminant_overdisp_raw;            // dispersion parameter for amount of contamination in true negative count observations
}
transformed parameters {
    vector[DRC+D] dataset_scales_log = dataset_scales_raw * mass_slow;               // overall scale of each dataset is adapted slowly so starting mode with no latent variables is not explored much
    vector[DRC+D] dataset_scales = exp(dataset_scales_log);                          // overall scale of each dataset
    vector[K] alpha_Z_log = alpha_Z_raw * mass_slow;                                 // concentration is slowly adapted because singularities in extremes cause overflow at beginning
    vector<lower=0>[K] alpha_Z = exp(alpha_Z_log);                                   // concentration
    vector<lower=0,upper=1>[K] beta_Z_prop = inv_logit(beta_Z_prop_raw * mass_slow); // skew is slowly adapted to help keep things from getting stuck in mode of less likely reflection
    vector[VOB_all+V_all+D] sds_log = sds_raw * mass_slow;                           // scales tend to blow up so this is a hack to adjust the starting mass matrix but should have no effect on model
    vector<lower=0>[VOB_all+V_all+D] sds = exp(sds_log);                             // scales tend to blow up so this is a hack to adjust the starting mass matrix but should have no effect on model
    vector[D] contaminant_overdisp_log = contaminant_overdisp_raw * mass_slow;
    vector<lower=0>[D] contaminant_overdisp = exp(contaminant_overdisp_log);
    matrix[VOB_all+V_all+D,K] W_norm;                                  // PCA variable loadings
    vector[K] latent_var_log = latent_var_raw;                   // initialize with raw values
    vector[K] latent_scales;                                           // overall scale of each axis
    matrix[DRC+D,K] weight_scales_log = weight_scales_raw * mass_slow; // per-dataset-and-axis scales
    matrix<lower=0>[DRC+D,K] weight_scales = exp(weight_scales_log);   // per-dataset-and-axis scales
    matrix[K,N] Z;
    matrix[K,N_var_groups] Z_higher;
    matrix[VOB,K] W;
    matrix[V,K] W_binary_counts;
    vector<lower=0>[VOB_all+V_all+D] var_scales = sds .* prior_scales;
    vector[VOB_all+V_all+D] var_scales_log = sds_log + log(prior_scales);
    corr_matrix[M[DRC]] corr_sites[K];
    vector[D] log_less_contamination = inv(inv_log_less_contamination);
    latent_var_log[idx_scales] = log_positive_ordered_transform(latent_var_log[idx_scales]);
    latent_var_log[1:K_linear] = log_positive_ordered_transform(latent_var_log[1:K_linear]);
    Z[1:K_linear,]
        = transform_MVN_kumaraswamy(Z_linear_raw,
                                    alpha_Z[1:K_linear],
                                    alpha_Z[1:K_linear] .* beta_Z_prop[1:K_linear]); // first axes are independent
    for(g in 1:KG) {
        matrix[N,N] L = L_cov_exp_quad_ARD(Z[1:K_linear,], rho_Z[,g], 1e-9)';
        int s = K_gp * (g-1) + 1;
        int f = K_gp * g;
        Z[(K_linear + s):(K_linear + f),]
            = transform_MVN_kumaraswamy(Z_gp_raw[s:f,] * L,
                                        alpha_Z[(K_linear + s):(K_linear + f)],
                                        alpha_Z[(K_linear + s):(K_linear + f)] .* beta_Z_prop[(K_linear + s):(K_linear + f)]);
        latent_var_log[(K_linear + s):(K_linear + f)] = log_positive_ordered_transform(latent_var_log[(K_linear + s):(K_linear + f)]);
    } // other axes are dependent on the first axes through gaussian processes
    latent_scales = exp(0.5 * latent_var_log); // sqrt(exp(latent_var_log)) but more stable
    W_norm = diag_post_multiply(W_raw, latent_scales);
    Z_higher = Z * samp2group;
    for(k in 1:K) {
        corr_sites[k]
            = fill_sym(site_prop[k]
                  //   * circular_matern(dist_sites, site_smoothness, inv(rho_sites[k]), ffKJ, chooseRJ),
                       * exp(-dist_sites / rho_sites[k]),
                       M[DRC],
                       1 + 1e-9);
    } // determine covariance of sites
    for(drc in 1:DRC) {
        W[(sum_M[drc] + 1):(sum_M[drc] + M[drc]),]
            = diag_pre_multiply(segment(var_scales, sum_M_all[drc] + 1, M[drc]),
                                W_norm[(sum_M_all[drc] + 1):(sum_M_all[drc] + M[drc]),]); // effects for each base variable
        if(drc == DRC) {
            for(k in 1:K) {
                W[(sum_M[drc] + 1):(sum_M[drc] + M[drc]),k]
                    = cholesky_decompose(corr_sites[k])
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
        if(log_less_contamination[d] > 0) log_less_contamination[d] = 0; // overflow when inverse is near zero makes this infinite instead of0?
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
    vector[VOB+V+D] intercept_scales;
    vector[VOB_all+V_all+D] dsv;
    matrix[VOB_all+V_all+D,K] wsn;
    vector[F_higher] sds_log_abundance_higher;
    vector[F_higher] sds_log_prevalence_higher;
    vector[H] P_predicted;
    vector[H] var_P_log;
    vector[H_higher] var_P_higher;
    matrix[sum(M[(DR+1):DRC]), N_all] Y_higher_matrix = rep_matrix(0,sum(M[(DR+1):DRC]), N_all);
    int i_X = 1;
    int i_X_higher = 1;
    int i_P = 1;
    int i_P_higher = 1;
    int i_Y = 1;
    int i_Y_higher = 1;
    int i_multinom = 1;
    for(drc in 1:DRC) {
        intercept_scales[(sum_M[drc] + 1):(sum_M[drc] + M[drc])] = segment(var_scales, sum_M_all[drc] + 1, M[drc]);
        if(M_all[drc] > M[drc]) {
            intercept_scales[(sum_M[drc] + 1):(sum_M[drc] + M[drc])] =
                sqrt(square(intercept_scales[(sum_M[drc] + 1):(sum_M[drc] + M[drc])])
                     + to_matrix(segment(mm, sum_MxM_all[drc] + 1, MxM_all[drc]), M[drc], M_higher[drc])
                       * square(segment(var_scales, sum_M_all[drc] + M[drc] + 1, M_higher[drc])));
        }
        for(k in 1:K) {
            wsn[(sum_M_all[drc] + 1):(sum_M_all[drc] + M_all[drc]), k] = rep_vector(weight_scales[drc,k], M_all[drc]);
        }
        dsv[(sum_M_all[drc] + 1):(sum_M_all[drc] + M_all[drc])] = rep_vector(dataset_scales_log[drc], M_all[drc]);
        if(drc <= D) {
            intercept_scales[(VOB + sum_M[drc] + 1):(VOB + sum_M[drc] + M[drc])] = segment(var_scales, VOB_all + sum_M_all[drc] + 1, M[drc]);
            intercept_scales[VOB + V + drc] = var_scales[VOB_all + V_all + drc];
            if(M_all[drc] > M[drc]) {
                intercept_scales[(VOB + sum_M[drc] + 1):(VOB + sum_M[drc] + M[drc])] =
                    sqrt(square(intercept_scales[(VOB + sum_M[drc] + 1):(VOB + sum_M[drc] + M[drc])])
                         + to_matrix(segment(mm, sum_MxM_all[drc] + 1, MxM_all[drc]), M[drc], M_higher[drc])
                           * square(segment(var_scales, VOB_all + sum_M_all[drc] + M[drc] + 1, M_higher[drc])));
            }
            for(k in 1:K) {
                wsn[(VOB_all + sum_M_all[drc] + 1):(VOB_all + sum_M_all[drc] + M_all[drc]), k] = rep_vector(weight_scales[DRC+drc,k], M_all[drc]);
                wsn[VOB_all + V_all + drc, k] = weight_scales[DRC+drc,k];
            }
            dsv[(VOB_all + sum_M_all[drc] + 1):(VOB_all + sum_M_all[drc] + M_all[drc])] = rep_vector(dataset_scales_log[DRC+drc], M_all[drc]);
            dsv[VOB_all + V_all + drc] = dataset_scales_log[DRC+drc];
        }
    }
    // end data wrangling
    // priors
    target += std_normal_lupdf(dataset_scales) + sum(dataset_scales_log);                                    // entire datasets may have uniformly biased difference in scale compared to priors
    target += normal_lupdf(sds_log | dsv, 1);                                                                // per-variable sigmas shrink toward dataset scales
    target += normal_lupdf(intercepts | prior_intercept_centers, 25 * intercept_scales);                     //
    target += inv_gamma_lupdf(latent_scales[1] | 5, 5 * global_scale_prior) + 0.5 * latent_var_log[1];       // first axis has minimum scale of global scale prior, shrunk toward that value, with jacobian correction for putting prior on transformed first scale
    target += generalized_normal_lpdf(latent_var_raw[2:] | 0, rep_vector(0.75, K-1), 10);                    // other axes have scale with logistic generalized normal between zero and the previous axis
    target += generalized_std_normal_lpdf(to_vector(weight_scales) | shape_gnorm) + sum(weight_scales_log);  // sparse selection of datasets per axis
    target += generalized_std_normal_lpdf(inv_log_less_contamination ./ inv_log_max_contam | shape_gnorm);   // shrink amount of contamination in 'true zeros' toward zero
    target += normal_lupdf(contaminant_overdisp_log | 0, 0.1);                                               // shrink overdispersion of contaminant counts in 'true zeros' toward zero
    target += generalized_std_normal_lpdf(alpha_Z_log | shape_gnorm);                                        // extreme values of concentration tend to blow up sampling, but want very flat prior within reasonable range
    target += generalized_std_normal_lpdf(beta_Z_prop_raw * mass_slow + 1 | shape_gnorm);                    // Z should significantly skewed to prevent reflections and orient categorical-like variables. Kuramaswamy skews to left by default so this emphasizes same direction
    target += student_t_lupdf(to_vector(W_raw) | 5, 0, to_vector(wsn));                                      // sparsity should help prevent arbitrary rotation and makes axes more interpretable
    target += std_normal_lupdf(to_vector(Z_linear_raw));                                                     // normal density required for raw parameter in order to transform to kuramaswamy
    target += std_normal_lupdf(to_vector(Z_gp_raw));                                                         // gp effects, prior to correlation with cholesky and transformation to kuramaswamy
    target += inv_gamma_lupdf(rho_sites | 5, 5 * rho_sites_prior);                                           // length scale for gaussian process on sites
    target += inv_gamma_lupdf(to_vector(rho_Z) | rho_Z_shape, rho_Z_scale);                                  // length scale for gaussian process on PCA axis scores
    // end priors
    // likelihoods
    // count likelihood
    for(d in 1:D) {
        matrix[M[d],sum_ID[d]] abundance_predicted
            = rep_matrix(intercepts[(sum_M[d] + 1):(sum_M[d] + M[d])], sum_ID[d])
              + W[(sum_M[d] + 1):(sum_M[d] + M[d]),]
                * Z_Z_higher[,idx_ID[d,1:sum_ID[d]]];
        matrix[M[d],sum_ID[d]] prevalence
            = rep_matrix(intercepts[VOB+V+d]
                         + segment(intercepts, VOB + sum_M[d] + 1, M[d]),
                         sum_ID[d])
              + W_binary_counts[(sum_M[d] + 1):(sum_M[d] + M[d]),]
                * Z_Z_higher[,idx_ID[d,1:sum_ID[d]]];
        matrix[M[d],sum_ID[d]] abundance_contam
            = rep_matrix(segment(intercepts, sum_M[d] + 1, M[d])
                         + log_inv_logit(intercepts[VOB+V+d]
                                      + segment(intercepts, VOB + sum_M[d] + 1, M[d]))
                         + log_less_contamination[d],
              sum_ID[d]);
        matrix[M[d],sum_ID[d]] abundance_observed
            = diag_pre_multiply(segment(scales_observed, sum_M[d] + 1, M[d]),
                                to_matrix(segment(abundance_observed_vector, i_X, M[d] * sum_ID[d]),
                                          M[d], sum_ID[d]));
        if(M_all[d] > M[d]) {
            matrix[M[d],M_higher[d]] MM
                = to_matrix(segment(mm, sum_MxM_all[d] + 1, MxM_all[d]),
                            M[d],
                            M_higher[d]);
            matrix[M[d],sum_ID[d]] higher_summed
                = diag_post_multiply(MM, segment(prior_scales, sum_M_all[d] + M[d] + 1, M_higher[d]))
                  * to_matrix(segment(abundance_higher_vector, i_X_higher, M_higher[d] * sum_ID[d]),
                             M_higher[d], sum_ID[d]);
            abundance_predicted += higher_summed;
            sds_log_abundance_higher[i_X_higher:(i_X_higher + sum_ID[d] * M_higher[d] - 1)]
                = to_vector(rep_matrix(segment(sds_log, sum_M_all[d] + M[d] + 1, M_higher[d]), sum_ID[d]));
            abundance_contam += higher_summed;
            prevalence
                += diag_post_multiply(MM, segment(prior_scales, VOB_all + sum_M_all[d] + M[d] + 1, M_higher[d]))
                   * to_matrix(segment(prevalence_higher_vector, i_X_higher, M_higher[d] * sum_ID[d]),
                               M_higher[d], sum_ID[d]);
            sds_log_prevalence_higher[i_X_higher:(i_X_higher + sum_ID[d] * M_higher[d] - 1)]
                = to_vector(rep_matrix(segment(sds_log, VOB_all + sum_M_all[d] + M[d] + 1, M_higher[d]), sum_ID[d]));
            i_X_higher += M_higher[d] * sum_ID[d];
        }
        target += poisson_log_lpmf(segment(X, i_X, M[d] * sum_ID[d]) |
                                   to_vector(abundance_observed + rep_matrix(segment(multinomial_nuisance, i_multinom, sum_ID[d]), M[d])));
        i_X += M[d] * sum_ID[d];
        i_multinom += sum_ID[d];
        for(n in 1:sum_ID[d]) {
            for(m in 1:M[d]) {
                target += log_sum_exp(log1m_inv_logit(prevalence[m,n])
                                      + student_t_log_lpdf(abundance_observed[m,n] |
                                                           nu_residuals,
                                                           abundance_contam[m,n],
                                                           contaminant_overdisp_log[d] + var_scales_log[sum_M_all[d] + m]), //estimated abundance if true negative
                                      log_inv_logit(prevalence[m,n])
                                      + student_t_log_lpdf(abundance_observed[m,n] |
                                                           nu_residuals,
                                                           log_sum_exp(abundance_contam[m,n], abundance_predicted[m,n]),
                                                           var_scales_log[sum_M_all[d] + m])); //estimated abundance if true positive
            }
        }
    }
    target += student_t_log_v0_lpdf(abundance_higher_vector | nu_residuals, sds_log_abundance_higher);
    target += student_t_log_v0_lpdf(prevalence_higher_vector | nu_residuals, sds_log_prevalence_higher);
    // end count likelihood
    // continuous likelihood
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
                var_P_log[i_P:(i_P + sum_IR[r,n] - 1)] = segment(var_scales_log, sum_M[D+r] + 1, M[D+r])[idx_IR[r,n,1:sum_IR[r,n]]];
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
    target += student_t_log_v_lpdf(P[idx_Pnm] |
                                   nu_residuals,
                                   P_predicted[idx_Pnm],
                                   var_P_log[idx_Pnm]);
    target += student_t_lcdf(P[idx_Pm] |
                             nu_residuals,
                             P_predicted[idx_Pm],
                             exp(var_P_log[idx_Pm]));
    // end continuous likelihood
    // categorical/binary likelihood allowing uncertain (real-valued) data
    for(c in 1:C) {
        matrix[M[DR+c], M_higher[DR+c]] MM
            = to_matrix(segment(mm, sum_MxM_all[DR+c] + 1, MxM_all[DR+c]),
                        M[DR+c],
                        M_higher[DR+c]);
        for(n in 1:N_all) {
            Y_higher_matrix[(sum_M[DR+c] - sum_M[DR+1] + 1):(sum_M[DR+c] - sum_M[DR+1] + M[DR+c]),n]
                = MM[, idx_IC_higher[c,n,1:sum_IC_higher[c,n]]]
                  * segment(Y_higher_vector, i_Y_higher, sum_IC_higher[c,n]);
            target += student_t_lupdf(segment(Y_higher_vector, i_Y_higher, sum_IC_higher[c,n]) |
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
                                      + Y_higher_matrix[(sum_Mc[cv] + 1):(sum_Mc[cv] + Mc[cv]),n]);
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
                          + Y_higher_matrix[sum_Mc[cv] + 1, n];
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
    // end categorical/binary likelihood
}
