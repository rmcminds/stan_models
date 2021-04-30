functions {
    real generalized_normal_lpdf(vector y, real mu, vector alpha, real beta) {
        return sum(log(beta) - log(2) - log(alpha) - lgamma(inv(beta)) - exp(beta * log(fabs(y-mu)./alpha)));
    }
    real multi_student_t_cholesky_lpdf(matrix y, real nu, matrix mu, matrix L) {
        int N = cols(y);
        int K = cols(L);
        real lp
            = N * (- 0.5 * K * log(nu)
                   + lgamma(0.5 * (nu + K))
                   - lgamma(0.5 * nu)
                   - sum(log(diagonal(L))));
        return(lp - 0.5 * (nu + K) * sum(log1p(columns_dot_self(L * (y - mu)) / nu)));
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
    vector[] to_vector_array(matrix x) {
        vector[rows(x)] y[cols(x)];
        for(j in 1:cols(x)) y[j] = x[,j];
        return(y);
    }
    matrix cholesky_like(matrix A) {
        vector[rows(A)] lambda = eigenvalues_sym(A);
        for(i in 1:rows(A)) if(lambda[i] < 0) lambda[i] = 0;
        return(qr_thin_R(diag_post_multiply(eigenvectors_sym(A), sqrt(lambda))')');
    } // https://math.stackexchange.com/questions/423138/cholesky-for-non-positive-definite-matrices
    matrix mix_skew_normal(matrix Z1, matrix Z2, vector alpha) {
        vector[rows(Z1)] delta = alpha ./ sqrt(1 + alpha^2);
        matrix[rows(Z1),cols(Z1)] Z
            = diag_pre_multiply(inv_sqrt(1 - square(delta) * 2 / pi()),
                                diag_pre_multiply(delta ./ alpha,
                                                  Z1 + diag_pre_multiply(alpha, Z2))
                                - rep_matrix(delta * sqrt(2 / pi()), cols(Z1)));
        return(Z);
    } // skew-normal matrix with mean 0 and sd 1, assuming Z1 is normal, Z2 is half-normal, and each have location 0 and scale 1
    vector transform_tMVN_vector_lp(matrix L, vector u) {
        int N = rows(u);
        vector[N] z;
            for(n in 1:N) {
                int nm1 = n - 1;
                real u_star = Phi((n > 1) ? L[n,1:nm1] * head(z,nm1) / L[n,n] : 0);
                target += log1m(u_star);
                z[n] = inv_Phi(u_star + u[n] - u_star * u[n]);
            }
            return z;
    } // simplified from https://discourse.mc-stan.org/t/multi-normal-lcdf/6652/6
    matrix transform_tMVN_lp(matrix u, matrix L) {
        int K = rows(u);
        int N = cols(u);
        matrix[K,N] z;
            for(n in 1:N) {
                int nm1 = n - 1;
                vector[K] u_star;
                if(n > 1) {
                    u_star = Phi(z[,1:nm1] * L[1:nm1,n] / L[n,n]);
                } else {
                    u_star = rep_vector(0.5,K);
                }
                z[,n] = inv_Phi(u_star + u[,n] - u_star .* u[,n]);
                target += log1m(u_star);
            }
        return z * L;
    } // vectorized and transposed from above
    matrix nearest_special(matrix x) {
        matrix[rows(x),cols(x)] U = svd_U(x);
        matrix[cols(x),cols(x)] V_prime = svd_V(x)';
        matrix[rows(x),cols(x)] O = U * V_prime;
        matrix[rows(x),cols(x)] Q = qr_thin_Q(O);
        if(prod(O[1,].* Q[1,]) < 0) {
            U[,cols(U)] = -U[,cols(U)];
            return(U * diag_post_multiply(V_prime, sqrt(columns_dot_self(x))));
        } else {
            return(diag_post_multiply(O, sqrt(columns_dot_self(x))));
        }
    } // https://math.stackexchange.com/a/3481936
    matrix nearest_ps(matrix x) {
        matrix[rows(x),cols(x)] O = svd_U(x) * svd_V(x)';
        matrix[cols(x),cols(x)] Q = qr_thin_Q(O);
        for(i in 1:cols(x)) {
            if((O[1,i] * Q[1,i]) < 0) {
                O[,i:] = rep_matrix(0, rows(O), cols(O)-i+1);
                break;
            }
        }
        return(diag_post_multiply(O, sqrt(columns_dot_self(x))));
    } // https://scicomp.stackexchange.com/questions/30631/how-to-find-the-nearest-a-near-positive-definite-from-a-given-matrix
    matrix mean_orthogonal_points(matrix x) {
        matrix[rows(x),cols(x)] O = svd_U(x) * svd_V(x)';
        matrix[rows(x),cols(x)] weighted_mean;
        for(i in 1:cols(x)) {
            vector[cols(x)] dists;
            for(j in 1:cols(x)) {
                dists[j] = squared_distance(x[,i], O[,j]);
            }
            weighted_mean[,i] = diag_post_multiply(O, inv(dists)) * rep_vector(inv(sum(inv(dists))), cols(x));
        }
        return(diag_post_multiply(weighted_mean, sqrt(columns_dot_self(x))));
    } // https://scicomp.stackexchange.com/questions/30631/how-to-find-the-nearest-a-near-positive-definite-from-a-given-matrix
    matrix mean_special_orthogonal_points(matrix x) {
        matrix[rows(x),cols(x)] Q = qr_thin_Q(svd_U(x) * svd_V(x)');
        matrix[rows(x),cols(x)] weighted_mean;
        for(i in 1:cols(x)) {
            vector[cols(x)] dists;
            for(j in 1:cols(x)) {
                dists[j] = squared_distance(x[,i], Q[,j]);
            }
            weighted_mean[,i] = diag_post_multiply(Q, inv(dists)) * rep_vector(inv(sum(inv(dists))), cols(x));
        }
        return(diag_post_multiply(weighted_mean, sqrt(columns_dot_self(x))));
    } // https://scicomp.stackexchange.com/questions/30631/how-to-find-the-nearest-a-near-positive-definite-from-a-given-matrix
}

functions {
    real generalized_normal_lpdf(vector y, real mu, vector alpha, real beta) {
        return sum(log(beta) - log(2) - log(alpha) - lgamma(inv(beta)) - exp(beta * log(fabs(y-mu)./alpha)));
    }
    real multi_student_t_cholesky_lpdf(matrix y, real nu, matrix mu, matrix L) {
        int N = cols(y);
        int K = cols(L);
        real lp
            = N * (- 0.5 * K * log(nu)
                   + lgamma(0.5 * (nu + K))
                   - lgamma(0.5 * nu)
                   - sum(log(diagonal(L))));
        return(lp - 0.5 * (nu + K) * sum(log1p(columns_dot_self(L * (y - mu)) / nu)));
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
    vector[] to_vector_array(matrix x) {
        vector[rows(x)] y[cols(x)];
        for(j in 1:cols(x)) y[j] = x[,j];
        return(y);
    }
    matrix cholesky_like(matrix A) {
        vector[rows(A)] lambda = eigenvalues_sym(A);
        for(i in 1:rows(A)) if(lambda[i] < 0) lambda[i] = 0;
        return(qr_thin_R(diag_post_multiply(eigenvectors_sym(A), sqrt(lambda))')');
    } // https://math.stackexchange.com/questions/423138/cholesky-for-non-positive-definite-matrices
    matrix mix_skew_normal(matrix Z1, matrix Z2, vector alpha) {
        vector[rows(Z1)] delta = alpha ./ sqrt(1 + alpha^2);
        matrix[rows(Z1),cols(Z1)] Z
            = diag_pre_multiply(inv_sqrt(1 - square(delta) * 2 / pi()),
                                diag_pre_multiply(delta ./ alpha,
                                                  Z1 + diag_pre_multiply(alpha, Z2))
                                - rep_matrix(delta * sqrt(2 / pi()), cols(Z1)));
        return(Z);
    } // skew-normal matrix with mean 0 and sd 1, assuming Z1 is normal, Z2 is half-normal, and each have location 0 and scale 1
    vector transform_tMVN_vector_lp(matrix L, vector u) {
        int N = rows(u);
        vector[N] z;
            for(n in 1:N) {
                int nm1 = n - 1;
                real u_star = Phi((n > 1) ? L[n,1:nm1] * head(z,nm1) / L[n,n] : 0);
                target += log1m(u_star);
                z[n] = inv_Phi(u_star + u[n] - u_star * u[n]);
            }
            return z;
    } // simplified from https://discourse.mc-stan.org/t/multi-normal-lcdf/6652/6
    matrix transform_tMVN_lp(matrix u, matrix L) {
        int K = rows(u);
        int N = cols(u);
        matrix[K,N] z;
            for(n in 1:N) {
                int nm1 = n - 1;
                vector[K] u_star;
                if(n > 1) {
                    u_star = Phi(z[,1:nm1] * L[1:nm1,n] / L[n,n]);
                } else {
                    u_star = rep_vector(0.5,K);
                }
                z[,n] = inv_Phi(u_star + u[,n] - u_star .* u[,n]);
                target += log1m(u_star);
            }
        return z * L;
    } // vectorized and transposed from above
    matrix nearest_special(matrix x) {
        matrix[rows(x),cols(x)] U = svd_U(x);
        matrix[cols(x),cols(x)] V_prime = svd_V(x)';
        matrix[rows(x),cols(x)] O = U * V_prime;
        matrix[rows(x),cols(x)] Q = qr_thin_Q(O);
        if(prod(O[1,].* Q[1,]) < 0) {
            U[,cols(U)] = -U[,cols(U)];
            return(U * diag_post_multiply(V_prime, sqrt(columns_dot_self(x))));
        } else {
            return(diag_post_multiply(O, sqrt(columns_dot_self(x))));
        }
    } // https://math.stackexchange.com/a/3481936
    matrix nearest_ps(matrix x) {
        matrix[rows(x),cols(x)] O = svd_U(x) * svd_V(x)';
        matrix[cols(x),cols(x)] Q = qr_thin_Q(O);
        for(i in 1:cols(x)) {
            if((O[1,i] * Q[1,i]) < 0) {
                O[,i:] = rep_matrix(0, rows(O), cols(O)-i+1);
                break;
            }
        }
        return(diag_post_multiply(O, sqrt(columns_dot_self(x))));
    } // https://scicomp.stackexchange.com/questions/30631/how-to-find-the-nearest-a-near-positive-definite-from-a-given-matrix
    matrix mean_orthogonal_points(matrix x) {
        matrix[rows(x),cols(x)] O = svd_U(x) * svd_V(x)';
        matrix[rows(x),cols(x)] weighted_mean;
        for(i in 1:cols(x)) {
            vector[cols(x)] inv_dists;
            for(j in 1:cols(x)) {
                inv_dists[j] = inv(squared_distance(x[,i], O[,j]));
            }
            weighted_mean[,i] = O * (inv_dists / sum(inv_dists));
        }
        return(diag_post_multiply(weighted_mean, sqrt(columns_dot_self(x))));
    } // https://scicomp.stackexchange.com/questions/30631/how-to-find-the-nearest-a-near-positive-definite-from-a-given-matrix
    matrix mean_special_orthogonal_points(matrix x) {
        matrix[rows(x),cols(x)] Q = diag_post_multiply(qr_thin_Q(svd_U(x) * svd_V(x)'), sqrt(columns_dot_self(x)));
        matrix[rows(x),cols(x)] weighted_mean;
        for(i in 1:cols(x)) {
            vector[cols(x)] inv_dists = inv(columns_dot_self(rep_matrix(x[,i],cols(x)) - Q))';
            weighted_mean[,i] = Q * (inv_dists / sum(inv_dists));
        }
        return(weighted_mean);
    } // https://scicomp.stackexchange.com/questions/30631/how-to-find-the-nearest-a-near-positive-definite-from-a-given-matrix
}

