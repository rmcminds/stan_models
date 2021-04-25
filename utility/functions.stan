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
        vector[rows(Z1)] delta = alpha ./ sqrt(1 - alpha^2);
        matrix[rows(Z1),cols(Z1)] Z
            = diag_pre_multiply(inv_sqrt(1 - square(delta) * 2 / pi()),
                                diag_pre_multiply(delta ./ alpha,
                                                  Z1 + diag_pre_multiply(alpha, Z2))
                                - rep_matrix(delta * sqrt(2 / pi()), cols(Z1)));
        return(Z);
    } // skew-normal matrix with mean 0 and sd 1, assuming Z1 is normal, Z2 is half-normal, and each have location 0 and scale 1
}

