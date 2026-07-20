functions {
  // Function returns a posterior predictive sample given hyperparameter values
  vector gp_pred_rng(vector x_star,
                     vector x,
                     vector y,
                     vector y_se,
                     int D,
                     array[] int ns,
                     array[] int ns_star,
                     vector w,
                     vector Sigma,
                     vector mu,
                     vector theta,
                     real epsilon) {

    int N = rows(y);
    int N_star = size(x_star);
    vector[N_star] f_star;

    {
      matrix[N, N] L;
      vector[N] alpha;

      matrix[N, N_star] v;
      vector[N_star] fstar_mu;
      matrix[N_star, N_star] fstar_Sigma;

      // N x N covariance KS = K(X,X) + Sigma_noise
      matrix[N, N] KS;
      KS = K_mat(x, x, D, ns, ns, w, Sigma, mu, theta);

      for (r in 1:N) {
        for (c in r:N) {
          KS[c,r] = KS[r,c];
        }
      }
      KS = add_diag(KS, square(y_se));
      KS = add_diag(KS, epsilon);

      // N x N* covariance K* = K(X,X*)
      matrix[N, N_star] K_star;
      K_star = K_mat(x, x_star, D, ns, ns_star, w, Sigma, mu, theta);

      // N* x N* covariance K** = K(X*,X*)
      matrix[N_star, N_star] K_starstar;
      K_starstar = K_mat(x_star, x_star, D, ns_star, ns_star, w, Sigma, mu, theta);

      for (r in 1:N_star) {
        for (c in r:N_star) {
          K_starstar[c,r] = K_starstar[r,c];
        }
      }
      K_starstar = add_diag(K_starstar, epsilon);

      L = cholesky_decompose(KS);
      alpha = mdivide_left_tri_low(L, y);
      alpha = mdivide_right_tri_low(alpha', L)';

      fstar_mu = K_star' * alpha;

      v = mdivide_left_tri_low(L, K_star);
      fstar_Sigma = K_starstar - v' * v;

      f_star = multi_normal_cholesky_rng(
        fstar_mu,
        cholesky_decompose(add_diag(fstar_Sigma, epsilon))
      );
    }
    return f_star;
  }
  // Returns the cross-covariance between x in band i and x' in band j
  real k_ij(real x,
            real x_prime,
            real w_ij,
            real Sigma_ij,
            real mu_ij,
            real theta_ij) {
    real tau = abs(x - x_prime);
    real alpha_ij = w_ij * (2.0*pi())^(1/2.0) * sqrt(abs(Sigma_ij));

    return(
      alpha_ij *
      exp(-0.5 * (tau + theta_ij)^2 * Sigma_ij) *
      cos((tau + theta_ij) * mu_ij)
    );
  }
  // Returns the autocovariance between x and x' within band i
  real k_ii(real x,
            real x_prime,
            real w_i,
            real Sigma_i,
            real mu_i) {
    real tau = abs(x - x_prime);
    real alpha_ii = w_i^2.0 * (2.0*pi())^(1/2.0) * sqrt(abs(Sigma_i));

    return( alpha_ii * exp(-0.5*tau^2 * Sigma_i) * cos(tau*mu_i) );
  }
  // Returns the cross-covariance matrix between locations x and x' in
  // bands i and j respectively
  matrix K_ij(vector x_i,
              vector x_j,
              int i,
              int j,
              vector ws,
              vector Sigmas,
              vector mus,
              vector thetas) {

    int n_i = size(x_i);
    int n_j = size(x_j);

    matrix[n_i, n_j] Kij;

    for (r in 1:n_i) {
      for (c in 1:n_j) {

        if (i == j) { // Covariance matrix n_i x n_i
          Kij[r][c] = k_ii(x_i[r], x_j[c], ws[i], Sigmas[i], mus[i]);
        }
        else { // Cross-covariance matrix n_i x n_j
          real w_ij = ws[i]*ws[j] *
            exp( -1.0/4*(mus[i]-mus[j]) *
                 (Sigmas[i]+Sigmas[j])^-1 *
                 (mus[i]-mus[j])
                );
          real Sigma_ij = 2*Sigmas[i] * (Sigmas[i] + Sigmas[j])^-1 * Sigmas[j];
          real mu_ij = (Sigmas[i] + Sigmas[j])^-1 *
                       (Sigmas[i]*mus[j] + Sigmas[j]*mus[i]);
          real theta_ij = thetas[i] - thetas[j];

          Kij[r][c] = k_ij(x_i[r], x_j[c], w_ij, Sigma_ij, mu_ij, theta_ij);
        }
      }
    }
    return(Kij);
  }
  // Returns the multi-band cross-variance matrix, comprising D x D
  // K_ij submatrices, evaluated between vectors of points x and x_prime
  //
  // if x == x_prime, then returns symmetrical cross-covariance matrix

  matrix K_mat(vector x,
               vector x_prime,
               int D,
               array[] int ns,
               array[] int ns_prime,
               vector ws,
               vector Sigmas,
               vector mus,
               vector thetas) {

    array[D] int row_end_idx = cumulative_sum(ns);
    array[D] int col_end_idx = cumulative_sum(ns_prime);
    array[D] int row_start_idx = to_int(
                                   to_array_1d(1 + (
                                     to_vector(cumulative_sum(ns)) -
                                     to_vector(ns))
                                   )
                                 );
    array[D] int col_start_idx = to_int(
                                    to_array_1d(1 + (
                                      to_vector(cumulative_sum(ns_prime)) -
                                      to_vector(ns_prime))
                                   )
                                 );

    int N = size(x);
    int N_prime = size(x_prime);

    matrix[N, N_prime] Kmat;

    for (i in 1:D) {
      for (j in 1:D) {

        int nrows = ns[i];
        int ncols = ns_prime[j];

        matrix[nrows, ncols] K_ij_mat;

        int start_row = row_start_idx[i];
        int end_row = row_end_idx[i];
        int start_col = col_start_idx[j];
        int end_col = col_end_idx[j];

        K_ij_mat = K_ij(x[start_row:end_row],
                        x_prime[start_col:end_col],
                        i, j, ws, Sigmas, mus, thetas);

        Kmat[start_row:end_row, start_col:end_col] = K_ij_mat;
      }
    }
    return(Kmat);
  }
}
data {
  int<lower=1> D; // no. bands
  int<lower=1> N; // no. observations
  array[D] int<lower=0> ns; // no. observations in each band
  array[N] int<lower=1, upper=D> d; // band of each observation
  vector[N] x;
  vector[N] y;
  vector[N] y_se;

  int<lower=1> N_star; // no. target locations
  array[D] int<lower=0> ns_star; // no. targets in each band
  array[N_star] int<lower=1, upper=D> d_star; // band of each target location
  vector[N_star] x_star; // target locations
}
transformed data {
  real epsilon = 1e-9; // jitter
  vector[D] Sigma = [0.02, 0.02]';
}
parameters {
  vector<lower=0>[D] w;
  positive_ordered[D] mu;
  vector[D] theta;
}
transformed parameters {
  real deltaTheta = theta[1] - theta[2];
}
model {
  w ~ normal(0, 1.0);
  mu ~ normal(0, 1.0);
  theta ~ normal(0, 1.0);

  // N x N covariance KS = K(X,X) + Sigma_noise
  matrix[N, N] KS;
  KS = K_mat(x, x, D, ns, ns, w, Sigma, mu, theta);

  for (r in 1:N) {
    for (c in r:N) {
      KS[c,r] = KS[r,c];
    }
  }
  KS = add_diag(KS, square(y_se));

  y ~ multi_normal_cholesky(
    rep_vector(0, N),
    cholesky_decompose(add_diag(KS, epsilon))
  );
}

generated quantities {
  vector[N_star] f_star;

  {
    matrix[N, N] L;
    vector[N] alpha;

    matrix[N, N_star] v;
    vector[N_star] fstar_mu;
    matrix[N_star, N_star] fstar_Sigma;

    // N x N covariance KS = K(X,X) + Sigma_noise
    matrix[N, N] KS;
    KS = K_mat(x, x, D, ns, ns, w, Sigma, mu, theta);

    for (r in 1:N) {
      for (c in r:N) {
        KS[c,r] = KS[r,c];
      }
    }
    KS = add_diag(KS, square(y_se));
    KS = add_diag(KS, epsilon);

    // N x N* covariance K* = K(X,X*)
    matrix[N, N_star] K_star;
    K_star = K_mat(x, x_star, D, ns, ns_star, w, Sigma, mu, theta);

    // N* x N* covariance K** = K(X*,X*)
    matrix[N_star, N_star] K_starstar;
    K_starstar = K_mat(x_star, x_star, D, ns_star, ns_star, w, Sigma, mu, theta);

    for (r in 1:N_star) {
      for (c in r:N_star) {
        K_starstar[c,r] = K_starstar[r,c];
      }
    }
    K_starstar = add_diag(K_starstar, epsilon);

    L = cholesky_decompose(KS);
    alpha = mdivide_left_tri_low(L, y);
    alpha = mdivide_right_tri_low(alpha', L)';

    fstar_mu = K_star' * alpha;

    v = mdivide_left_tri_low(L, K_star);
    fstar_Sigma = K_starstar - v' * v;

    f_star = multi_normal_cholesky_rng(
        fstar_mu,
        cholesky_decompose(add_diag(fstar_Sigma, epsilon))
    );
  }
}
