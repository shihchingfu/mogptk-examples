functions {
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
  matrix K_ij(vector x,
              vector x_prime,
              int i,
              int j,
              vector ws,
              vector Sigmas,
              vector mus,
              vector thetas) {

    int n = size(x);
    int n_prime = size(x_prime);

    matrix[n, n_prime] Kij;

    for (r in 1:n) {
      for (c in 1:n_prime) {

        if (i == j) { // Covariance matrix n_i x n_i
          Kij[r,c] = k_ii(x[r], x_prime[c], ws[i], Sigmas[i], mus[i]);
        }
        else { // Cross-covariance matrix n_i x n_j
          real w_ij = ws[i]*ws[j] * exp( -1.0/4*(mus[i]-mus[j]) *
                 (Sigmas[i]+Sigmas[j])^-1 * (mus[i]-mus[j]) );
          real Sigma_ij = 2*Sigmas[i] * (Sigmas[i] + Sigmas[j])^-1 * Sigmas[j];
          real mu_ij = (Sigmas[i] + Sigmas[j])^-1 *
                       (Sigmas[i]*mus[j] + Sigmas[j]*mus[i]);
          real theta_ij = thetas[i] - thetas[j];

          Kij[r,c] = k_ij(x[r], x_prime[c], w_ij, Sigma_ij, mu_ij, theta_ij);
        }
      }
    }
    return(Kij);
  }
  // Returns the square multi-band cross-variance matrix, comprising D x D
  // K_ij submatrices, evaluated between vectors of points x with itself
  matrix Kxx_mat(vector x,
               int D,
               array[] int ns,
               vector ws,
               vector Sigmas,
               vector mus,
               vector thetas) {

    array[D] int end_idx = cumulative_sum(ns);
    array[D] int start_idx = to_int(
      to_array_1d(1 + (
        to_vector(cumulative_sum(ns)) -
        to_vector(ns))
      )
    );

    int N = size(x);

    matrix[N, N] Kmat;

    for (i in 1:D) {
      for (j in 1:D) {

        int nrows = ns[i];
        int ncols = ns[j];

        matrix[nrows, ncols] K_ij_mat;

        int start_row = start_idx[i];
        int end_row = end_idx[i];
        int start_col = start_idx[j];
        int end_col = end_idx[j];

        K_ij_mat = K_ij(x[start_row:end_row],
                        x[start_col:end_col],
                        i, j, ws, Sigmas, mus, thetas);

        Kmat[start_row:end_row, start_col:end_col] = K_ij_mat;
      }
    }

    for (r in 1:N) {
      for (c in r:N) {
        Kmat[c,r] = Kmat[r,c];
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
}
transformed data {
  real epsilon = 1e-9; // jitter
}
parameters {
  vector<lower=0>[D] w;
  vector<lower=0>[D] Sigma;
  positive_ordered[D] mu;
  vector[D] theta;
}
transformed parameters {
  real deltaTheta = theta[1] - theta[2];
}
model {
  w ~ normal(0, 10.0);
  Sigma ~ normal(0, 1.0);
  mu ~ normal(0, 1.0);
  theta ~ normal(0, 10.0);

  // N x N covariance KS = K(X,X) + Sigma_noise
  matrix[N, N] KS;
  KS = Kxx_mat(x, D, ns, w, Sigma, mu, theta);
  KS = add_diag(KS, square(y_se));

  y ~ multi_normal_cholesky(
    rep_vector(0, N),
    cholesky_decompose(add_diag(KS, epsilon))
  );
}

generated quantities {
}
