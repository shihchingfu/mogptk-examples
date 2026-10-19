// Two-band, single-component multi-output spectral mixture GP (MOSM)
// One-dimensional heterotopic time series
// Covariance of the stacked vector [y1; y2] (length N1 + N2):
//   K11 (N1 x N1), K22 (N2 x N2),
//   K12 (N1 x N2) with K12[n, m] = k_12(t1[n] - t2[m]),
//   K21 = K12'.
functions {
  real mosm_k(real tau, int i, int j, vector w, vector mu, vector Sigma,
              vector theta, vector phi) {
    real theta_ij = theta[i] - theta[j];
    real phi_ij = phi[i] - phi[j];
    real mu_ij = (Sigma[j] * mu[i] + Sigma[i] * mu[j]) / (Sigma[i] + Sigma[j]);
    real Sigma_ij = 2 * Sigma[i] * Sigma[j] / (Sigma[i] + Sigma[j]);
    real w_ij = w[i] * w[j]
                * exp(-square(mu[i] - mu[j]) / (4 * (Sigma[i] + Sigma[j])));
    real alpha_ij = w_ij * sqrt(2 * pi() * Sigma_ij);
    real tt = tau + theta_ij;
    return alpha_ij * exp(-0.5 * Sigma_ij * square(tt)) * cos(mu_ij * tt + phi_ij);
  }

  matrix K_mat(vector t1, vector t2, vector w, vector mu, vector Sigma,
               vector theta, vector phi) {

    int N1 = rows(t1);
    int N2 = rows(t2);
    int M = N1 + N2;
    matrix[M, M] K;

    // Block K11 (N1 x N1), symmetric
    for (n in 1:N1) {
      for (m in n:N1) {
        real k = mosm_k(t1[n] - t1[m], 1, 1, w, mu, Sigma, theta, phi);
        K[n, m] = k;
        K[m, n] = k;
      }
    }
    // Block K22 (N2 x N2), symmetric
    for (n in 1:N2) {
      for (m in n:N2) {
        real k = mosm_k(t2[n] - t2[m], 2, 2, w, mu, Sigma, theta, phi);
        K[N1 + n, N1 + m] = k;
        K[N1 + m, N1 + n] = k;
      }
    }
    // Cross-block K12 (N1 x N2) and its transpose K21
    for (n in 1:N1) {
      for (m in 1:N2) {
        real k = mosm_k(t1[n] - t2[m], 1, 2, w, mu, Sigma, theta, phi);
        K[n, N1 + m] = k;
        K[N1 + m, n] = k;
      }
    }
    return(K);
  }
}

data {
  int<lower=1> N1;
  int<lower=1> N2;
  vector[N1] t1;
  vector[N2] t2;
  vector[N1] y1;
  vector[N2] y2;
  vector[N1] y1_se;
  vector[N2] y2_se;
}

transformed data {
  int M = N1 + N2;
  vector[M] y = append_row(y1, y2);
  vector[M] y_se = append_row(y1_se, y2_se);
  vector[M] zero = rep_vector(0, M);
  real epsilon = 1e-6;
}

parameters {
  vector<lower=0>[2] w; // amplitudes
  vector<lower=0>[2] mu; // spectral means (radians/unit time)
  vector<lower=0>[2] Sigma; // spectral variances
  real theta2; // delay of band 2 relative to 1
  real<lower=-pi(), upper=pi()> phi2; // phase of band 2 relative to 1
}

transformed parameters {
  vector[2] theta = [0, theta2]';
  vector[2] phi = [0, phi2]';
}

model {
  w ~ uniform(0, 10);
  mu ~ uniform(0, 5);
  Sigma ~ uniform(0, 5);
  theta2 ~ uniform(0, 1);
  phi2 ~ uniform(-pi(), pi());

  matrix[M, M] K = K_mat(t1, t2, w, mu, Sigma, theta, phi);

  K = add_diag(K, square(y_se));
  K = add_diag(K, epsilon);

  y ~ multi_normal_cholesky(zero, cholesky_decompose(K));

}
