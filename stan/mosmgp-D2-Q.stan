// Two-band, Q-component MOSM GP, heterotopic time-series
//
// k_ij(tau) = sum_q alpha_ij^(q) exp(-Sigma_ij^(q) (tau + theta_ij^(q))^2 / 2)
//                    * cos(mu_ij^(q) (tau + theta_ij^(q)) + phi_ij^(q))
//
// For each component q:
//   theta_ij = theta_i - theta_j,  phi_ij = phi_i - phi_j
//   mu_ij    = (Sigma_j mu_i + Sigma_i mu_j) / (Sigma_i + Sigma_j)
//   Sigma_ij = 2 Sigma_i Sigma_j / (Sigma_i + Sigma_j)
//   w_ij     = w_i w_j exp(-(mu_i - mu_j)^2 / (4 (Sigma_i + Sigma_j)))
//   alpha_ij = w_ij (2 pi)^(1/2) Sigma_ij^(1/2)
//
// Covariance of the stacked vector [y1; y2] (length N1 + N2):
//   K11 (N1 x N1), K22 (N2 x N2), K12 (N1 x N2) with K12[n, m] = k_12(t1[n] - t2[m]),
//   K21 = K12'.
//
// Identifiability:
//  - Only theta_ij, phi_ij matter, so theta_1^(q) = phi_1^(q) = 0 for all q.
//  - Label switching between components is prevented by ordering channel-1
//    means: mu1[1] < mu1[2] < ... < mu1[Q].

functions {
  // Per-component kernel parameters for the channel pair (i, j).
  // Columns: 1 alpha_ij, 2 Sigma_ij, 3 mu_ij, 4 theta_ij, 5 phi_ij
  matrix pair_pars(int i, int j, matrix w, matrix mu, matrix Sigma,
                   matrix theta, matrix phi) {
    int Q = rows(w);
    matrix[Q, 5] p;
    for (q in 1:Q) {
      real Sigma_ij = 2 * Sigma[q, i] * Sigma[q, j] / (Sigma[q, i] + Sigma[q, j]);
      real mu_ij = (Sigma[q, j] * mu[q, i] + Sigma[q, i] * mu[q, j])
                   / (Sigma[q, i] + Sigma[q, j]);
      real w_ij = w[q, i] * w[q, j]
                  * exp(-square(mu[q, i] - mu[q, j]) / (4 * (Sigma[q, i] + Sigma[q, j])));
      p[q, 1] = w_ij * sqrt(2 * pi() * Sigma_ij);
      p[q, 2] = Sigma_ij;
      p[q, 3] = mu_ij;
      p[q, 4] = theta[q, i] - theta[q, j];
      p[q, 5] = phi[q, i] - phi[q, j];
    }
    return p;
  }

  // Sum over components of the cross-covariance at lag tau
  real mosm_k(real tau, matrix p) {
    real k = 0;
    for (q in 1:rows(p)) {
      real tt = tau + p[q, 4];
      k += p[q, 1] * exp(-0.5 * p[q, 2] * square(tt)) * cos(p[q, 3] * tt + p[q, 5]);
    }
    return k;
  }
}

data {
  int<lower=1> Q; // no. spectral components
  int<lower=1> N1;
  int<lower=1> N2;
  vector[N1] t1;
  vector[N2] t2;
  vector[N1] y1;
  vector[N2] y2;
  vector<lower=0>[N1] y1_se;
  vector<lower=0>[N2] y2_se;

  real<lower=0> ell_min; // min. envelope width
  real<lower=ell_min> ell_max; // max. envelope width
  real ell_loc; // lognormal prior location for ell (log scale)
  real<lower=0> ell_scale; // lognormal prior scale for ell
  real<lower=0> v_scale; // prior scale for component variances
  real<lower=0> mu_scale; // prior scale for frequencies (radians/unit time)
}

transformed data {
  int M = N1 + N2;
  vector[M] y = append_row(y1, y2);
  vector[M] y_se = append_row(y1_se, y2_se);
  vector[M] zero = rep_vector(0, M);
  real epsilon = 1e-6;
}

parameters {
  matrix<lower=0>[Q, 2] v; // variance of component q in band i
  positive_ordered[Q] mu1; // band 1 spectral means (ordered)
  vector<lower=0>[Q] mu2;
  matrix<lower=ell_min, upper=ell_max>[Q, 2] ell; // envelope widths
  vector[Q] theta2; // delay of band 2 relative to 1
  vector<lower=-pi(), upper=pi()>[Q] phi2; // phase of band 2 relative to 1
}

transformed parameters {
  matrix[Q, 2] Sigma = inv_square(ell); // spectral variances
  matrix[Q, 2] mu = append_col(mu1, mu2);
  matrix[Q, 2] theta = append_col(rep_vector(0, Q), theta2);
  matrix[Q, 2] phi = append_col(rep_vector(0, Q), phi2);
  matrix[Q, 2] w;

  // weights derived from component variances
  for (i in 1:2) {
    for (q in 1:Q) {
      w[q, i] = sqrt(v[q, i] / sqrt(2 * pi() * Sigma[q, i]));
    }
  }
}

model {
  to_vector(v) ~ normal(0, v_scale); // half normal
  mu1 ~ normal(0, mu_scale);
  mu2 ~ normal(0, mu_scale);
  to_vector(ell) ~ lognormal(ell_loc, ell_scale);
  theta2 ~ normal(0, 1);
  phi2 ~ uniform(-pi(), pi());

  // Per-component kernel parameters for each distinct channel pair
  matrix[Q, 5] p11 = pair_pars(1, 1, w, mu, Sigma, theta, phi);
  matrix[Q, 5] p22 = pair_pars(2, 2, w, mu, Sigma, theta, phi);
  matrix[Q, 5] p12 = pair_pars(1, 2, w, mu, Sigma, theta, phi);

  matrix[M, M] K;

  // Block K11 (N1 x N1), symmetric
  for (n in 1:N1) {
    for (m in n:N1) {
      real k = mosm_k(t1[n] - t1[m], p11);
      K[n, m] = k;
      K[m, n] = k;
    }
  }
  // Block K22 (N2 x N2), symmetric
  for (n in 1:N2) {
    for (m in n:N2) {
      real k = mosm_k(t2[n] - t2[m], p22);
      K[N1 + n, N1 + m] = k;
      K[N1 + m, N1 + n] = k;
    }
  }
  // Cross-block K12 (N1 x N2) and its transpose K21
  for (n in 1:N1) {
    for (m in 1:N2) {
      real k = mosm_k(t1[n] - t2[m], p12);
      K[n, N1 + m] = k;
      K[N1 + m, n] = k;
    }
  }

  K = add_diag(K, square(y_se));
  K = add_diag(K, epsilon);

  y ~ multi_normal_cholesky(zero, cholesky_decompose(K));
}
