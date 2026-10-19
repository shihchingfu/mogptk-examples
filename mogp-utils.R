library(mvtnorm)
library(tictoc)
library(stringr)

# Cross-band covariance
#
# Within-band i parameters
#   w_i: weight
#   mu_i: spectral mean (frequency) [radians per unit time, i.e., period = 2*pi/mu]
#   Sigma_i: spectral variance
#   theta_i: delay
#   phi_i: phase
#
# Cross-band quantities (i, j):
#   theta_ij = theta_i - theta_j [offset]
#   phi_ij   = phi_i - phi_j [phase offset, between -/+ pi]
#   mu_ij    = (Sigma_j mu_i + Sigma_i mu_j) / (Sigma_i + Sigma_j) [mean frequency of cross-covariance]
#   Sigma_ij = 2 Sigma_i Sigma_j / (Sigma_i + Sigma_j) [envelope width  ]
#   w_ij     = w_i w_j exp(-(mu_i - mu_j)^2 / (4 (Sigma_i + Sigma_j)))
#   alpha_ij = w_ij (2 pi)^(1/2) Sigma_ij^(1/2) [amplitude]
#
# Kernel (tau = t - t'):
#   k_ij(tau) = alpha_ij exp(-1/2 Sigma_ij (tau + theta_ij)^2) cos(mu_ij (tau + theta_ij) + phi_ij)
k_ij <- function(tau, i, j, w, mu, Sigma, theta, phi) {
  theta_ij <- theta[i] - theta[j]
  phi_ij <- phi[i] - phi[j]
  mu_ij <- (Sigma[j] * mu[i] + Sigma[i] * mu[j]) / (Sigma[i] + Sigma[j])
  Sigma_ij <- 2 * Sigma[i] * Sigma[j] / (Sigma[i] + Sigma[j])
  w_ij <- w[i] * w[j] * exp(-(mu[i] - mu[j])^2 / (4 * (Sigma[i] + Sigma[j])))
  alpha_ij <- w_ij * sqrt(2 * pi * Sigma_ij)
  alpha_ij *
    exp(-0.5 * Sigma_ij * (tau + theta_ij)^2) *
    cos(mu_ij * (tau + theta_ij) + phi_ij)
}

# Submatrix where k_ij(ta[n] - tb[m]) = cov(f_i(ta[n]), f_j(tb[m]))
Kmat_ij <- function(ta, tb, i, j, w, mu, Sigma, theta, phi) {
  k_ij(outer(ta, tb, "-"), i, j, w, mu, Sigma, theta, phi)
}

# Two-band, single component MOSK GP, covariance matrix with heterotopic inputs
K_D2_Q1 <- function(t1, t2, w, mu, Sigma, theta, phi) {
  N1 <- length(t1)
  N2 <- length(t2)

  K11 <- Kmat_ij(t1, t1, 1, 1, w, mu, Sigma, theta, phi) # N1 x N1
  K22 <- Kmat_ij(t2, t2, 2, 2, w, mu, Sigma, theta, phi) # N2 x N2
  K12 <- Kmat_ij(t1, t2, 1, 2, w, mu, Sigma, theta, phi) # N1 x N2

  K <- rbind(
    cbind(K11, K12),
    cbind(t(K12), K22)
  )
  K <- (K + t(K)) / 2
}

# Simulate a single draw from an MOSK GP (D = 2, Q = 1)
sim_D2_Q1 <- function(
  tstart = 0,
  tend = 10,
  n_per_band = 100,
  masked_pct = 0,
  w = NULL,
  mu = NULL,
  Sigma = NULL,
  theta = NULL,
  phi = NULL,
  noise = c(0, 0), # observation noise std
  seed = round(runif(1, min = 0, max = 9999))
) {
  if (is.null(w)) {
    w <- round(sqrt(runif(2, min = 0, max = 5)), 2)
  } # amplitude (acf scale is w_i^2)

  if (is.null(mu)) {
    mu1 <- round((runif(1, min = 0, max = 10)), 2) # period = 2*pi/mu
    mu2 <- round(rnorm(1, mu1, sd = 1), 2)
    mu <- c(mu1, mu2)
  }

  if (is.null(Sigma)) {
    Sigma <- round(runif(2, min = 0, max = 5), 2)^2
  } # spectral variance

  if (is.null(theta)) {
    theta <- c(0.0, round(runif(1, min = 0, max = 5), 2))
  }

  if (is.null(phi)) {
    phi <- c(0.0, round(runif(1, min = -pi, max = pi), 2))
  }

  # Covariance matrix of the stacked vector [y1(t); y2(t)]
  t1 <- seq(tstart, tend, length.out = n_per_band)
  t2 <- seq(tstart, tend, length.out = n_per_band)
  N1 <- length(t1)
  N2 <- length(t2)

  K <- K_D2_Q1(t1, t2, w, mu, Sigma, theta, phi) +
    diag(c(rep(noise[1]^2, N1), rep(noise[2]^2, N2))) +
    diag(1e-6, N1 + N2) # jitter

  # Sample
  set.seed(seed)
  y <- t(chol(K)) %*% rnorm(N1 + N2)

  output_df <- data.frame(
    t = c(t1, t2),
    y = y,
    y_se = rep(noise, times = c(N1, N2)),
    d = factor(rep(c(1, 2), times = c(N1, N2))),
    masked = create_block_missing_mask(N1 + N2, masked_pct, seed)
  )

  params_df = data.frame(
    seed = seed,
    w = w,
    mu = mu,
    Sigma = Sigma,
    theta = theta,
    phi = phi
  )

  return(list(params_df, output_df))
}

# Draw single posterior predictive sample
pred_D2_Q1_draw <- function(
  posterior_draw, # 1-row matrix
  t1,
  t2,
  y1,
  y2,
  y1_se,
  y2_se,
  ts1,
  ts2,
  jitter = 1e-6
) {
  N1 <- length(t1)
  N2 <- length(t2)
  S1 <- length(ts1)
  S2 <- length(ts2)

  w <- posterior_draw[c("w[1]", "w[2]")]
  mu <- posterior_draw[c("mu[1]", "mu[2]")]
  Sigma <- posterior_draw[c("Sigma[1]", "Sigma[2]")]
  theta <- posterior_draw[c("theta[1]", "theta[2]")]
  phi <- posterior_draw[c("phi[1]", "phi[2]")]

  # Training covariance of [y1; y2], including observation noise
  Kxx12 <- Kmat_ij(t1, t2, 1, 2, w, mu, Sigma, theta, phi)
  Kxx11 <- Kmat_ij(t1, t1, 1, 1, w, mu, Sigma, theta, phi)
  Kxx22 <- Kmat_ij(t2, t2, 2, 2, w, mu, Sigma, theta, phi)

  Kxx <- rbind(
    cbind(Kxx11, Kxx12),
    cbind(t(Kxx12), Kxx22)
  )
  Kxx <- (Kxx + t(Kxx)) / 2 + diag(c(y1_se, y2_se) + jitter)

  # Test-train cross-covariance: rows [f1*; f2*], columns [y1; y2]
  Kstar11 <- Kmat_ij(ts1, t1, 1, 1, w, mu, Sigma, theta, phi)
  Kstar12 <- Kmat_ij(ts1, t2, 1, 2, w, mu, Sigma, theta, phi)
  Kstar21 <- Kmat_ij(ts2, t1, 2, 1, w, mu, Sigma, theta, phi)
  Kstar22 <- Kmat_ij(ts2, t2, 2, 2, w, mu, Sigma, theta, phi)
  Kstar <- rbind(
    cbind(Kstar11, Kstar12),
    cbind(Kstar21, Kstar22)
  )

  # Test prior covariance of [f1*; f2*] (noise-free)
  Kstarstar11 <- Kmat_ij(ts1, ts1, 1, 1, w, mu, Sigma, theta, phi)
  Kstarstar12 <- Kmat_ij(ts1, ts2, 1, 2, w, mu, Sigma, theta, phi)
  Kstarstar22 <- Kmat_ij(ts2, ts2, 2, 2, w, mu, Sigma, theta, phi)
  Kstarstar <- rbind(
    cbind(Kstarstar11, Kstarstar12),
    cbind(t(Kstarstar12), Kstarstar22)
  )

  # Conditional mean and covariance via Cholesky solves
  L <- t(chol(Kxx))
  A <- forwardsolve(L, t(Kstar)) # (N1+N2) x (S1+S2)
  v <- forwardsolve(L, c(y1, y2))
  m <- as.vector(crossprod(A, v)) # conditional mean
  C <- Kstarstar - crossprod(A) # conditional covariance
  C <- (C + t(C)) / 2

  f <- rmvn(m, C) # joint draw of [f1*; f2*]

  data.frame(
    ts = c(ts1, ts2),
    fs = f,
    d = rep(c(1, 2), times = c(S1, S2))
  )
}

# Draw multiple posterior predictive samples
pred_D2_Q1 <- function(
  draws, # matrix
  t1,
  t2,
  y1,
  y2,
  y1_se,
  y2_se,
  ts1,
  ts2 = ts1,
  n_draws = 200,
  seed = NULL
) {
  if (!is.null(seed)) {
    set.seed(seed)
  }

  #draws <- as.data.frame(draws)
  idx <- if (n_draws >= nrow(draws)) {
    seq_len(nrow(draws))
  } else {
    sort(sample(nrow(draws), n_draws))
  }

  S1 <- length(ts1)
  S2 <- length(ts2)
  D <- length(idx)

  out_list <- vector("list", n_draws)

  for (d in seq_along(idx)) {
    r <- draws[idx[d], ]

    pr <- pred_D2_Q1_draw(
      r,
      t1,
      t2,
      y1,
      y2,
      y1_se,
      y2_se,
      ts1,
      ts2
    ) |>
      mutate(.draw = d)

    out_list[[d]] <- pr
  }
  out_df <- bind_rows(out_list, .id = ".draw")

  return(out_df)
}

# Sample from multi-variate normal where C might be numerically singular with small
# negative eigenvalues
rmvn <- function(mu, covmat) {
  for (jitter in 1e-6 * 10^(0:3)) {
    L <- tryCatch(
      {
        chol(covmat + diag(jitter, nrow(covmat))) # add jitter to diagonal to make stable
      },
      error = function(e) NULL
    )
    if (!is.null(L)) {
      return(mu + as.vector(t(L) %*% rnorm(length(mu))))
    }
  }
  e <- eigen(covmat, symmetric = TRUE)
  # replace negative eigenvalues with zero
  mu + as.vector(e$vectors %*% (sqrt(pmax(e$values, 0)) * rnorm(length(mu))))
}

# Returns a logical mask for dropping values in contiguous blocks
create_block_missing_mask <- function(n, prop_missing = 0.5, seed = NULL) {
  set.seed(seed)

  n_missing <- floor(prop_missing * n)
  idx <- 1:n

  for (i in 1:(prop_missing * n)) {
    lo_cutpoint <- runif(1, min = 1, max = n_missing)
    up_cutpoint <- mean(c(lo_cutpoint, n))
    offset <- ceiling(up_cutpoint - lo_cutpoint)

    inc_idx <- which(idx >= lo_cutpoint & idx < up_cutpoint)
    dec_idx <- which(idx >= up_cutpoint)

    idx[inc_idx] <- idx[inc_idx] + offset
    idx[dec_idx] <- idx[dec_idx] - offset
  }

  missing_mask <- rep(c(TRUE, FALSE), times = c(n_missing, n - n_missing))
  missing_mask <- missing_mask[idx]

  return(missing_mask)
}
