# Two-band, Q-component MOSM GP ############

# k_ij(tau) = \sum_q alpha_ij^(q) exp(-1/2 Sigma_ij^(q) (tau + theta_ij^(q))^2)
#                      * cos(mu_ij^(q) (tau + theta_ij^(q)) + phi_ij^(q))
#
# Parameters expressed as Q x 2 matrices, each row a component, each column a band.).
#
# Per component q and channel pair (i, j):
#   theta_ij = theta_i - theta_j,   phi_ij = phi_i - phi_j
#   mu_ij    = (Sigma_j mu_i + Sigma_i mu_j) / (Sigma_i + Sigma_j)
#   Sigma_ij = 2 Sigma_i Sigma_j / (Sigma_i + Sigma_j)
#   w_ij     = w_i w_j exp(-(mu_i - mu_j)^2 / (4 (Sigma_i + Sigma_j)))
#   alpha_ij = w_ij (2 pi)^(1/2) Sigma_ij^(1/2)
k_ij_D2_Q <- function(tau, i, j, p) {
  Q <- nrow(p$w)

  k <- 0
  for (q in seq_len(Q)) {
    theta_ij <- p$theta[q, i] - p$theta[q, j]
    phi_ij   <- p$phi[q, i] - p$phi[q, j]
    mu_ij    <- (p$Sigma[q, j] * p$mu[q, i] + p$Sigma[q, i] * p$mu[q, j]) /
      (p$Sigma[q, i] + p$Sigma[q, j])
    Sigma_ij <- 2 * p$Sigma[q, i] * p$Sigma[q, j] / (p$Sigma[q, i] + p$Sigma[q, j])
    w_ij     <- p$w[q, i] * p$w[q, j] *
      exp(-(p$mu[q, i] - p$mu[q, j])^2 / (4 * (p$Sigma[q, i] + p$Sigma[q, j])))
    alpha_ij <- w_ij * sqrt(2 * pi * Sigma_ij)
    k <- k + alpha_ij * exp(-0.5 * Sigma_ij * (tau + theta_ij)^2) *
      cos(mu_ij * (tau + theta_ij) + phi_ij)
  }
  k
}

Kmat_ij_D2_Q <- function(ta, tb, i, j, params) {
  tau <- outer(ta, tb, "-")
  k_ij_D2_Q(tau, i, j, params)
}

# Two-band, multi-component MOSK GP, covariance matrix with heterotopic inputs
K_D2_Q <- function(t1, t2, params) {
  K11 <- Kmat_ij_D2_Q(t1, t1, 1, 1, params)
  K22 <- Kmat_ij_D2_Q(t2, t2, 2, 2, params)
  K12 <- Kmat_ij_D2_Q(t1, t2, 1, 2, params)

  K <- rbind(
    cbind(K11, K12),
    cbind(t(K12), K22)
  )
  K <- (K + t(K)) / 2
}

# Q x 2 matrix from columns named e.g. "w[1,1]", "w[2,1]", "w[1,2]", ...
get_param_mat <- function(r, name, Q) {
  cols <- paste0(name, "[", rep(seq_len(Q), 2), ",", rep(1:2, each = Q), "]")
  matrix(as.numeric(r[cols]), Q, 2)
}

# Draw single posterior predictive sample
pred_D2_Q_draw <- function(
    posterior_draw, # 1 row df
    Q,
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

  params <- list(w = get_param_mat(posterior_draw, "w", Q),
                 mu = get_param_mat(posterior_draw, "mu", Q),
                 Sigma = get_param_mat(posterior_draw, "Sigma", Q),
                 theta = get_param_mat(posterior_draw, "theta", Q),
                 phi =  get_param_mat(posterior_draw, "phi", Q))

  # Training covariance of [y1; y2], including observation noise
  Kxx12 <- Kmat_ij_D2_Q(t1, t2, 1, 2, params)
  Kxx11 <- Kmat_ij_D2_Q(t1, t1, 1, 1, params)
  Kxx22 <- Kmat_ij_D2_Q(t2, t2, 2, 2, params)

  Kxx <- rbind(
    cbind(Kxx11, Kxx12),
    cbind(t(Kxx12), Kxx22)
  )
  Kxx <- (Kxx + t(Kxx)) / 2 + diag(c(y1_se, y2_se) + jitter)

  # Test-train cross-covariance: rows [f1*; f2*], columns [y1; y2]
  Kstar11 <- Kmat_ij_D2_Q(ts1, t1, 1, 1, params)
  Kstar12 <- Kmat_ij_D2_Q(ts1, t2, 1, 2, params)
  Kstar21 <- Kmat_ij_D2_Q(ts2, t1, 2, 1, params)
  Kstar22 <- Kmat_ij_D2_Q(ts2, t2, 2, 2, params)
  Kstar <- rbind(
    cbind(Kstar11, Kstar12),
    cbind(Kstar21, Kstar22)
  )

  # Test prior covariance of [f1*; f2*] (noise-free)
  Kstarstar11 <- Kmat_ij_D2_Q(ts1, ts1, 1, 1, params)
  Kstarstar12 <- Kmat_ij_D2_Q(ts1, ts2, 1, 2, params)
  Kstarstar22 <- Kmat_ij_D2_Q(ts2, ts2, 2, 2, params)
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
    d = factor(rep(c(1, 2), times = c(S1, S2)))
  )
}

# Draw multiple posterior predictive samples
pred_D2_Q <- function(
    draws, # df
    Q,
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
  if (!is.null(seed))
    set.seed(seed)

  draws <- as.data.frame(draws)

  idx <- if (n_draws >= nrow(draws)) {
    seq_len(nrow(draws))
  } else {
    sort(sample(nrow(draws), n_draws))
  }

  S1 <- length(ts1);
  S2 <- length(ts2);
  D <- length(idx)

  out_list <- vector("list", n_draws)

  for (d in seq_along(idx)) {
    r <- as.data.frame( draws[idx[d], ] )

    pr <- pred_D2_Q_draw(r, Q, t1, t2, y1, y2, y1_se, y2_se, ts1, ts2) |>
      mutate(.draw = d)

    out_list[[d]] <- pr
  }

  bind_rows(out_list, .id = ".draw")
}

# Simulate a single draw from an MOSM GP (D = 2, Q)
sim_D2_Q <- function(
    tstart = 0,
    tend = 10,
    n_per_band = 100,
    masked_pct = 0,
    Q,
    params = list(
      w = NULL,
      mu = NULL,
      Sigma = NULL,
      theta = NULL,
      phi = NULL),
    noise = c(0, 0), # observation noise std
    seed = round(runif(1, min = 0, max = 9999))
) {
  set.seed(seed)
  p <- params

  if (is.null(p$w)) {
    p$w <- round(
      sqrt(
        matrix(
          runif(2 * Q, min = 0, max = 5),
          ncol = 2,
          dimnames = list(NULL, c("mu1", "mu2")))
      ),
      2
    )
  }

  if (is.null(p$mu)) {
    mu1 <- round((runif(Q, min = 0, max = 5)), 2) # period = 2*pi/mu
    mu2 <- round(rnorm(Q, mu1, sd = 1), 2)

    while(any(mu2 <= 0)) {
      mu2 <- round(rnorm(Q, mu1, sd = 1), 2)
    }

    p$mu <- cbind(mu1, mu2)
  }

  if (is.null(p$Sigma)) {
    Sigma1 <- round(runif(Q, min = 0, max = 5), 2)^2
    Sigma2 <- round(rnorm(Q, mean = Sigma1, sd = 0.5), 2)

    while(any(Sigma2 <= 0)) {
      Sigma2 <- round(rnorm(Q, mean = Sigma1, sd = 0.5), 2)
    }

    p$Sigma <- cbind(Sigma1, Sigma2)
  }

  if (is.null(p$theta)) {
    p$theta <- cbind(theta1 = 0,
                     theta2 = round(runif(Q, min = 0, max = (tend - tstart)/4), 2))
  }

  if (is.null(p$phi)) {
    p$phi <- cbind(phi1 = 0,
                   phi2 = round(runif(Q, min = -pi/2, max = pi/2), 2))
  }

  # Covariance matrix of the stacked vector [y1(t); y2(t)]
  t1 <- seq(tstart, tend, length.out = n_per_band)
  t2 <- seq(tstart, tend, length.out = n_per_band)
  N1 <- length(t1)
  N2 <- length(t2)

  K <- K_D2_Q(t1, t2, p) +
    diag(c(rep(noise[1]^2, N1), rep(noise[2]^2, N2))) +
    diag(1e-6, N1 + N2) # jitter

  # Sample

  y <- t(chol(K)) %*% rnorm(N1 + N2)

  output_df <- data.frame(
    t = c(t1, t2),
    y = y,
    y_se = rep(noise, times = c(N1, N2)),
    d = factor(rep(c(1, 2), times = c(N1, N2))),
    masked = create_block_missing_mask(N1 + N2, masked_pct, seed)
  )

  params_df <- list(
    seed = seed,
    w = p$w,
    mu = p$mu,
    Sigma = p$Sigma,
    theta = p$theta,
    phi = p$phi
  )

  return(list(params_df, output_df))
}


# Two-band, single-component MOSM GP ##################

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
k_ij_D2_Q1 <- function(tau, i, j, p) {

  theta_ij <- p$theta[i] - p$theta[j]
  phi_ij   <- p$phi[i] - p$phi[j]
  mu_ij    <- (p$Sigma[j] * p$mu[i] + p$Sigma[i] * p$mu[j]) / (p$Sigma[i] + p$Sigma[j])
  Sigma_ij <- 2 * p$Sigma[i] * p$Sigma[j] / (p$Sigma[i] + p$Sigma[j])
  w_ij     <- p$w[i] * p$w[j] *
    exp(-(p$mu[i] - p$mu[j])^2 / (4 * (p$Sigma[i] + p$Sigma[j])))
  alpha_ij <- w_ij * sqrt(2 * pi * Sigma_ij)

  alpha_ij * exp(-0.5 * Sigma_ij * (tau + theta_ij)^2) *
    cos(mu_ij * (tau + theta_ij) + phi_ij)
}

# Submatrix where k_ij(ta[n] - tb[m]) = cov(f_i(ta[n]), f_j(tb[m]))
Kmat_ij_D2_Q1 <- function(ta, tb, i, j, params) {
  tau <- outer(ta, tb, "-")
  k_ij_D2_Q1(tau, i, j, params)
}

# Two-band, single component MOSK GP, covariance matrix with heterotopic inputs
K_D2_Q1 <- function(t1, t2, params) {
  K11 <- Kmat_ij_D2_Q1(t1, t1, 1, 1, params)
  K22 <- Kmat_ij_D2_Q1(t2, t2, 2, 2, params)
  K12 <- Kmat_ij_D2_Q1(t1, t2, 1, 2, params)

  K <- rbind(
    cbind(K11, K12),
    cbind(t(K12), K22)
  )
  K <- (K + t(K)) / 2
}

# Draw single posterior predictive sample
pred_D2_Q1_draw <- function(
  posterior_draw, # 1 row df
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

  params <- list(w = as.numeric( posterior_draw[c("w[1]", "w[2]")] ),
                 mu = as.numeric( posterior_draw[c("mu[1]", "mu[2]")] ),
                 Sigma = as.numeric( posterior_draw[c("Sigma[1]", "Sigma[2]")] ),
                 theta = as.numeric( posterior_draw[c("theta[1]", "theta[2]")] ),
                 phi = as.numeric( posterior_draw[c("phi[1]", "phi[2]")] ))

  # Training covariance of [y1; y2], including observation noise
  Kxx11 <- Kmat_ij_D2_Q1(t1, t1, 1, 1, params)
  Kxx12 <- Kmat_ij_D2_Q1(t1, t2, 1, 2, params)
  Kxx22 <- Kmat_ij_D2_Q1(t2, t2, 2, 2, params)

  Kxx <- rbind(
    cbind(Kxx11, Kxx12),
    cbind(t(Kxx12), Kxx22)
  )
  Kxx <- (Kxx + t(Kxx)) / 2 + diag(c(y1_se, y2_se) + jitter)

  # Test-train cross-covariance: rows [f1*; f2*], columns [y1; y2]
  Kstar11 <- Kmat_ij_D2_Q1(ts1, t1, 1, 1, params)
  Kstar12 <- Kmat_ij_D2_Q1(ts1, t2, 1, 2, params)
  Kstar21 <- Kmat_ij_D2_Q1(ts2, t1, 2, 1, params)
  Kstar22 <- Kmat_ij_D2_Q1(ts2, t2, 2, 2, params)
  Kstar <- rbind(
    cbind(Kstar11, Kstar12),
    cbind(Kstar21, Kstar22)
  )

  # Test prior covariance of [f1*; f2*] (noise-free)
  Kstarstar11 <- Kmat_ij_D2_Q1(ts1, ts1, 1, 1, params)
  Kstarstar12 <- Kmat_ij_D2_Q1(ts1, ts2, 1, 2, params)
  Kstarstar22 <- Kmat_ij_D2_Q1(ts2, ts2, 2, 2, params)
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
    d = factor(rep(c(1, 2), times = c(S1, S2)))
  )
}

# Draw multiple posterior predictive samples
pred_D2_Q1 <- function(
  draws, # df
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
    r <- as.data.frame( draws[idx[d], ] )

    pr <- pred_D2_Q1_draw(r, t1, t2, y1, y2, y1_se, y2_se, ts1, ts2) |>
      mutate(.draw = d)

    out_list[[d]] <- pr
  }

  bind_rows(out_list, .id = ".draw")
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

  params <- list(
    w = w,
    mu = mu,
    Sigma = Sigma,
    theta = theta,
    phi = phi
  )

  K <- K_D2_Q1(t1, t2, params) +
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

  params_df <- data.frame(
    seed = seed,
    w = w,
    mu = mu,
    Sigma = Sigma,
    theta = theta,
    phi = phi
  )

  return(list(params_df, output_df))
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
