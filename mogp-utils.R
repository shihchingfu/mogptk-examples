library(mvtnorm)

# Cross-covariance between x in band i and x' in band j
k_ij <- function(x, x_prime, w_ij, Sigma_ij, mu_ij, theta_ij, phi_ij) {
  alpha_ij <- w_ij * (2.0*pi)^(1/2.0) * sqrt(abs(Sigma_ij))
  tau <- abs(x - x_prime)

  return( alpha_ij *
            exp(-0.5 * (tau + theta_ij)^2 * Sigma_ij) *
            cos((tau + theta_ij) * mu_ij + phi_ij) )
}

# Autocovariance between x and x' within band i
k_ii <- function(x, x_prime, w_i, Sigma_i, mu_i) {
  alpha_ii <- w_i^2.0 * (2.0*pi)^(1/2.0) * sqrt(abs(Sigma_i))
  tau <- abs(x - x_prime)

  return( alpha_ii * exp(-0.5*tau^2 * Sigma_i) * cos(tau*mu_i) )
}

# the n x n' covariance matrix between elements _x_ in band i and elements _x'_ in band j
K_ij <- function(xs, xs_prime, i, j, ws, Sigmas, mus, thetas, phis) {

  n <- length(xs)
  n_prime <- length(xs_prime)

  Kij <- matrix(NA, nrow = n, ncol = n_prime)

  for (r in 1:n) {
    for (c in 1:n_prime) {
      if (i == j) {
        Kij[r,c] <- k_ii(xs[r], xs_prime[c], ws[i], Sigmas[i], mus[i])
      }
      else {
        w_ij <- ws[i]*ws[j] * exp( -1.0/4*(mus[i]-mus[j]) * (Sigmas[i]+Sigmas[j])^-1 * (mus[i]-mus[j]) )
        Sigma_ij <- 2*Sigmas[i] * (Sigmas[i] + Sigmas[j])^-1 * Sigmas[j]
        mu_ij <- (Sigmas[i] + Sigmas[j])^-1 * (Sigmas[i]*mus[j] + Sigmas[j]*mus[i])
        theta_ij <- thetas[i] - thetas[j]
        phi_ij <- phis[i] - phis[j]

        Kij[r,c] = k_ij(xs[r], xs_prime[c], w_ij, Sigma_ij, mu_ij, theta_ij, phi_ij)
      }
    }
  }
  return(Kij)
}

# Square cross-covariance matrix of locations x with itself
Kxx_mat <- function(x, D, n, ws, Sigmas, mus, thetas, phis) {

  end_idx <- cumsum(n)
  start_idx <- 1 + end_idx - n

  N <- length(x)

  Kxxmat <- matrix(NA, nrow = N, ncol = N)

  for (r in 1:D) {
    for (c in 1:D) {

      nrows <- n[r]
      rcols <- n[c]

      start_row <- start_idx[r];
      end_row <- end_idx[r];
      start_col <- start_idx[c];
      end_col <- end_idx[c];

      K_ij_mat <- K_ij(x[start_row:end_row],
                       x[start_col:end_col],
                       r, c, ws, Sigmas, mus, thetas, phis)

      Kxxmat[start_row:end_row, start_col:end_col] <- K_ij_mat
    }
  }

  for (ii in 1:N) {
    for (jj in 1:N) {
      Kxxmat[jj,ii] <- Kxxmat[ii,jj]
    }
  }

  return(Kxxmat)
}

# the cross-covariance matrix between the D bands comprising a D x D arrangement
# of covariance matrices between i and j bands.
K_mat_v1 <- function(xs, xs_star, ds, ds_star, D, ws, Sigmas, mus, thetas, phis) {

  N_rows <- length(xs)
  N_cols <- length(xs_star)
  Kmat <- matrix(NA, nrow = N_rows, ncol = N_cols)

  for (i in 1:D) { # for each combination of bands
    for (j in 1:D) {

      i_indices <- which(ds == i)
      j_indices <- which(ds_star == j)

      x_i <- xs[i_indices]
      x_j <- xs_star[j_indices]

      Kij_mat <- K_ij(x_i, x_j, i, j, ws, Sigmas, mus, thetas, phis)

      Kmat[i_indices, j_indices] <- Kij_mat
    }
  }

  return(Kmat)
}


# the cross-covariance matrix between the D bands comprising a D x D arrangement of covariance matrices between i and j bands.
K_mat_v2 <- function(xs, xs_star, ds, ds_star, D, ws, Sigmas, mus, thetas, phis) {

  N <- length(xs)
  N_star <- length(xs_star)

  Kmat <- matrix(NA, nrow = N, ncol = N_star)

  mgrid <- meshgrid(ds, ds_star)

  r_grid <- mgrid$Y
  c_grid <- mgrid$X

  for (r in 1:D) {
    for (c in 1:D) {

      r_indices <- which(ds == r)
      c_indices <- which(ds_star == c)

      xs <- xs[r_indices]
      xs_prime <- xs_star[c_indices]

      K_rc <- K_ij(xs, xs_prime, r, c, ws, Sigmas, mus, thetas, phis)

      mask <- r_grid == r & c_grid == c

      Kmat[mask] <- K_rc
    }
  }

  return(Kmat)
}

# Simulate a draw from an MOGP
simulate_mogp <- function(D, Q, ns = 100, noise_sigma = 0.25, masked_pct = 0.2, seed) {

  xs <- rep( seq(from = 0, to = 1, length.out = ns), times = D)
  ds <- rep(1:D, each = ns)
  N <- length(xs)

  set.seed(seed)

  ws <-     round( abs( rnorm(D, mean = 0, sd = 5.0) ), 2 )
  Sigmas <- round( abs( rnorm(D, mean = 0, sd = 5.0) ), 2 )
  mus <-    round( sort( rnorm(D, mean = 0, sd = 5.0) ), 2 )
  thetas <- round( rnorm(D, mean = 0, sd = 5.0), 2 )
  phis <-   rep(0.0, D)

  KK <- Kxx_mat(xs, D, rep(ns, D), ws, Sigmas, mus, thetas, phis)

  set.seed(seed)
  output_df <- data.frame(
    t = xs,
    f = c(rmvnorm(1, sigma = KK)),
    y = rep(NA, N),
    y_se = rep(noise_sigma, N),
    d = factor(ds),
    masked = runif(N) < masked_pct
  ) |>
    mutate(y = f + rnorm(length(xs), mean = 0, sd = noise_sigma))

  cat("Seed =", seed)
  cat("\nws =", ws)
  cat("\nSigmas =", Sigmas)
  cat("\nmus =", mus)
  cat("\nthetas =", thetas)
  cat("\nphis =", phis)

  return(output_df)
}

# Generate a single variate
postpred_draw <- function(
    x, y, y_se, d,
    D,
    ns, # no. observations by band
    x_star,
    d_star,
    ns_star,
    w, Sigma, mu, theta, phi,
    seed = NULL,
    epsilon = 1e-9) {

  N <- length(x);
  N_star <- length(x_star);

  # N x N covariance KS = K(X,X) + Sigma_noise
  K_xx <- Kxx_mat(x, D, ns, w, Sigma, mu, theta, phi)
  diag(K_xx) <- diag(K_xx) + y_se^2

  # N x N* covariance K* = K(X,X*)
  K_star = K_mat_v1(x, x_star, d, d_star, D, w, Sigma, mu, theta, phi)

  # N* x N* covariance K** = K(X*,X*)
  K_starstar <- Kxx_mat(x_star, D, ns_star, w, Sigma, mu, theta, phi)

  fstar_mu <- t(K_star) %*% solve(K_xx) %*% y
  fstar_Sigma <- K_starstar - t(K_star) %*% solve(K_xx) %*% K_star

  diag(fstar_Sigma) <- diag(fstar_Sigma) + epsilon

  set.seed(seed)
  f_star <- rmvnorm(
    n = 1,
    mean = fstar_mu,
    sigma = fstar_Sigma)

  result_df <- data.frame(
    d_star = factor(d_star),
    x_star,
    f_star = t(f_star)
  )

  return(result_df)
}

postpred_draws <- function(n_draws = 1,
                           x,
                           y,
                           y_se,
                           d,
                           D,
                           ns, # no. observations by band
                           x_star,
                           d_star,
                           ns_star,
                           w, Sigma, mu, theta, phi,
                           seed = NULL,
                           epsilon = 1e-9) {

  l <- vector("list", n_draws)

  for (i in 1:n_draws) {

    one_draw_df <- postpred_draw(
      x = x, y = y, y_se = y_se,
      d = d, D = D, ns = ns, x_star = x_star,
      d_star = d_star, ns_star = ns_star,
      w = w, Sigma = Sigma, mu = mu, theta = theta,
      phi = phi, seed = seed, epsilon = epsilon)

    l[[i]] <- one_draw_df
  }

  return( bind_rows(l, .id = ".draw"))
}

postpred_valid_draws <- function(
    x,
    y,
    y_se,
    d,
    D,
    ns,
    x_star,
    d_star,
    ns_star,
    ws,
    Sigmas,
    mus,
    thetas,
    phis,
    seed = NULL,
    epsilon = 1e-9) {

  n_draws <- NROW(ws)
  pp_list <- vector("list", n_draws)
  valid_pps <- rep(FALSE, n_draws)

  for (r in 1:n_draws) {
    if (r %% 10 == 0)
      cat(r,"of", n_draws,"\n")

    this_w <- ws[r,]
    this_Sigma <- Sigmas[r,]
    this_mu <- mus[r,]
    this_theta <- thetas[r,]
    this_phi <- phis[r,]

    tryCatch(
      warning = function(cnd) {
        cat("Warning:", r, conditionMessage(cnd), "\n")
        pp_list[[r]] <- NA
      },
      error = function(cnd) {
        cat("Error:", r, conditionMessage(cnd), "\n")
        pp_list[[r]] <- NA
      },
      {
        new_draw <- postpred_draw(
          x = x,
          y = y,
          y_se = y_se,
          d = d,
          D = D,
          ns = ns,
          x_star = x_star,
          d_star = d_star,
          ns_star = ns_star,
          this_w,
          this_Sigma,
          this_mu,
          this_theta,
          this_phi
        ) |>
          mutate(.draw = r)

        pp_list[[r]] <- new_draw
        valid_pps[r] <- TRUE
      }
    )
  }

  pp_df <- bind_rows(pp_list, .id = ".draw")

  return( list(pp_df, valid_pps) )

}
