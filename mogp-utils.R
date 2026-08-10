library(mvtnorm)
library(tictoc)

# Returns the cross-covariance between scalar location x in band i and
# scalar location x_prime in band j, conditioned on hyperparameter values.
k_ij <- function(x, x_prime,
                 w_ij, Sigma_ij, mu_ij, theta_ij, phi_ij) {

  tau <- abs(x - x_prime)
  alpha_ij <- w_ij * sqrt( 2.0*pi * abs(Sigma_ij) )

  return(
    alpha_ij *
      exp(-0.5 * (tau + theta_ij)^2 * Sigma_ij) *
      cos(2*pi*(tau + theta_ij) * mu_ij + phi_ij) # NB: 2*pi factor
  )
}

# Returns the autocovariance between scalar location x and scalar location
# x_prime within band i.
k_ii <- function(x, x_prime,
                 w_i, Sigma_i, mu_i) {

  tau <- abs(x - x_prime)
  alpha_ii <- w_i^2.0 * sqrt( 2.0*pi * abs(Sigma_i) )

  return(
    alpha_ii *
      exp(-0.5*tau^2 * Sigma_i) *
      cos(2*pi*tau*mu_i) # NB: 2*pi factor
  )
}

# Returns a [n x n_prime] matrix containing the covariance between the elements
# of the vector of locations x in band i with the elements of the vector of
# locations x_prime in band j.
K_ij <- function(xs, xs_prime, i, j,
                 ws, Sigmas, mus, thetas, phis) {

  n <- length(xs)
  n_prime <- length(xs_prime)

  Kij <- matrix(NA, nrow = n, ncol = n_prime)

  for (r in 1:n) {
    for (c in 1:n_prime) {
      if (i == j) {
        Kij[r,c] <- k_ii(xs[r], xs_prime[c], ws[i], Sigmas[i], mus[i])
      }
      else {
        #w_ij <- ws[i]*ws[j] * exp( -1.0/4*(mus[i]-mus[j]) * (Sigmas[i]+Sigmas[j])^-1 * (mus[i]-mus[j]) )
        w_ij <- ws[i]*ws[j] * exp( -pi^2*(mus[i]-mus[j]) * (Sigmas[i]+Sigmas[j])^-1 * (mus[i]-mus[j]) )
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

# Returns a square cross-covariance matrix of elements of the vector of
# locations x with itself in the same band.
Kxx_mat <- function(x, D, n,
                    ws, Sigmas, mus, thetas, phis) {

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
      Kxxmat[jj,ii] <- Kxxmat[ii,jj] # enforce symmetry
    }
  }

  return(Kxxmat)
}

# Returns a [n x n_star] cross-covariance matrix between the D bands comprising
# a D x D arrangement of covariance matrices of elements of the vector of
# locations xs in band i with the elements of the vector of locations xs_star
# in band j.
K_mat <- function(xs, xs_star, ds, ds_star, D,
                  ws, Sigmas, mus, thetas, phis) {

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

# Simulate a single draw from an MOSK GP (Q = 1)
simulate_Q1_moskgp <- function(D,
                               ns = 100,
                               noise_sigma = 0.25,
                               masked_pct = 0.2,
                               zero_phi = FALSE,
                               seed) {

  xs <- rep( seq(from = 0, to = 1, length.out = ns), times = D)
  ds <- rep(1:D, each = ns)
  N <- length(xs)

  set.seed(seed)

  ws <-     round( abs( rnorm(D, mean = 0, sd = 1.0) ), 2 )
  Sigmas <- round( abs( rnorm(D, mean = 0, sd = 1.0) ), 2 )
  mus <-    round( abs( rnorm(D, mean = 5, sd = 1.0) ), 2 )
  thetas <- sort( round( rnorm(D, mean = 0, sd = 1.0), 2 ) )

  if (zero_phi) {
    phis <- rep(0.0, D)
  } else {
    phis <- round( rnorm(D, mean = 0, sd = 1.0), 2 )
  }

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

  params_df = data.frame(
    seed = seed,
    w = ws,
    Sigma = Sigmas,
    mu = mus,
    theta = thetas,
    phi = phis
  )

  return(list(params_df, output_df))
}

# Generate a single variate from a MOSK GP conditioned on hyperparameter values
postpred_Q1_draw <- function(
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
  K_star = K_mat(x, x_star, d, d_star, D, w, Sigma, mu, theta, phi)

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

# Generate a n_draw variates conditioned on hyperparameter values
postpred_Q1_draws <- function(n_draws = 1,
                              x, y, y_se, d, D,
                              ns, # no. observations by band
                              x_star, d_star, ns_star,
                              w, Sigma, mu, theta, phi,
                              seed = NULL,
                              epsilon = 1e-9) {

  l <- vector("list", n_draws)

  for (i in 1:n_draws) {

    one_draw_df <- postpred_Q1_draw(
      x = x, y = y, y_se = y_se,
      d = d, D = D, ns = ns, x_star = x_star,
      d_star = d_star, ns_star = ns_star,
      w = w, Sigma = Sigma, mu = mu, theta = theta,
      phi = phi, seed = seed, epsilon = epsilon)

    l[[i]] <- one_draw_df
  }

  return( bind_rows(l, .id = ".draw"))
}

# Check which posterior predictive draws generate a proper covariance matrix.
check_valid_postpred_draws <- function(
    x, y, y_se, d, D, ns,
    x_star, d_star, ns_star,
    ws, Sigmas, mus, thetas, phis,
    seed = NULL,
    epsilon = 1e-9) {

  n_draws <- NROW(ws)
  pp_list <- vector("list", n_draws)
  valid_pps <- rep(FALSE, n_draws)

  warning_count <- 0
  error_count <- 0
  valid_count <- 0


  for (r in 1:n_draws) {
    if (r == 1) {
      format(start_time <- Sys.time())
    }

    this_w <- ws[r,]
    this_Sigma <- Sigmas[r,]
    this_mu <- mus[r,]
    this_theta <- thetas[r,]
    this_phi <- phis[r,]

    return_string <- tryCatch(
      warning = function(cnd) {
        pp_list[[r]] <- NA
        return("warning")
      },
      error = function(cnd) {
        pp_list[[r]] <- NA
        return("error")
      },
      {
        new_draw <- postpred_Q1_draw(
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

    if (return_string == "warning") {
      warning_count <- warning_count + 1
    } else if (return_string == "error") {
      error_count <- error_count + 1
    } else {
      valid_count <- valid_count + 1
    }

    if (r %% 100 == 0) {

      cat(
        paste0(
          r,"/", n_draws,
          ":\tWarnings = ", warning_count,
          ",\tErrors = ", error_count,
          ",\tValid = ", valid_count,
          "\t(", format( difftime(Sys.time(), start_time), digits = 3 ), ")",
          "\n")
      )
    }

  }

  pp_df <- bind_rows(pp_list, .id = ".draw")

  cat(
    paste0(
      "Total draws: ", n_draws,
      "\nTotal warnings: ", warning_count,
      "\nTotal errors: ", error_count, "\n")
  )

  return( list(pp_df, valid_pps) )
}



check_valid_postpred_draws_parallel <- function(
    x, y, y_se, d, D, ns,
    x_star, d_star, ns_star,
    ws, Sigmas, mus, thetas, phis,
    seed = NULL,
    epsilon = 1e-9) {

  n_draws <- NROW(ws)
  pp_list <- vector("list", n_draws)
  valid_pps <- rep(FALSE, n_draws)

  warning_count <- 0
  error_count <- 0
  valid_count <- 0


  for (r in 1:n_draws) {
    if (r == 1) {
      format(start_time <- Sys.time())
    }

    this_w <- ws[r,]
    this_Sigma <- Sigmas[r,]
    this_mu <- mus[r,]
    this_theta <- thetas[r,]
    this_phi <- phis[r,]

    return_string <- tryCatch(
      warning = function(cnd) {
        pp_list[[r]] <- NA
        return("warning")
      },
      error = function(cnd) {
        pp_list[[r]] <- NA
        return("error")
      },
      {
        new_draw <- postpred_Q1_draw(
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

    if (return_string == "warning") {
      warning_count <- warning_count + 1
    } else if (return_string == "error") {
      error_count <- error_count + 1
    } else {
      valid_count <- valid_count + 1
    }

    if (r %% 100 == 0) {

      cat(
        paste0(
          r,"/", n_draws,
          ":\tWarnings = ", warning_count,
          ",\tErrors = ", error_count,
          ",\tValid = ", valid_count,
          "\t(", format( difftime(Sys.time(), start_time), digits = 3 ), ")",
          "\n")
      )
    }

  }

  pp_df <- bind_rows(pp_list, .id = ".draw")

  cat(
    paste0(
      "Total draws: ", n_draws,
      "\nTotal warnings: ", warning_count,
      "\nTotal errors: ", error_count, "\n")
  )

  return( list(pp_df, valid_pps) )
}

