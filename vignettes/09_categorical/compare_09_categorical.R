## R cross-check for 09_categorical.qmd (categorical composition)

# --- Control parameters ---
ctrl <- list(
  g_int     = 5.781,
  g_slope   = 0.988,
  g_sd      = 20.55699,
  s_int     = -0.352,
  s_slope   = 0.122,
  s_slope_2 = -0.000213,
  f_r_int   = -11.46,
  f_r_slope = 0.0835,
  f_s_int   = 2.6204,
  f_s_slope = 0.01256,
  f_d_mu    = 5.6655,
  f_d_sd    = 2.0734,
  e_p       = 0.15,
  g_i       = 0.5067
)

# --- Domain ---
L <- 1.02; U <- 624; n <- 500

# --- Vital rates ---
survival <- function(z, p) {
  1 / (1 + exp(-(p$s_int + p$s_slope * z + p$s_slope_2 * z^2)))
}

growth <- function(z_prime, z, p) {
  mu <- p$g_int + p$g_slope * z
  ev <- pnorm(U, mu, p$g_sd) - pnorm(L, mu, p$g_sd)
  ifelse(ev > 0, dnorm(z_prime, mu, p$g_sd) / ev, dnorm(z_prime, mu, p$g_sd))
}

P_fn <- function(z_new, z) survival(z, ctrl) * growth(z_new, z, ctrl)

# --- Left Kan extension ---
left_kan_extension <- function(kernel_fn, lower, upper, n_bins) {
  h <- (upper - lower) / n_bins
  midpoints <- lower + ((1:n_bins) - 0.5) * h
  A <- matrix(0, n_bins, n_bins)
  for (j in 1:n_bins) {
    for (i in 1:n_bins) {
      A[i, j] <- h * kernel_fn(midpoints[i], midpoints[j])
    }
  }
  list(A = A, midpoints = midpoints, h = h)
}

# --- Stratify ---
stratify <- function(A_local, dispersal) {
  n_bins <- nrow(A_local)
  n_patches <- nrow(dispersal)
  n_total <- n_bins * n_patches
  A_strat <- matrix(0, n_total, n_total)
  for (p_to in 1:n_patches) {
    for (p_from in 1:n_patches) {
      rows <- ((p_to - 1) * n_bins + 1):(p_to * n_bins)
      cols <- ((p_from - 1) * n_bins + 1):(p_from * n_bins)
      A_strat[rows, cols] <- dispersal[p_to, p_from] * A_local
    }
  }
  A_strat
}

# --- Coarsen ---
coarsen <- function(A, n_coarse) {
  n_fine <- nrow(A)
  bins_per_coarse <- n_fine %/% n_coarse
  coarsening_map <- ((1:n_fine) - 1) %/% bins_per_coarse + 1
  A_coarse <- matrix(0, n_coarse, n_coarse)
  for (j_fine in 1:n_fine) {
    c_j <- coarsening_map[j_fine]
    for (i_fine in 1:n_fine) {
      c_i <- coarsening_map[i_fine]
      A_coarse[c_i, c_j] <- A_coarse[c_i, c_j] + A[i_fine, j_fine] / bins_per_coarse
    }
  }
  A_coarse
}

# === Tests ===
cat("=== Categorical Composition Cross-Check (R) ===\n\n")

# 1. Left Kan extension convergence
ref <- left_kan_extension(P_fn, L, U, 500)
lambda_ref <- Re(eigen(ref$A)$values[1])
cat("Reference CC lambda (n=500):", lambda_ref, "\n\n")

cat("Convergence:\n")
for (nb in c(5, 10, 20, 50, 100, 200, 500)) {
  res <- left_kan_extension(P_fn, L, U, nb)
  lam <- Re(eigen(res$A)$values[1])
  cat(sprintf("  n=%3d: lambda=%.6f, delta=%.2e\n", nb, lam, abs(lam - lambda_ref)))
}

# 2. Stratification
D_sym <- matrix(c(0.9, 0.1, 0.1, 0.9), 2, 2)
K_strat <- stratify(ref$A, D_sym)
lambda_strat <- Re(eigen(K_strat)$values[1])
cat(sprintf("\nSymmetric stratification: lambda=%.6f (ratio=%.6f)\n",
            lambda_strat, lambda_strat / lambda_ref))

# 3. Coarsening
res_100 <- left_kan_extension(P_fn, L, U, 100)
lambda_100 <- Re(eigen(res_100$A)$values[1])
A_coarsened <- coarsen(res_100$A, 20)
lambda_coarsened <- Re(eigen(A_coarsened)$values[1])
cat(sprintf("\nCoarsening: fine(100)=%.6f, coarsened(20)=%.6f\n",
            lambda_100, lambda_coarsened))
