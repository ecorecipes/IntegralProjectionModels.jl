## R cross-check for 06_stochastic.qmd (Stochastic IPM)
## Source: ipmr/tests/testthat/test-simple_di_stoch_kern.R (simplified)

# --- Fixed Parameters ---
s_int <- 1.03;    s_slope <- 2.2
g_int <- 8.0;     g_slope <- 0.92;  sd_g <- 0.9
f_r_int <- 0.09;  f_r_slope <- 0.05
f_s_int <- 0.1;   f_s_slope <- 0.005
mu_fd <- 9.0;     sd_fd <- 2.0

# --- Random intercepts ---
set.seed(50127)
g_r <- rnorm(5, 0, 0.3)
s_r <- rnorm(5, 0, 0.7)
f_s_r <- rnorm(5, 0, 0.2)

# --- Domain ---
L <- 0.2; U <- 40; n <- 100
b <- seq(L, U, length.out = n + 1)
sv1 <- (b[2:(n+1)] + b[1:n]) * 0.5
h <- sv1[2] - sv1[1]

# --- Vital rate functions ---
f_r <- function(sv1, params) {
  1/(1 + exp(-(params[1] + params[2] * sv1)))
}

s <- function(sv1, params, r_effect) {
  1/(1 + exp(-(params[1] + params[2] * sv1 + r_effect))) *
    (1 - f_r(sv1, c(f_r_int, f_r_slope)))
}

g <- function(sv1, sv2, params, r_effect, L, U) {
  mu <- params[1] + params[2] * sv1 + r_effect
  ev <- pnorm(U, mu, params[3]) - pnorm(L, mu, params[3])
  dnorm(sv2, mean = mu, sd = params[3]) / ev
}

f_d <- function(sv2, params, L, U) {
  ev <- pnorm(U, params[1], params[2]) - pnorm(L, params[1], params[2])
  dnorm(sv2, mean = params[1], sd = params[2]) / ev
}

fec <- function(sv1, sv2, params, r_effect, L, U) {
  f_r(sv1, params[1:2]) * exp(params[3] + params[4] * sv1 + r_effect) *
    f_d(sv2, params[5:6], L, U)
}

# --- Build 5 year-specific kernels ---
domains <- expand.grid(d1 = sv1, d2 = sv1)
K_list <- list()

for (i in 1:5) {
  G <- h * g(domains$d2, domains$d1,
             params = c(g_int, g_slope, sd_g),
             r_effect = g_r[i], L = L, U = U)
  s_vals <- s(sv1, c(s_int, s_slope), s_r[i])
  P_i <- t(s_vals * t(matrix(G, nrow = n, ncol = n, byrow = TRUE)))
  F_i <- h * fec(domains$d2, domains$d1,
                 params = c(f_r_int, f_r_slope, f_s_int, f_s_slope, mu_fd, sd_fd),
                 r_effect = f_s_r[i], L = L, U = U)
  F_i <- matrix(F_i, nrow = n, ncol = n, byrow = TRUE)
  K_list[[i]] <- P_i + F_i
}

# --- Eigenvalues ---
lambdas <- sapply(K_list, function(K) Re(eigen(K)$values[1]))

# --- Mean kernel ---
K_mean <- Reduce("+", K_list) / 5
lam_det <- Re(eigen(K_mean)$values[1])

cat("=== Stochastic IPM Cross-Check (R) ===\n")
cat("Year-specific lambdas:", lambdas, "\n")
cat("Deterministic lambda (mean kernel):", lam_det, "\n")
cat("\nNote: stochastic lambda requires simulation and is not computed here.\n")
cat("Tuljapurkar's inequality: lambda_s <= lambda_det\n")
