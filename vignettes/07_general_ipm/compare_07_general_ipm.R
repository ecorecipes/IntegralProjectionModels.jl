## R cross-check for 07_general_ipm.qmd (Ligustrum with seed bank)
## Source: ipmr/tests/testthat/test-general_di_det.R (exact reproduction)

# --- Control parameters ---
data_list_control <- list(
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

# --- CR parameters ---
data_list_cr <- list(
  g_int     = 7.229,
  g_slope   = 0.988,
  g_sd      = 21.72262,
  s_int     = 0.0209,
  s_slope   = 0.0831,
  s_slope_2 = -0.00012999,
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

# --- Functions ---
s_x <- function(int, slope1, slope2, sv1) {
  1/(1 + exp(-(int + slope1 * sv1 + slope2 * sv1^2)))
}

f_r_x <- function(int, slope1, sv1) {
  1/(1 + exp(-(int + slope1 * sv1)))
}

g_x <- function(sv2, sv1, int, slope, gsd, L, U) {
  mu <- int + slope * sv1
  ev <- pnorm(U, mu, gsd) - pnorm(L, mu, gsd)
  dnorm(sv2, mean = mu, sd = gsd) / ev
}

# --- Build mega-kernel (following exact test code) ---
build_mega <- function(dl) {
  b <- seq(L, U, length.out = n + 1)
  d1 <- (b[2:(n+1)] + b[1:n]) * 0.5
  d2 <- d1
  h <- d1[3] - d1[2]

  # expand.grid with d2 varying fastest
  domains <- expand.grid(list(d2 = d1, d1 = d1))

  # Growth: g_x(sv2=domains$d1, sv1=domains$d2, ...)
  # Note argument swap: sv2 gets d1 (slow-varying), sv1 gets d2 (fast-varying)
  G <- h * g_x(domains$d1, domains$d2,
               int = dl$g_int, slope = dl$g_slope,
               gsd = dl$g_sd, L = L, U = U)

  # Survival: function of d1 (500 values)
  s <- s_x(dl$s_int, dl$s_slope, dl$s_slope_2, d1)

  # P = s * G  (s recycled across 250000 elements)
  P <- s * G

  # CC block: fill by row
  cc <- matrix(P, nrow = n, ncol = n, byrow = TRUE)

  # CD: continuous -> discrete (1 x n row vector)
  cd <- f_r_x(dl$f_r_int, dl$f_r_slope, d1) *
    exp(dl$f_s_int + dl$f_s_slope * d1) *
    dl$g_i
  cd <- matrix(cd, nrow = 1, ncol = n)

  # DC: discrete -> continuous (n x 1 column vector)
  dc <- (dnorm(d2, mean = dl$f_d_mu, sd = dl$f_d_sd) /
         (pnorm(U, dl$f_d_mu, dl$f_d_sd) - pnorm(L, dl$f_d_mu, dl$f_d_sd))) *
    dl$e_p * h
  dc <- matrix(dc, nrow = n, ncol = 1)

  # DD: 0
  dd <- 0

  # Assemble: discrete state first (row 1), then continuous (rows 2:501)
  K <- rbind(
    cbind(dd, cd),
    cbind(dc, cc)
  )

  list(K = K, meshpts = d1, h = h)
}

# --- Compute ---
ipm_ctrl <- build_mega(data_list_control)
ipm_cr <- build_mega(data_list_cr)

lam_ctrl <- Re(eigen(ipm_ctrl$K)$values[1])
lam_cr <- Re(eigen(ipm_cr$K)$values[1])

cat("=== General IPM (Ligustrum) Cross-Check (R) ===\n")
cat("Control lambda:", lam_ctrl, "\n")
cat("CR lambda:     ", lam_cr, "\n")
cat("Published control: ~1.22, Published CR: ~1.26\n")
