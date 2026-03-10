## R cross-check for 03_eviction.qmd (Eviction correction)
## Source: first-edition/Rcode/c2/MonocarpGrowthEviction.R

# --- Parameters ---
m.par.true <- c(
  surv.int = -0.65, surv.z = 0.75,
  flow.int = -18.0, flow.z = 6.9,
  grow.int = 0.96, grow.z = 0.59, grow.sd = 0.67,
  rcsz.int = -0.08, rcsz.sd = 0.76,
  seed.int = 1.0, seed.z = 2.2,
  p.r = 0.007
)

# --- Vital rates ---
s_z <- function(z, m.par) 1/(1 + exp(-(m.par["surv.int"] + m.par["surv.z"] * z)))
p_bz <- function(z, m.par) 1/(1 + exp(-(m.par["flow.int"] + m.par["flow.z"] * z)))
G_z1z <- function(z1, z, m.par) dnorm(z1, m.par["grow.int"] + m.par["grow.z"] * z, m.par["grow.sd"])
b_z <- function(z, m.par) exp(m.par["seed.int"] + m.par["seed.z"] * z)
c_0z1 <- function(z1, m.par) dnorm(z1, m.par["rcsz.int"], m.par["rcsz.sd"])

P_z1z <- function(z1, z, m.par) (1 - p_bz(z, m.par)) * s_z(z, m.par) * G_z1z(z1, z, m.par)
F_z1z <- function(z1, z, m.par) p_bz(z, m.par) * b_z(z, m.par) * m.par["p.r"] * c_0z1(z1, m.par)

# --- No correction ---
mk_K <- function(m, m.par, L, U) {
  h <- (U - L)/m
  meshpts <- L + ((1:m) - 1/2) * h
  P <- h * outer(meshpts, meshpts, P_z1z, m.par = m.par)
  F <- h * outer(meshpts, meshpts, F_z1z, m.par = m.par)
  K <- P + F
  list(K = K, meshpts = meshpts, P = P, F = F, h = h)
}

# --- Truncated distributions ---
mk_K_trunc <- function(m, m.par, L, U) {
  h <- (U - L)/m
  meshpts <- L + ((1:m) - 1/2) * h

  P <- matrix(0, m, m)
  for (j in 1:m) {
    mu <- m.par["grow.int"] + m.par["grow.z"] * meshpts[j]
    sd <- m.par["grow.sd"]
    ev <- pnorm(U, mu, sd) - pnorm(L, mu, sd)
    surv_eff <- s_z(meshpts[j], m.par) * (1 - p_bz(meshpts[j], m.par))
    for (i in 1:m) {
      P[i, j] <- h * surv_eff * dnorm(meshpts[i], mu, sd) / ev
    }
  }
  F <- h * outer(meshpts, meshpts, F_z1z, m.par = m.par)
  K <- P + F
  list(K = K, meshpts = meshpts, P = P, F = F, h = h)
}

# --- Compute ---
L <- -2.65; U <- 4.5; m <- 250

IPM.nocorr <- mk_K(m, m.par.true, L, U)
IPM.trunc <- mk_K_trunc(m, m.par.true, L, U)

# Column sum diagnostics
meshpts <- IPM.nocorr$meshpts
true_surv <- s_z(meshpts, m.par.true) * (1 - p_bz(meshpts, m.par.true))
col_sums_nocorr <- colSums(IPM.nocorr$P)
col_sums_trunc <- colSums(IPM.trunc$P)

lam_nocorr <- Re(eigen(IPM.nocorr$K)$values[1])
lam_trunc <- Re(eigen(IPM.trunc$K)$values[1])

cat("=== Eviction Correction Cross-Check (R) ===\n")
cat("Lambda (no correction):", lam_nocorr, "\n")
cat("Lambda (truncated):    ", lam_trunc, "\n")
cat("Difference:            ", lam_trunc - lam_nocorr, "\n")
cat("\nMax column sum error (no correction):", max(abs(col_sums_nocorr - true_surv)), "\n")
cat("Max column sum error (truncated):    ", max(abs(col_sums_trunc - true_surv)), "\n")
