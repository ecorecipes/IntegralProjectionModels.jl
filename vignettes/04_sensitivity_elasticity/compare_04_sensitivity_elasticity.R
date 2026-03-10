## R cross-check for 04_sensitivity_elasticity.qmd
## Source: first-edition/Rcode/c2/Ungulate model

# --- Parameters ---
m.par.true <- c(
  surv.int = -9.65, surv.z = 3.77,
  grow.int = 1.41, grow.z = 0.557, grow.sd = 0.0799,
  repr.int = -7.23, repr.z = 2.60,
  recr.int = 1.93,
  rcsz.int = 0.362, rcsz.z = 0.709, rcsz.sd = 0.159
)

# --- Vital rates ---
s_z <- function(z, m.par) 1/(1 + exp(-(m.par["surv.int"] + m.par["surv.z"] * z)))
g_z1z <- function(z1, z, m.par) dnorm(z1, m.par["grow.int"] + m.par["grow.z"] * z, m.par["grow.sd"])
pb_z <- function(z, m.par) 1/(1 + exp(-(m.par["repr.int"] + m.par["repr.z"] * z)))
pr_z <- function(m.par) 1/(1 + exp(-m.par["recr.int"]))
c_z1z <- function(z1, z, m.par) dnorm(z1, m.par["rcsz.int"] + m.par["rcsz.z"] * z, m.par["rcsz.sd"])

P_z1z <- function(z1, z, m.par) s_z(z, m.par) * g_z1z(z1, z, m.par)
F_z1z <- function(z1, z, m.par) s_z(z, m.par) * pb_z(z, m.par) * 0.5 * pr_z(m.par) * c_z1z(z1, z, m.par)

mk_K <- function(m, m.par, L, U) {
  h <- (U - L)/m
  meshpts <- L + ((1:m) - 1/2) * h
  P <- h * outer(meshpts, meshpts, P_z1z, m.par = m.par)
  F <- h * outer(meshpts, meshpts, F_z1z, m.par = m.par)
  K <- P + F
  list(K = K, meshpts = meshpts, P = P, F = F, h = h)
}

# --- Compute ---
IPM <- mk_K(100, m.par.true, 1.6, 3.7)

eig <- eigen(IPM$K)
lam <- Re(eig$values[1])
w <- Re(eig$vectors[, 1]); w <- w / sum(w)

eig.left <- eigen(t(IPM$K))
v <- Re(eig.left$vectors[, 1]); v <- v / v[1]

# Sensitivity
vw <- sum(v * w)
S <- outer(v, w) / vw

# Elasticity
E <- (IPM$K * S) / lam
E_P <- (IPM$P * S) / lam
E_F <- (IPM$F * S) / lam

cat("=== Sensitivity/Elasticity Cross-Check (R) ===\n")
cat("Lambda:", lam, "\n")
cat("Sum of elasticities:", sum(E), "\n")
cat("Elasticity of P:", sum(E_P), "\n")
cat("Elasticity of F:", sum(E_F), "\n")
