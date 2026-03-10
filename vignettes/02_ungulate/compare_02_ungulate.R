## R cross-check for 02_ungulate.qmd (Soay sheep DI deterministic IPM)
## Source: first-edition/Rcode/c2/Ungulate Demog Funs.R + Ungulate Calculations.R

# --- Parameters ---
m.par.true <- c(
  surv.int  = -9.65,
  surv.z    =  3.77,
  grow.int  =  1.41,
  grow.z    =  5.57e-1,
  grow.sd   =  7.99e-2,
  repr.int  = -7.23,
  repr.z    =  2.60,
  recr.int  =  1.93,
  rcsz.int  =  3.62e-1,
  rcsz.z    =  7.09e-1,
  rcsz.sd   =  1.59e-1
)

# --- Vital rate functions ---
s_z <- function(z, m.par) {
  linear.p <- m.par["surv.int"] + m.par["surv.z"] * z
  1/(1 + exp(-linear.p))
}

g_z1z <- function(z1, z, m.par) {
  mu <- m.par["grow.int"] + m.par["grow.z"] * z
  dnorm(z1, mean = mu, sd = m.par["grow.sd"])
}

pb_z <- function(z, m.par) {
  linear.p <- m.par["repr.int"] + m.par["repr.z"] * z
  1/(1 + exp(-linear.p))
}

pr_z <- function(m.par) {
  linear.p <- m.par["recr.int"]
  1/(1 + exp(-linear.p))
}

c_z1z <- function(z1, z, m.par) {
  mu <- m.par["rcsz.int"] + m.par["rcsz.z"] * z
  dnorm(z1, mean = mu, sd = m.par["rcsz.sd"])
}

# --- Kernel functions ---
P_z1z <- function(z1, z, m.par) {
  s_z(z, m.par) * g_z1z(z1, z, m.par)
}

F_z1z <- function(z1, z, m.par) {
  s_z(z, m.par) * pb_z(z, m.par) * 0.5 * pr_z(m.par) * c_z1z(z1, z, m.par)
}

# --- Build kernel ---
mk_K <- function(m, m.par, L, U) {
  h <- (U - L)/m
  meshpts <- L + ((1:m) - 1/2) * h
  P <- h * outer(meshpts, meshpts, P_z1z, m.par = m.par)
  F <- h * outer(meshpts, meshpts, F_z1z, m.par = m.par)
  K <- P + F
  list(K = K, meshpts = meshpts, P = P, F = F, h = h)
}

# --- Compute ---
nBigMatrix <- 100
L <- 1.6
U <- 3.7

IPM.true <- mk_K(nBigMatrix, m.par.true, L, U)

eig <- eigen(IPM.true$K)
lambda.true <- Re(eig$values[1])
w.true <- Re(eig$vectors[, 1])
w.true <- w.true / sum(w.true)

cat("=== Ungulate DI Deterministic IPM (R) ===\n")
cat("Lambda:", lambda.true, "\n")
cat("Mean log body mass:", sum(w.true * IPM.true$meshpts), "\n")
cat("Var log body mass:", sum(w.true * IPM.true$meshpts^2) - sum(w.true * IPM.true$meshpts)^2, "\n")

# Age-specific size distributions
lam <- lambda.true
a0 <- IPM.true$F %*% w.true / lam
a1 <- IPM.true$P %*% a0 / lam
a2 <- IPM.true$P %*% a1 / lam
a3 <- IPM.true$P %*% a2 / lam

cat("Age 0 mean size:", sum((a0/sum(a0)) * IPM.true$meshpts), "\n")
cat("Age 1 mean size:", sum((a1/sum(a1)) * IPM.true$meshpts), "\n")
cat("Age 2 mean size:", sum((a2/sum(a2)) * IPM.true$meshpts), "\n")
cat("Age 3 mean size:", sum((a3/sum(a3)) * IPM.true$meshpts), "\n")
