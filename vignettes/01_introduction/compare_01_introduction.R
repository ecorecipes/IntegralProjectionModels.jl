## R cross-check for 01_introduction.qmd (Monocarp deterministic IPM)
## Source: first-edition/Rcode/c2/Monocarp Demog Funs.R + Monocarp Calculations.R

# --- Parameters ---
m.par.true <- c(
  surv.int  = -0.65,
  surv.z    =  0.75,
  flow.int  = -18.00,
  flow.z    =  6.9,
  grow.int  =  0.96,
  grow.z    =  0.59,
  grow.sd   =  0.67,
  rcsz.int  = -0.08,
  rcsz.sd   =  0.76,
  seed.int  =  1.00,
  seed.z    =  2.20,
  p.r       =  0.007
)

# --- Vital rate functions ---
s_z <- function(z, m.par) {
  linear.p <- m.par["surv.int"] + m.par["surv.z"] * z
  1/(1 + exp(-linear.p))
}

p_bz <- function(z, m.par) {
  linear.p <- m.par["flow.int"] + m.par["flow.z"] * z
  1/(1 + exp(-linear.p))
}

G_z1z <- function(z1, z, m.par) {
  mu <- m.par["grow.int"] + m.par["grow.z"] * z
  sig <- m.par["grow.sd"]
  dnorm(z1, mean = mu, sd = sig)
}

b_z <- function(z, m.par) {
  exp(m.par["seed.int"] + m.par["seed.z"] * z)
}

c_0z1 <- function(z1, m.par) {
  dnorm(z1, mean = m.par["rcsz.int"], sd = m.par["rcsz.sd"])
}

# --- Kernel functions ---
P_z1z <- function(z1, z, m.par) {
  (1 - p_bz(z, m.par)) * s_z(z, m.par) * G_z1z(z1, z, m.par)
}

F_z1z <- function(z1, z, m.par) {
  p_bz(z, m.par) * b_z(z, m.par) * m.par["p.r"] * c_0z1(z1, m.par)
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
nBigMatrix <- 250
L <- -2.65
U <- 4.5

IPM.true <- mk_K(nBigMatrix, m.par.true, L, U)

eig <- eigen(IPM.true$K)
lambda.true <- Re(eig$values[1])
w.true <- Re(eig$vectors[, 1])
w.true <- w.true / sum(w.true)

# Left eigenvector (reproductive value)
eig.left <- eigen(t(IPM.true$K))
v.true <- Re(eig.left$vectors[, 1])
v.true <- v.true / v.true[1]

cat("=== Monocarp Deterministic IPM (R) ===\n")
cat("Lambda:", lambda.true, "\n")
cat("Mean size (stable dist):", sum(w.true * IPM.true$meshpts), "\n")
cat("Max reproductive value:", max(v.true), "\n")

# Sensitivity and elasticity
vw <- sum(v.true * w.true)
S <- outer(v.true, w.true) / vw
E <- (IPM.true$K * S) / lambda.true
cat("Sum of elasticities:", sum(E), "\n")
cat("Elasticity of P:", sum(IPM.true$P * S / lambda.true), "\n")
cat("Elasticity of F:", sum(IPM.true$F * S / lambda.true), "\n")
