## R cross-check for 05_density_dependent.qmd (Soay DD IPM)
## Source: first-edition/Rcode/c5/Soay DD Demog Funs.R + Soay DD IPM calcs.R

# --- Parameters ---
m.par.true <- c(
  surv.int  =  1.06,
  surv.z    =  2.09,
  surv.Nt   = -0.018,
  grow.int  =  1.41,
  grow.z    =  0.557,
  grow.sd   =  0.0799,
  repr.int  = -7.23,
  repr.z    =  2.60,
  recr.int  =  4.43,
  recr.Nt   = -0.00919,
  rcsz.int  =  0.540,
  rcsz.z    =  0.710,
  rcsz.Nt   = -0.000642,
  rcsz.sd   =  0.159
)

# --- Vital rates ---
s_z <- function(z, Nt, m.par) {
  1/(1 + exp(-(m.par["surv.int"] + m.par["surv.z"] * z + m.par["surv.Nt"] * Nt)))
}
g_z1z <- function(z1, z, m.par) {
  dnorm(z1, m.par["grow.int"] + m.par["grow.z"] * z, m.par["grow.sd"])
}
pb_z <- function(z, m.par) {
  1/(1 + exp(-(m.par["repr.int"] + m.par["repr.z"] * z)))
}
pr_z <- function(Nt, m.par) {
  1/(1 + exp(-(m.par["recr.int"] + m.par["recr.Nt"] * Nt)))
}
c_z1z <- function(z1, z, Nt, m.par) {
  dnorm(z1, m.par["rcsz.int"] + m.par["rcsz.z"] * z + m.par["rcsz.Nt"] * Nt, m.par["rcsz.sd"])
}

# --- Kernel ---
P_z1z <- function(z1, z, Nt, m.par) s_z(z, Nt, m.par) * g_z1z(z1, z, m.par)
F_z1z <- function(z1, z, Nt, m.par) 0.5 * s_z(z, Nt, m.par) * pb_z(z, m.par) * pr_z(Nt, m.par) * c_z1z(z1, z, Nt, m.par)

mk_K <- function(Nt, m, m.par, L, U) {
  h <- (U - L)/m
  meshpts <- L + ((1:m) - 1/2) * h
  P <- h * outer(meshpts, meshpts, P_z1z, Nt = Nt, m.par = m.par)
  F <- h * outer(meshpts, meshpts, F_z1z, Nt = Nt, m.par = m.par)
  K <- P + F
  list(K = K, meshpts = meshpts, P = P, F = F, h = h)
}

# --- Domain ---
L <- 0.5; U <- 3.8; m <- 165

# --- Lambda as function of N ---
lamFun <- function(Nt) {
  IPM.dd <- mk_K(Nt, m, m.par.true, L, U)
  Re(eigen(IPM.dd$K)$values[1])
}

# --- Find equilibrium ---
Nbar <- uniroot(function(N) lamFun(N) - 1, lower = 200, upper = 400)$root

cat("=== Density-Dependent IPM Cross-Check (R) ===\n")
cat("Equilibrium N:", Nbar, "females\n")
cat("Lambda at equilibrium:", lamFun(Nbar), "\n")
cat("Lambda at N=50:", lamFun(50), "\n")
cat("Lambda at N=450:", lamFun(450), "\n")
