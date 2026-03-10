## R cross-check for 08_age_size.qmd (Age x Size Soay sheep)
## Source: first-edition/Rcode/c6/Ungulate Age Demog Funs.R

# --- Parameters ---
m.par.true <- c(
  surv.int  = -17.0,
  surv.z    =  6.68,
  surv.a    = -0.334,
  grow.int  =  1.27,
  grow.z    =  0.612,
  grow.a    = -0.00724,
  grow.sd   =  0.0787,
  repr.int  = -7.88,
  repr.z    =  3.11,
  repr.a    = -0.078,
  recr.int  =  1.11,
  recr.a    =  0.184,
  rcsz.int  =  0.362,
  rcsz.z    =  0.709,
  rcsz.sd   =  0.159
)

# --- Vital rates (age-dependent) ---
s_z <- function(z, a, m.par) {
  1/(1 + exp(-(m.par["surv.int"] + m.par["surv.z"] * z + m.par["surv.a"] * a)))
}

g_z1z <- function(z1, z, a, m.par) {
  mu <- m.par["grow.int"] + m.par["grow.z"] * z + m.par["grow.a"] * a
  dnorm(z1, mean = mu, sd = m.par["grow.sd"])
}

pb_z <- function(z, a, m.par) {
  if (a == 0) return(rep(0, length(z)))  # age 0 cannot reproduce
  1/(1 + exp(-(m.par["repr.int"] + m.par["repr.z"] * z + m.par["repr.a"] * a)))
}

pr_z <- function(a, m.par) {
  1/(1 + exp(-(m.par["recr.int"] + m.par["recr.a"] * a)))
}

c_z1z <- function(z1, z, m.par) {
  dnorm(z1, m.par["rcsz.int"] + m.par["rcsz.z"] * z, m.par["rcsz.sd"])
}

# --- Domain ---
L <- 1.6; U <- 3.7; m <- 100; max_age <- 10

h <- (U - L)/m
meshpts <- L + ((1:m) - 1/2) * h

# --- Build age-specific P and F kernels ---
mk_P <- function(a, m.par) {
  P <- matrix(0, m, m)
  for (j in 1:m) {
    for (i in 1:m) {
      P[i, j] <- h * s_z(meshpts[j], a, m.par) * g_z1z(meshpts[i], meshpts[j], a, m.par)
    }
  }
  P
}

mk_F <- function(a, m.par) {
  F <- matrix(0, m, m)
  for (j in 1:m) {
    for (i in 1:m) {
      F[i, j] <- h * s_z(meshpts[j], a, m.par) * pb_z(meshpts[j], a, m.par) *
        0.5 * pr_z(a, m.par) * c_z1z(meshpts[i], meshpts[j], m.par)
    }
  }
  F
}

# --- Build block matrix ---
# Block structure:
# - F kernels in first row block (all ages produce age-1 offspring)
# - P kernels on sub-diagonal (aging from a to a+1)
# - P at max age on diagonal (absorbing)

total <- m * max_age
K <- matrix(0, total, total)

for (a in 1:max_age) {
  row_start <- (a - 1) * m + 1
  row_end <- a * m

  # Fecundity: age a produces offspring in age class 1
  F_a <- mk_F(a, m.par.true)
  K[1:m, row_start:row_end] <- K[1:m, row_start:row_end] + F_a

  # Survival-growth
  P_a <- mk_P(a, m.par.true)
  if (a < max_age) {
    # Age a -> age a+1 (sub-diagonal block)
    next_start <- a * m + 1
    next_end <- (a + 1) * m
    K[next_start:next_end, row_start:row_end] <- P_a
  } else {
    # Max age stays at max age (diagonal block)
    K[row_start:row_end, row_start:row_end] <- P_a
  }
}

# --- Eigenanalysis ---
eig <- eigen(K)
lam <- Re(eig$values[1])

cat("=== Age x Size IPM Cross-Check (R) ===\n")
cat("Lambda (age x size):", lam, "\n")
cat("Block matrix size:", nrow(K), "x", ncol(K), "\n")

# Stable age distribution
w <- Re(eig$vectors[, 1])
w <- w / sum(w)

age_totals <- numeric(max_age)
for (a in 1:max_age) {
  start <- (a - 1) * m + 1
  end <- a * m
  age_totals[a] <- sum(w[start:end])
}
cat("\nStable age distribution:\n")
for (a in 1:max_age) {
  cat(sprintf("  Age %2d: %.4f%%\n", a, age_totals[a] * 100))
}
