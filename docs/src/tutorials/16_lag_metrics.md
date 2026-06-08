# Time-Lag Metrics: Augmentation, R0, and Generation Time

## Overview

The `12_time_lag` vignette built a single-lag `LaggedKernel` and computed
its dominant eigenvalue via the augmented block matrix.  This vignette
focuses on the **lower-level helpers** that underpin lagged analysis:

- `expand_lag_matrix`  — assemble the $(L{+}1)m \times (L{+}1)m$ augmented
  block matrix from component kernels;
- `extract_lag_components` — split an augmented matrix back into its
  per-lag blocks;
- `augment_population` / `extract_population` — replicate / project the
  population state in / out of the augmented representation;
- `net_repro_rate_lagged` — net reproductive rate $R_0$ from
  $\rho(F_{\text{aug}} (I - U_{\text{aug}})^{-1})$ (Kuss et al. 2008);
- `generation_time_lagged` — generation time $T = \log R_0 / \log \lambda$
  for a lagged model.

These helpers are useful when you want to skip the high-level
`LaggedKernel` wrapper and work directly with assembled matrices, for
example inside larger composed models or sensitivity analyses.

## Setup

```@example ipm
using IntegralProjectionModels
using LinearAlgebra

```

## Build P and F matrices for a 4-bin model

```@example ipm
m = 4
P = [0.30 0.10 0.00 0.00;
     0.40 0.50 0.20 0.00;
     0.00 0.30 0.50 0.30;
     0.00 0.00 0.20 0.60]
F = [0.00 0.20 0.80 1.50;
     0.00 0.00 0.00 0.00;
     0.00 0.00 0.00 0.00;
     0.00 0.00 0.00 0.00]
println("P size = ", size(P), "   F size = ", size(F))
```

## Augmented matrix from P and F

`expand_lag_matrix(P, F)` assembles the standard single-lag block matrix
$\begin{bmatrix} P & F \\ I & 0 \end{bmatrix}$.

```@example ipm
K_aug = expand_lag_matrix(P, F)
println("K_aug size = ", size(K_aug))
K_aug
```

## Multi-lag construction

`expand_lag_matrix(lag_kernels, lag_structure)` is the general form for
$L>1$.  `lag_kernels[1]` is the immediate (lag-0) kernel and
`lag_kernels[k+1]` acts on $n(t-k)$.

```@example ipm
F1 = F .* 0.4   # fecundity at lag 1
F2 = F .* 0.6   # fecundity at lag 2
lag_struct = TimeLagStructure(2)
K_aug2 = expand_lag_matrix([P, F1, F2], lag_struct)
println("K_aug2 size = ", size(K_aug2), "   (expected $(3m)×$(3m))")
```

## Decomposing an augmented matrix

`extract_lag_components` is the inverse of `expand_lag_matrix`: given the
block matrix and the per-population dimension `m`, it returns a NamedTuple
with the per-lag kernels.

```@example ipm
parts = extract_lag_components(K_aug2, m, lag_struct)
println("typeof(parts).name.name = ", typeof(parts).name.name)
println("# kernels recovered     = ", length(parts.kernels))
println("kernels[1] ≈ P          = ", isapprox(parts.kernels[1], P))
println("kernels[2] ≈ F1         = ", isapprox(parts.kernels[2], F1))
println("kernels[3] ≈ F2         = ", isapprox(parts.kernels[3], F2))
```

## Population augmentation

`augment_population` replicates a physical state vector $n_0$ into the
$(L{+}1)m$-dimensional augmented representation (every history slot
initialised with $n_0$).  `extract_population` is the projection back to
the first $m$ entries.

```@example ipm
n0 = [10.0, 5.0, 2.0, 1.0]
n_aug = augment_population(n0, lag_struct)
println("length(n_aug)  = ", length(n_aug))
println("extract back  = ", extract_population(n_aug, m))
println("matches n0?   = ", extract_population(n_aug, m) == n0)
```

A single iteration of the augmented system reproduces the lagged update:

```@example ipm
n_next_aug = K_aug2 * n_aug
println("n_{t+1} (physical) = ", round.(extract_population(n_next_aug, m); digits=3))
```

## R0 and generation time

For the simple single-lag model:

```@example ipm
R0 = net_repro_rate_lagged(P, F)
T  = generation_time_lagged(P, F)
λ  = lambda(K_aug)
println("λ  = ", round(λ,  digits=4))
println("R0 = ", round(R0, digits=4))
println("T  = ", round(T,  digits=4))
```

Sanity check — the textbook identity $\log R_0 = T \log \lambda$:

```@example ipm
println("log(R0)         = ", round(log(R0),    digits=5))
println("T * log(λ)      = ", round(T * log(λ), digits=5))
```

For a multi-lag model, pass the kernel vector and the lag structure:

```@example ipm
R0_multi = net_repro_rate_lagged([P, F1, F2], lag_struct)
println("R0 (multi-lag) = ", round(R0_multi, digits=4))
```

## Summary

- `expand_lag_matrix` and `extract_lag_components` form an inverse pair for
  switching between the per-lag block representation and the augmented
  matrix used by direct iteration / eigenanalysis.
- `augment_population` / `extract_population` keep the population state
  consistent between the two representations.
- `net_repro_rate_lagged` (single-lag and multi-lag) and
  `generation_time_lagged` adapt classical fundamental-matrix demography
  to lagged models.  They satisfy the identity
  $\log R_0 = T \log \lambda$ as expected.
