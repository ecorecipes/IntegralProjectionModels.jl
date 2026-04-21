# IPM ↔ MPM Discretization Bridge


## Introduction

A core insight of structured population modelling is that an integral
projection model (IPM) and a matrix projection model (MPM) are two ends
of the same mathematical object. An IPM describes population dynamics
via a continuous kernel $K(z', z)$ that integrates over a continuous
state variable $z$:

$$n(z', t+1) = \int_\Omega K(z', z)\, n(z, t)\, dz$$

When we discretise this integral using the midpoint rule on a mesh of
$m$ points, we obtain an $m \times m$ projection matrix $\mathbf{A}$ —
exactly the object a matrix projection model operates on.

This vignette makes that bridge explicit using
`IntegralProjectionModels.jl` and `MatrixProjectionModels.jl`. We:

1.  Define a biologically realistic IPM (survival/growth + fecundity).
2.  Solve it with `EigenAnalysis()` to materialise the kernel matrix.
3.  Wrap that matrix as a `MatrixProjectionModel`.
4.  Show that eigenanalysis results agree to machine precision.
5.  Demonstrate mesh-convergence to the true continuous-kernel
    eigenvalue.

``` julia
using IntegralProjectionModels
using MatrixProjectionModels
using StructuredPopulationCore
using LinearAlgebra

const IPM = IntegralProjectionModels
const MPM = MatrixProjectionModels
const SPC = StructuredPopulationCore
```

    StructuredPopulationCore

## Model specification

We use a classic size-structured perennial plant model (following
@easterling2000, @ellner2016). The state variable is log-size $z$.

**Vital-rate functions:**

| Process | Model | Parameters |
|----|----|----|
| Survival | Logistic | $\text{logit}(s) = \beta_0^s + \beta_1^s z$ |
| Growth | Normal | $\mu = \beta_0^g + \beta_1^g z$, $\sigma_g$ |
| Pr(flowering) × Seeds | Logistic × Exp | $\text{logit}(p_f) = a_0 + a_1 z$; $\mu_s = \exp(b_0 + b_1 z)$ |
| Recruit size | Normal | $\mu_r, \sigma_r$ |

``` julia
# --- Vital-rate parameters (plausible perennial herb) ---
# Survival
β0_s, β1_s = -1.5, 0.8     # intercept, slope on logit scale
# Growth
β0_g, β1_g = 0.5, 0.75     # growth regression
σ_g = 0.4                   # growth variance
# Fecundity (flowering × seed production × establishment)
# LogisticFecundityRate(repr_int, repr_slope, seed_int, seed_slope,
#                       recruit_mean, recruit_sd)
repr_int, repr_slope = -3.0, 1.2     # logistic flowering probability
seed_int, seed_slope = 0.0, 0.4      # log-linear seed count
μ_r, σ_r = 1.0, 0.3                  # recruit size distribution

# State domain
L, U_bound = 0.0, 7.0
```

    (0.0, 7.0)

## Building the IPM

``` julia
function build_ipm(n_mesh)
    domain = ContinuousDomain(L, U_bound, n_mesh)

    # Survival × Growth kernel (P)
    surv = IPM.LinearSurvival(β0_s, β1_s)
    grow = IPM.NormalGrowth(β0_g, β1_g, σ_g)
    P = IPM.PKernel(surv, grow, domain; eviction = IPM.TruncatedDistributions)

    # Fecundity kernel (F)
    fec = IPM.LogisticFecundityRate(repr_int, repr_slope, seed_int, seed_slope,
                                     μ_r, σ_r)
    F = IPM.FKernel(fec, domain)

    # Initial population
    n0 = IPM.normal_population(domain, 3.0, 1.0)

    K_composed = P + F

    prob = IPM.IPMProblem(SPC.SimpleIPM(), SPC.DensityIndependent(),
                          SPC.Deterministic(), K_composed, domain, n0, (0, 1))
    return prob, P, F
end

prob100, P100, F100 = build_ipm(100)
```

    (IPMProblem{SimpleIPM, DensityIndependent, Deterministic, ComposedKernel{Tuple{PKernel{LinearSurvival{Float64}, NormalGrowth{Float64}, ContinuousDomain{Float64}}, FKernel{LogisticFecundityRate{Float64}, ContinuousDomain{Float64}}}}, ContinuousDomain{Float64}, Vector{Float64}, Nothing, Nothing}(SimpleIPM(), DensityIndependent(), Deterministic(), ComposedKernel{Tuple{PKernel{LinearSurvival{Float64}, NormalGrowth{Float64}, ContinuousDomain{Float64}}, FKernel{LogisticFecundityRate{Float64}, ContinuousDomain{Float64}}}}((PKernel{LinearSurvival{Float64}, NormalGrowth{Float64}, ContinuousDomain{Float64}}(LinearSurvival{Float64}(-1.5, 0.8), NormalGrowth{Float64}(0.5, 0.75, 0.4), ContinuousDomain{Float64}(0.0, 7.0, 100), TruncatedDistributions), FKernel{LogisticFecundityRate{Float64}, ContinuousDomain{Float64}}(LogisticFecundityRate{Float64}(-3.0, 1.2, 0.0, 0.4, 1.0, 0.3), ContinuousDomain{Float64}(0.0, 7.0, 100), NoCorrection))), ContinuousDomain{Float64}(0.0, 7.0, 100), [0.0003448395330298932, 0.00042334150221797634, 0.0005171739028904952, 0.0006287157125329741, 0.0007605783929403077, 0.0009155996329340513, 0.0010968296786339213, 0.001307509141074552, 0.001551037261671132, 0.0018309297687676357  …  0.00010750949692738818, 8.491741425434166e-5, 6.674498658883192e-5, 5.220504758144273e-5, 4.0632944854401784e-5, 3.147139969399101e-5, 2.425636815348011e-5, 1.8604047920512736e-5, 1.4199108588460835e-5, 1.0784168587657666e-5], (0, 1), nothing, nothing, false, false), PKernel{LinearSurvival{Float64}, NormalGrowth{Float64}, ContinuousDomain{Float64}}(LinearSurvival{Float64}(-1.5, 0.8), NormalGrowth{Float64}(0.5, 0.75, 0.4), ContinuousDomain{Float64}(0.0, 7.0, 100), TruncatedDistributions), FKernel{LogisticFecundityRate{Float64}, ContinuousDomain{Float64}}(LogisticFecundityRate{Float64}(-3.0, 1.2, 0.0, 0.4, 1.0, 0.3), ContinuousDomain{Float64}(0.0, 7.0, 100), NoCorrection))

## Solve via EigenAnalysis

The `EigenAnalysis()` algorithm materialises the kernel matrix and
computes the full eigendecomposition in a single call.

``` julia
sol = IPM.solve(prob100, IPM.EigenAnalysis())

println("IPM λ = ", round(sol.eigenanalysis.lambda; digits = 8))
println("Retcode: ", sol.retcode)
```

    IPM λ = 0.78948836
    Retcode: Success

The materialised kernel matrix is stored in `sol.kernel_matrices`:

``` julia
K_matrix = sol.kernel_matrices
println("Type: ", typeof(K_matrix))
println("Size: ", size(K_matrix))
```

    Type: Matrix{Float64}
    Size: (100, 100)

## Wrapping as a Matrix Projection Model

We need the decomposition $\mathbf{A} = \mathbf{U} + \mathbf{F}$ for
life-history metrics. We materialise $P$ and $F$ separately:

``` julia
U_matrix = IPM.materialize(P100)   # survival/growth
F_matrix = IPM.materialize(F100)   # fecundity
A_check  = U_matrix + F_matrix

# Verify they sum to the total kernel
@assert K_matrix ≈ A_check "P + F should equal K"

mpm = MPM.MatrixProjectionModel(U_matrix, F_matrix)
println("MPM type: ", typeof(mpm))
```

    MPM type: MatrixProjectionModel{Float64, Matrix{Float64}}

## Eigenanalysis comparison

Since both the IPM and MPM operate on the **same discretised matrix**,
their eigenanalysis results should agree to machine precision.

``` julia
λ_ipm = sol.eigenanalysis.lambda
λ_mpm = MPM.lambda(mpm)
λ_spc = SPC.lambda(K_matrix)

println("λ (IPM EigenAnalysis):  ", round(λ_ipm; digits = 10))
println("λ (MPM):                ", round(λ_mpm; digits = 10))
println("λ (SPC on raw matrix):  ", round(λ_spc; digits = 10))
println()
println("|λ_ipm − λ_mpm| = ", abs(λ_ipm - λ_mpm))
println("|λ_ipm − λ_spc| = ", abs(λ_ipm - λ_spc))
```

    λ (IPM EigenAnalysis):  0.7894883596
    λ (MPM):                0.7894883596
    λ (SPC on raw matrix):  0.7894883596

    |λ_ipm − λ_mpm| = 0.0
    |λ_ipm − λ_spc| = 0.0

### Stable distribution and reproductive value

``` julia
w_ipm = sol.eigenanalysis.stable_dist
w_mpm = SPC.stable_distribution(K_matrix)

v_ipm = sol.eigenanalysis.repro_value
v_mpm = SPC.reproductive_value(K_matrix)

println("max |w_ipm − w_mpm| = ", maximum(abs.(w_ipm .- w_mpm)))
println("max |v_ipm − v_mpm| = ", maximum(abs.(v_ipm .- v_mpm)))
```

    max |w_ipm − w_mpm| = 0.0
    max |v_ipm − v_mpm| = 0.0

### Sensitivity and elasticity

``` julia
sens_ipm = SPC.sensitivity(K_matrix)
sens_mpm = MPM.sensitivity(mpm)
elas_ipm = SPC.elasticity(K_matrix)
elas_mpm = MPM.elasticity(mpm)

println("max |sens diff| = ", maximum(abs.(sens_ipm .- sens_mpm)))
println("max |elas diff| = ", maximum(abs.(elas_ipm .- elas_mpm)))
println("∑ elasticity (IPM) = ", round(sum(elas_ipm); digits = 6))
println("∑ elasticity (MPM) = ", round(sum(elas_mpm); digits = 6))
```

    max |sens diff| = 0.0
    max |elas diff| = 0.0
    ∑ elasticity (IPM) = 1.0
    ∑ elasticity (MPM) = 1.0

## Life-history metrics (MPM-specific)

The decomposition $\mathbf{A} = \mathbf{U} + \mathbf{F}$ unlocks metrics
that require separating survival from reproduction — these are available
in `MatrixProjectionModels.jl` but not in the IPM or SPC generic
interface.

``` julia
R0  = MPM.net_repro_rate(U_matrix, F_matrix)
T_g = MPM.gen_time(U_matrix, F_matrix)
println("Net reproductive rate R₀ = ", round(R0; digits = 4))
println("Generation time T        = ", round(T_g; digits = 4))
```

    Net reproductive rate R₀ = 0.0992
    Generation time T        = 9.7737

## Mesh convergence

The IPM is a numerical approximation of a continuous-kernel operator. As
we increase the mesh resolution, the eigenvalue should converge.

``` julia
mesh_sizes = [25, 50, 100, 200, 400]
lambdas = Float64[]

for m in mesh_sizes
    prob_m, _, _ = build_ipm(m)
    sol_m = IPM.solve(prob_m, IPM.EigenAnalysis())
    push!(lambdas, sol_m.eigenanalysis.lambda)
end

println("Mesh convergence:")
println("  m    λ            Δλ from m=400")
println("  " * "-"^40)
for (m, λ) in zip(mesh_sizes, lambdas)
    Δ = abs(λ - lambdas[end])
    println("  $(lpad(m, 3))  $(round(λ; digits=10))  $(round(Δ; sigdigits=3))")
end
```

    Mesh convergence:
      m    λ            Δλ from m=400
      ----------------------------------------
       25  0.78950298  1.57e-5
       50  0.7894917012  4.44e-6
      100  0.7894883596  1.09e-6
      200  0.7894874868  2.21e-7
      400  0.7894872662  0.0

### Convergence rate

The midpoint rule has error $O(h^2)$ where $h = (U - L) / m$. Doubling
the mesh should reduce the error by a factor of ≈4:

``` julia
for i in 2:length(mesh_sizes)-1
    Δ_prev = abs(lambdas[i-1] - lambdas[end])
    Δ_curr = abs(lambdas[i]   - lambdas[end])
    ratio  = Δ_prev / max(Δ_curr, eps())
    m_prev = mesh_sizes[i-1]
    m_curr = mesh_sizes[i]
    println("  m=$m_prev → m=$m_curr: error ratio = $(round(ratio; digits=2))" *
            "  (expected ≈ $(round((m_curr/m_prev)^2; digits=1)))")
end
```

      m=25 → m=50: error ratio = 3.54  (expected ≈ 4.0)
      m=50 → m=100: error ratio = 4.06  (expected ≈ 4.0)
      m=100 → m=200: error ratio = 4.96  (expected ≈ 4.0)

## Summary

| Quantity | IPM | MPM | Agreement |
|----|----|----|----|
| λ | `solve(prob, EigenAnalysis())` | `MPM.lambda(mpm)` | Machine precision |
| $\mathbf{w}$ | `sol.eigenanalysis.stable_dist` | `SPC.stable_distribution(A)` | Machine precision |
| $\mathbf{v}$ | `sol.eigenanalysis.repro_value` | `SPC.reproductive_value(A)` | Machine precision |
| Sensitivity | `SPC.sensitivity(A)` | `MPM.sensitivity(mpm)` | Machine precision |
| Elasticity | `SPC.elasticity(A)` | `MPM.elasticity(mpm)` | Machine precision |
| $R_0$, $T$ | — | `MPM.net_repro_rate(U, F)` | MPM-specific (needs U/F decomposition) |

**Key insight:** The IPM–MPM bridge is mathematically exact at any given
mesh resolution. The midpoint-rule discretisation is the *only*
approximation, and it is shared by both representations. The choice
between IPM and MPM is one of modelling convenience, not of analytical
accuracy.

## References

- Easterling, M. R., Ellner, S. P., & Dixon, P. M. (2000). Size-specific
  sensitivity: applying a new structured population model. *Ecology*,
  81, 694–708.
- Ellner, S. P., Childs, D. Z., & Rees, M. (2016). *Data-driven
  Modelling of Structured Populations*. Springer.
