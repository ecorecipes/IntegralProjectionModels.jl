# Categorical Composition of IPMs

## Introduction

Integral Projection Models can become complex: survival, growth, fecundity, seed banks, age structure, spatial patches. **Applied category theory** provides a principled framework for composing such models from smaller, reusable parts.

Key ideas:

- **Sub-kernels as morphisms**: each vital rate process (survival-growth, fecundity, seed bank transitions) is a morphism in a Markov category
- **Operadic composition**: an undirected wiring diagram (UWD) specifies *how* sub-processes share trait axes; the composition rule sums their contributions
- **Left/Right Kan extensions**: translate between IPMs (continuous kernels) and MPMs (discrete matrices), with quantified approximation error
- **Stratification**: extend a local model to multiple spatial patches via a pullback construction
- **Functorial semantics**: derive multiple analyses (deterministic projection, sensitivity, stochastic simulation) from a single model specification

This vignette demonstrates these ideas using the Ligustrum general IPM from [vignette 07](07_general_ipm.qmd), which includes a discrete seed bank.

## Setup

```@example ipm
using IntegralProjectionModels
using Catlab
using Catlab.CategoricalAlgebra
using Catlab.WiringDiagrams
using Catlab.Programs
using Distributions
using LinearAlgebra

```

## The Ligustrum Model

We reuse the *Ligustrum obtusifolium* (border privet) model from vignette 07. The population has two states:

- **Height** (`ht`): continuous, domain $[1.02, 624]$ cm
- **Seed bank** (`b`): discrete, single pool

### Parameters

```@example ipm
ctrl = (
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
    g_i       = 0.5067,
)

cr = (
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
    g_i       = 0.5067,
)
```

### Domain and Vital Rates

```@example ipm
L = 1.02
U = 624.0
n = 500

ht_domain = ContinuousDomain(L, U, n)
b_domain = DiscreteDomain([:seeds])

z = meshpoints(ht_domain)
h = step_size(ht_domain)

# Vital rate functions
function survival(z, p)
    1.0 / (1.0 + exp(-(p.s_int + p.s_slope * z + p.s_slope_2 * z^2)))
end

function growth(z_prime, z, p)
    mu = p.g_int + p.g_slope * z
    d = Normal(mu, p.g_sd)
    ev = cdf(d, U) - cdf(d, L)
    return ev > 0 ? pdf(d, z_prime) / ev : pdf(d, z_prime)
end

function repr_prob(z, p)
    1.0 / (1.0 + exp(-(p.f_r_int + p.f_r_slope * z)))
end

function seed_prod(z, p)
    exp(p.f_s_int + p.f_s_slope * z)
end

function offspring_size(z_prime, p)
    d = Normal(p.f_d_mu, p.f_d_sd)
    ev = cdf(d, U) - cdf(d, L)
    return ev > 0 ? pdf(d, z_prime) / ev : pdf(d, z_prime)
end
```

## Sub-Kernels as Morphisms

In the Markov category **Krn**, each vital rate process is a morphism. We define the Ligustrum sub-kernels as standalone functions:

```@example ipm
# CC: survival-growth kernel P(z', z)
P_kernel(z_new, z, p) = survival(z, p) * growth(z_new, z, p)

# CC: fecundity kernel F(z', z) — direct recruitment (bypassing seed bank)
# (Not used in the general model, but included for the simple-IPM demo)
F_kernel(z_new, z, p) = repr_prob(z, p) * seed_prod(z, p) * (1 - p.g_i) *
    p.e_p * offspring_size(z_new, p)
```

## Operadic Composition via UWD (Simple IPM)

Before tackling the full general IPM, consider the simpler case where $K = P + F$ — a single continuous state with no seed bank. The decomposition into survival-growth and fecundity sub-kernels is captured by an **undirected wiring diagram** (UWD): two boxes share the same trait axes $(z, z')$, and the operad's composition rule sums their contributions.

```@example ipm
# Define the wiring pattern
demographic_uwd = @relation (z, z_new) begin
    survive_grow(z, z_new)
    reproduce(z, z_new)
end
```

```@example ipm
println("UWD for K = P + F:")
println("  Boxes:     $(length(boxes(demographic_uwd)))  (survive_grow, reproduce)")
println("  Junctions: $(length(junctions(demographic_uwd)))  (z, z_new)")
```

Evaluate by plugging in concrete sub-kernels:

```@example ipm
sub_kernels_simple = Dict(
    :survive_grow => (z_new, z) -> P_kernel(z_new, z, ctrl),
    :reproduce    => (z_new, z) -> F_kernel(z_new, z, ctrl),
)

K_simple = compose_from_uwd(demographic_uwd, sub_kernels_simple, ht_domain)

# Verify: same as compose_kernels without Catlab
K_simple_check = compose_kernels(sub_kernels_simple, ht_domain)
println("UWD vs direct composition match: ", isapprox(K_simple, K_simple_check))

λ_simple = maximum(real.(eigvals(K_simple)))
println("Simple IPM (P + F, no seed bank) λ = ", round(λ_simple, digits=4))
```

## General IPM as Block-Structured Composition

The full Ligustrum model has four transition blocks forming the mega-kernel:

$$
K = \begin{pmatrix} K_{CC} & K_{DC} \\ K_{CD} & K_{DD} \end{pmatrix}
$$

Each block is a morphism between different state spaces. We build this using the existing `MegaKernel` infrastructure, but the sub-kernels are defined as standalone functions — ready for categorical composition.

```@example ipm
function build_categorical_model(p)
    # CC: survival × growth (500×500)
    cc_kernel = PKernel(
        CustomVitalRate(z -> survival(z, p)),
        CustomVitalRate((z_prime, z) -> growth(z_prime, z, p)),
        ht_domain
    )

    # CD: continuous → discrete — seeds entering bank (1×500)
    cd_func(z) = repr_prob(z, p) * seed_prod(z, p) * p.g_i
    cd_values = [cd_func(zi) for zi in z]
    cd_matrix = reshape(cd_values, 1, n)

    # DC: discrete → continuous — seedlings emerging (500×1)
    dc_values = [p.e_p * offspring_size(zi, p) * h for zi in z]
    dc_matrix = reshape(dc_values, n, 1)

    # DD: discrete → discrete — stay in bank (1×1, = 0)
    dd_matrix = zeros(1, 1)

    states = (ht=ht_domain, b=b_domain)
    mega = MegaKernel(;
        states=states,
        ht_to_ht=cc_kernel,
        ht_to_b=MatrixKernel(cd_matrix),
        b_to_ht=MatrixKernel(dc_matrix),
        b_to_b=MatrixKernel(dd_matrix)
    )
    return mega
end

mega_ctrl = build_categorical_model(ctrl)
mega_cr = build_categorical_model(cr)
```

## Lowering to Concrete Matrices

Materialise the mega-kernels and verify against published values:

```@example ipm
K_ctrl = materialize(mega_ctrl)
K_cr = materialize(mega_cr)

e_ctrl = eigen(K_ctrl)
e_cr = eigen(K_cr)

λ_ctrl = real(e_ctrl.values[argmax(real.(e_ctrl.values))])
λ_cr = real(e_cr.values[argmax(real.(e_cr.values))])
```

```@example ipm
println("Control:  λ = ", round(λ_ctrl, digits=4), "  (published ≈ 1.22)")
println("CR:       λ = ", round(λ_cr, digits=4), "  (published ≈ 1.26)")
```

## Left Kan Extension: IPM → MPM

The **left Kan extension** discretises a continuous kernel function into a matrix using the midpoint rule. This is the standard IPM-to-MPM translation.

We demonstrate convergence of $\lambda$ as the mesh refines, using just the CC (survival-growth) block:

```@example ipm
# Reference: high-resolution CC block
cc_fn(z_new, z) = P_kernel(z_new, z, ctrl)
ref_domain = ContinuousDomain(L, U, 500)
A_ref = left_kan_extension(cc_fn, ref_domain)
λ_ref = maximum(real.(eigvals(A_ref)))

println("Convergence of λ for CC block (survival-growth):\n")
println("  n_bins │    λ      │  |Δλ| from ref")
println(" ────────┼───────────┼─────────────────")
for nb in [5, 10, 20, 50, 100, 200, 500]
    dom_nb = ContinuousDomain(L, U, nb)
    A_nb = left_kan_extension(cc_fn, dom_nb)
    λ_nb = maximum(real.(eigvals(A_nb)))
    Δ = abs(λ_nb - λ_ref)
    marker = nb == 500 ? " ← reference" : ""
    println("  $(lpad(nb, 5)) │ $(lpad(round(λ_nb, digits=6), 9)) │ $(lpad(round(Δ, sigdigits=2), 14))$marker")
end
```

## Right Kan Extension: MPM → IPM

The **right Kan extension** constructs a piecewise-constant kernel from a matrix — the most conservative continuous representation of an MPM:

$$K_{\mathrm{pw}}(z', z) = \frac{A_{ij}}{h} \quad \text{where } z \in B_j,\; z' \in B_i$$

```@example ipm
# Discretise at n=20, then reconstruct
dom_20 = ContinuousDomain(L, U, 20)
A_20 = left_kan_extension(cc_fn, dom_20)
K_pw = right_kan_extension(A_20, dom_20)

# Compare original vs piecewise at sample points
println("Round trip: IPM → MPM(n=20) → piecewise IPM\n")
println("  (z', z)          │ K_original  │ K_piecewise │  |Δ|")
println(" ──────────────────┼─────────────┼─────────────┼───────")
for (zn, zv) in [(50.0, 50.0), (100.0, 200.0), (300.0, 300.0), (500.0, 400.0)]
    k_orig = cc_fn(zn, zv)
    k_pw = K_pw(zn, zv)
    Δ = abs(k_orig - k_pw)
    zn_s = lpad(round(zn, digits=0), 5)
    zv_s = lpad(round(zv, digits=0), 5)
    println("  ($zn_s, $zv_s)     │ $(lpad(round(k_orig, digits=6), 11)) │ $(lpad(round(k_pw, digits=6), 11)) │ $(round(Δ, digits=6))")
end
```

### Adjunction Errors

The adjunction $\mathrm{Lan}_D \dashv D^* \dashv \mathrm{Ran}_D$ has a **unit** (MPM self-consistency) and **counit** (round-trip approximation quality):

```@example ipm
function adjunction_errors(kernel_fn, domain_n)
    A = left_kan_extension(kernel_fn, domain_n)
    K_pw = right_kan_extension(A, domain_n)
    A_roundtrip = left_kan_extension(K_pw, domain_n)

    # Unit error: should be ~0 for midpoint rule
    unit_err = norm(A_roundtrip - A) / norm(A)

    # Counit error: L2 on fine grid
    n_test = 200
    test_dom = ContinuousDomain(L, U, n_test)
    z_test = meshpoints(test_dom)
    h_test = step_size(test_dom)
    counit_err = 0.0
    for z in z_test, z_new in z_test
        counit_err += (K_pw(z_new, z) - kernel_fn(z_new, z))^2
    end
    counit_err = sqrt(counit_err * h_test^2)

    λ_n = maximum(real.(eigvals(A)))
    λ_err = abs(λ_n - λ_ref)

    return unit_err, counit_err, λ_err
end

println("  n_bins │ Unit (η)    │ Counit (ε)  │ |λ - λ_ref|")
println(" ────────┼─────────────┼─────────────┼────────────")
for nb in [5, 10, 20, 50, 100, 200]
    dom_nb = ContinuousDomain(L, U, nb)
    η, ε, λe = adjunction_errors(cc_fn, dom_nb)
    println("  $(lpad(nb, 5)) │ $(lpad(round(η, sigdigits=2), 11)) │ $(lpad(round(ε, sigdigits=2), 11)) │ $(lpad(round(λe, sigdigits=2), 10))")
end

println()
println("Unit η ≈ 0:  midpoint rule is self-consistent under round-trip")
println("Counit ε → 0: piecewise IPM converges to original as bins refine")
```

## Stratification: Spatial Extension

**Stratification** extends a single-patch model to multiple spatial patches. In the categorical framework, this is a **pullback** in a slice category: local demography is replicated across patches and weighted by dispersal probabilities.

```@example ipm
# Materialise the control CC block
K_cc_ctrl = materialize(PKernel(
    CustomVitalRate(z -> survival(z, ctrl)),
    CustomVitalRate((z_prime, z) -> growth(z_prime, z, ctrl)),
    ht_domain
))

# Symmetric 2-patch dispersal: 90% stay, 10% move
D_sym = [0.9 0.1; 0.1 0.9]
K_strat_sym = stratify(K_cc_ctrl, D_sym)
λ_strat_sym = maximum(real.(eigvals(K_strat_sym)))
```

```@example ipm
println("Symmetric dispersal (90% stay, 10% move):")
println("  Single-patch λ = ", round(λ_ref, digits=6))
println("  Two-patch λ    = ", round(λ_strat_sym, digits=6))
println("  Ratio: ", round(λ_strat_sym / λ_ref, digits=6),
    "  (≈ 1.0 for identical patches)")
```

```@example ipm
# Asymmetric dispersal: source-sink dynamics
D_asym = [0.95 0.20; 0.05 0.80]
K_strat_asym = stratify(K_cc_ctrl, D_asym)
λ_strat_asym = maximum(real.(eigvals(K_strat_asym)))
```

```@example ipm
println("Asymmetric dispersal (source-sink):")
println("  Two-patch λ = ", round(λ_strat_asym, digits=6))
```

## Coarsening via Pushforward

Given a fine discretisation, **coarsening** aggregates bins to produce a lower-resolution matrix. This is the pushforward (left Kan extension) along a coarsening map.

```@example ipm
# Fine: 100 bins, coarse: 20 bins
fine_domain = ContinuousDomain(L, U, 100)
coarse_domain = ContinuousDomain(L, U, 20)

A_fine = left_kan_extension(cc_fn, fine_domain)
λ_fine = maximum(real.(eigvals(A_fine)))

# Coarsen via domain-based pushforward
A_coarsened = coarsen(A_fine, fine_domain, coarse_domain)
λ_coarsened = maximum(real.(eigvals(A_coarsened)))

# Direct discretisation at n=20
A_direct_20 = left_kan_extension(cc_fn, coarse_domain)
λ_direct_20 = maximum(real.(eigvals(A_direct_20)))
```

```@example ipm
println("Fine (n=100):              λ = ", round(λ_fine, digits=6))
println("Coarsened via Lan (n=20):  λ = ", round(λ_coarsened, digits=6))
println("Direct at n=20:            λ = ", round(λ_direct_20, digits=6))
println("Reference (n=500):         λ = ", round(λ_ref, digits=6))
```

### Coarsening with Catlab FinFunction

The same coarsening can be expressed categorically using a Catlab `FinFunction`:

```@example ipm
n_fine = 100
n_coarse = 20
bins_per_coarse = n_fine ÷ n_coarse

# Coarsening map: FinSet(100) → FinSet(20)
coarsening_map = FinFunction(
    [div(i - 1, bins_per_coarse) + 1 for i in 1:n_fine],
    n_coarse
)

A_coarsened_ff = coarsen(A_fine, coarsening_map)
```

```@example ipm
println("Coarsened via FinFunction matches domain-based: ",
    isapprox(A_coarsened, A_coarsened_ff))
```

## Functorial Semantics

The same abstract demographic model can be interpreted via different **functors** — each producing a different analysis from the same specification. This is "one model, multiple analyses".

### F_det: Deterministic Projection

```@example ipm
# Dominant eigenvalue and stable distribution
idx = argmax(real.(e_ctrl.values))
w = abs.(real.(e_ctrl.vectors[:, idx]))
w ./= sum(w)

println("F_det (deterministic projection):")
println("  λ = ", round(λ_ctrl, digits=4))
println("  Population ", λ_ctrl > 1.0 ? "growing" : "declining")

# Find peak in height distribution (first n entries)
w_ht = w[1:n]
peak_idx = argmax(w_ht)
println("  Stable height distribution peak at z ≈ ", round(z[peak_idx], digits=1), " cm")
```

### F_sens: Sensitivity and Elasticity

```@example ipm
# Sensitivity matrix
e_ctrl_t = eigen(transpose(K_ctrl))
idx_t = argmax(real.(e_ctrl_t.values))
v = abs.(real.(e_ctrl_t.vectors[:, idx_t]))

S = (v * w') / dot(v, w)
E = (K_ctrl ./ λ_ctrl) .* S

# Decompose elasticity by block
n_ht = n
n_b = 1
E_CC = E[1:n_ht, 1:n_ht]
E_CD = E[(n_ht+1):(n_ht+n_b), 1:n_ht]
E_DC = E[1:n_ht, (n_ht+1):(n_ht+n_b)]
E_DD = E[(n_ht+1):(n_ht+n_b), (n_ht+1):(n_ht+n_b)]
```

```@example ipm
println("F_sens (elasticity analysis):")
println("  Total elasticity:  ", round(sum(E), digits=4), "  (should ≈ 1.0)")
println("  CC (survive-grow): ", round(sum(E_CC), digits=4),
    "  (", round(100*sum(E_CC), digits=1), "%)")
println("  CD (to seedbank):  ", round(sum(E_CD), digits=4),
    "  (", round(100*sum(E_CD), digits=1), "%)")
println("  DC (from seedbank):", round(sum(E_DC), digits=4),
    "  (", round(100*sum(E_DC), digits=1), "%)")
println("  DD (stay in bank): ", round(sum(E_DD), digits=4),
    "  (", round(100*sum(E_DD), digits=1), "%)")
```

### F_stoch: Stochastic Simulation

```@example ipm
# Stochastic projection using the full mega-kernel
n0 = vcat(uniform_population(ht_domain), [20.0])
prob = IPMProblem(GeneralIPM(), mega_ctrl, (ht=ht_domain, b=b_domain),
    n0, (0, 50))
sol = solve(prob)

# Per-step growth rates
mean_λ = exp(mean(log.(sol.lambdas[5:end])))
```

```@example ipm
println("F_stoch (stochastic projection):")
println("  Geometric mean λ (after burn-in) = ", round(mean_λ, digits=4))
println("  Asymptotic λ                     = ", round(λ_ctrl, digits=4))
```

## Commutation: Discretise then Stratify = Stratify then Discretise

Because the Kan extension is functorial, it commutes with stratification:

```@example ipm
n_comm = 50
dom_comm = ContinuousDomain(L, U, n_comm)

# Path 1: discretise first, then stratify
A_p1_local = left_kan_extension(cc_fn, dom_comm)
A_path1 = stratify(A_p1_local, D_sym)

# Path 2: stratify the continuous kernel, then discretise
z_comm = meshpoints(dom_comm)
h_comm = step_size(dom_comm)
n_patches = 2
n_total = n_comm * n_patches
A_path2 = zeros(n_total, n_total)
for p_from in 1:n_patches, p_to in 1:n_patches
    for j in 1:n_comm, i in 1:n_comm
        row = (p_to - 1) * n_comm + i
        col = (p_from - 1) * n_comm + j
        A_path2[row, col] = h_comm * D_sym[p_to, p_from] * cc_fn(z_comm[i], z_comm[j])
    end
end

comm_error = norm(A_path1 - A_path2) / norm(A_path1)
```

```@example ipm
println("Commutation: ||Path1 - Path2|| / ||Path1|| = ", round(comm_error, sigdigits=2))
println("Commutes: ", comm_error < 1e-12)
```

## Summary of Categorical Constructions

| Construction | Category theory | Function | What it does |
|---|---|---|---|
| Kernel composition | Operadic algebra (UWD) | `compose_from_uwd` / `compose_kernels` | Sum sub-kernels sharing trait axes |
| IPM → MPM | Left Kan extension | `left_kan_extension` | Discretise kernel → matrix |
| MPM → IPM | Right Kan extension | `right_kan_extension` | Matrix → piecewise kernel |
| Spatial model | Pullback (stratification) | `stratify` | Local dynamics + dispersal |
| Bin aggregation | Pushforward | `coarsen` | Fine → coarse matrix |
| Multiple analyses | Functorial semantics | `solve`, `sensitivity`, ... | One model, many outputs |

The `IntegralProjectionModels.categorical` module provides these operations without requiring Catlab. When Catlab is loaded, the extension adds UWD evaluation (`compose_from_uwd`) and `FinFunction`-based coarsening.
