"""
Categorical composition tools for Integral Projection Models.

Provides general-purpose "lowering" functions inspired by applied category theory:

- **Left Kan extension**: kernel function → matrix (IPM → MPM discretisation)
- **Right Kan extension**: matrix → piecewise-constant kernel (MPM → IPM refinement)
- **Stratification**: local kernel + dispersal → block-structured spatial matrix
- **Coarsening**: fine matrix → coarse matrix via bin aggregation (pushforward)
- **Kernel composition**: sum sub-kernels additively (UWD evaluation without Catlab)

These functions have no dependency on Catlab.jl; the Catlab extension
(`IntegralProjectionModelsCatlabExt`) adds UWD- and FinFunction-aware methods.
"""

"""
    left_kan_extension(kernel_fn, domain::ContinuousDomain)

Discretise a continuous kernel function `kernel_fn(z_new, z)` into a matrix
using the midpoint rule on the given domain.

This is the **left Kan extension** along the discretisation functor:
it maps a morphism in **Krn** (Markov kernels on measurable spaces) to a
morphism in **FinStoch** (stochastic matrices on finite sets).

``A_{ij} = h \\cdot K(z_i, z_j)``

where ``h`` is the bin width and ``z_i`` are the midpoints.

# Arguments
- `kernel_fn`: a callable `(z_new, z) -> Real` giving the kernel density
- `domain`: a `ContinuousDomain` specifying bounds and mesh resolution

# Returns
- `A::Matrix{Float64}`: the discretised transition matrix
"""
function left_kan_extension(kernel_fn, domain::ContinuousDomain)
    z = meshpoints(domain)
    h = step_size(domain)
    m = length(z)
    A = zeros(m, m)
    for j in 1:m
        for i in 1:m
            A[i, j] = h * kernel_fn(z[i], z[j])
        end
    end
    return A
end

"""
    right_kan_extension(A::Matrix, domain::ContinuousDomain)

Construct a piecewise-constant kernel function from a matrix and domain.

This is the **right Kan extension**: the most conservative continuous
representation of a discrete transition matrix. The resulting kernel is
constant within each bin pair:

``K_{\\mathrm{pw}}(z', z) = A_{ij} / h \\quad \\text{where } z \\in B_j,\\; z' \\in B_i``

# Arguments
- `A`: an `n × n` transition matrix
- `domain`: a `ContinuousDomain` with `n` meshpoints matching `size(A, 1)`

# Returns
- A callable `(z_new, z) -> Real` giving the piecewise-constant kernel density
"""
function right_kan_extension(A::Matrix, domain::ContinuousDomain)
    z_min = domain.lower
    h = step_size(domain)
    n = size(A, 1)
    function piecewise_kernel(z_new, z)
        j = clamp(Int(ceil((z - z_min) / h)), 1, n)
        i = clamp(Int(ceil((z_new - z_min) / h)), 1, n)
        return A[i, j] / h
    end
    return piecewise_kernel
end

"""
    stratify(A_local::Matrix, dispersal::Matrix)

Construct a block-structured spatial matrix from a local transition matrix
and a dispersal (connectivity) matrix.

In the categorical framework, this is a **pullback** in a slice category:
the local demography is replicated across patches and weighted by the
dispersal kernel.

``A_{\\mathrm{strat}}[(p_{\\mathrm{to}}, i),\\, (p_{\\mathrm{from}}, j)]
  = D[p_{\\mathrm{to}}, p_{\\mathrm{from}}] \\cdot A_{\\mathrm{local}}[i, j]``

# Arguments
- `A_local`: an `n × n` local (single-patch) transition matrix
- `dispersal`: a `p × p` matrix of dispersal probabilities between patches

# Returns
- `A_strat::Matrix{Float64}`: an `(n·p) × (n·p)` block-structured matrix
"""
function stratify(A_local::Matrix, dispersal::Matrix)
    n_bins = size(A_local, 1)
    n_patches = size(dispersal, 1)
    n_total = n_bins * n_patches
    A_strat = zeros(n_total, n_total)
    for p_to in 1:n_patches
        for p_from in 1:n_patches
            rows = ((p_to - 1) * n_bins + 1):(p_to * n_bins)
            cols = ((p_from - 1) * n_bins + 1):(p_from * n_bins)
            A_strat[rows, cols] = dispersal[p_to, p_from] * A_local
        end
    end
    return A_strat
end

"""
    coarsen(A::Matrix, from_domain::ContinuousDomain, to_domain::ContinuousDomain)

Coarsen a fine transition matrix to a coarser resolution by aggregating bins.

This is the **pushforward** (left Kan extension) along the coarsening map
``f\\colon \\mathrm{FinSet}(n_{\\mathrm{fine}}) \\to \\mathrm{FinSet}(n_{\\mathrm{coarse}})``.
Each coarse bin collects contributions from the fine bins it contains.

The `from_domain` and `to_domain` must share the same `[lower, upper]` bounds.
The number of fine bins must be an integer multiple of the number of coarse bins.

# Arguments
- `A`: an `n_fine × n_fine` transition matrix
- `from_domain`: the fine `ContinuousDomain`
- `to_domain`: the coarse `ContinuousDomain`

# Returns
- `A_coarse::Matrix{Float64}`: an `n_coarse × n_coarse` transition matrix
"""
function coarsen(A::Matrix, from_domain::ContinuousDomain, to_domain::ContinuousDomain)
    n_fine = from_domain.n_meshpoints
    n_coarse = to_domain.n_meshpoints
    n_fine == size(A, 1) || throw(DimensionMismatch(
        "Matrix size $(size(A, 1)) does not match from_domain meshpoints $n_fine"))
    n_fine % n_coarse == 0 || throw(ArgumentError(
        "n_fine ($n_fine) must be a multiple of n_coarse ($n_coarse)"))
    bins_per_coarse = n_fine ÷ n_coarse

    # Build coarsening map: fine bin i → coarse bin ceil(i / bins_per_coarse)
    A_coarse = zeros(n_coarse, n_coarse)
    for j_fine in 1:n_fine
        c_j = (j_fine - 1) ÷ bins_per_coarse + 1
        for i_fine in 1:n_fine
            c_i = (i_fine - 1) ÷ bins_per_coarse + 1
            A_coarse[c_i, c_j] += A[i_fine, j_fine] / bins_per_coarse
        end
    end
    return A_coarse
end

"""
    compose_kernels(sub_kernels, domain::ContinuousDomain)

Compose sub-kernel functions additively and discretise onto the given domain.

This is the concrete evaluation of an **undirected wiring diagram** (UWD)
without requiring Catlab: each sub-kernel shares the same trait axes, and
the composition rule sums their contributions.

``K_{ij} = \\sum_k h \\cdot K_k(z_i, z_j)``

# Arguments
- `sub_kernels`: a `NamedTuple`, `Dict`, or iterable of `(name, fn)` pairs
  where each `fn` is a callable `(z_new, z) -> Real`
- `domain`: a `ContinuousDomain`

# Returns
- `K::Matrix{Float64}`: the summed, discretised transition matrix
"""
function compose_kernels(sub_kernels, domain::ContinuousDomain)
    z = meshpoints(domain)
    h = step_size(domain)
    m = length(z)
    K = zeros(m, m)
    for (_, kfn) in pairs(sub_kernels)
        for j in 1:m
            for i in 1:m
                K[i, j] += h * kfn(z[i], z[j])
            end
        end
    end
    return K
end

"""
    compose_from_uwd(uwd, sub_kernels::Dict, domain::ContinuousDomain)

Evaluate a Catlab undirected wiring diagram (UWD) with concrete sub-kernels.

This function is extended by `IntegralProjectionModelsCatlabExt` when Catlab
is loaded. Without Catlab, calling this function will throw an error.

# Arguments
- `uwd`: a Catlab UWD (from `@relation`)
- `sub_kernels`: a `Dict{Symbol, Function}` mapping box names to kernel functions
- `domain`: a `ContinuousDomain`

# Returns
- `K::Matrix{Float64}`: the composed, discretised transition matrix

See also: [`compose_kernels`](@ref) for a Catlab-free alternative.
"""
function compose_from_uwd end
