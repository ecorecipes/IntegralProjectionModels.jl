"""
Kernel materialization: discretize continuous kernels into matrices using the midpoint rule.

For a kernel K(z', z) on domain [L, U] with n meshpoints:
- Bin edges: linspace(L, U, n+1)
- Midpoints: z[i] = (edges[i] + edges[i+1]) / 2
- Step size: h = (U - L) / n
- Matrix entry: M[i,j] = h * K(z[i], z[j])
"""

"""
    materialize(k::PKernel) -> Matrix

Discretize a survival-growth kernel using the midpoint rule.
`M[i,j] = h * survival(z[j]) * growth(z[i], z[j])`
"""
function materialize(k::PKernel)
    z = meshpoints(k.domain)
    h = step_size(k.domain)
    m = length(z)
    T = typeof(h * k.survival(z[1]) * k.growth(z[1], z[1]))
    K = zeros(T, m, m)

    if k.eviction == TruncatedDistributions
        for j in 1:m
            s_j = k.survival(z[j])
            for i in 1:m
                g_ij = truncated_growth(k.growth, z[i], z[j], k.domain)
                K[i, j] = s_j * g_ij * h
            end
        end
    else
        for j in 1:m
            s_j = k.survival(z[j])
            for i in 1:m
                K[i, j] = s_j * k.growth(z[i], z[j]) * h
            end
        end
    end

    if k.eviction == DiscreteExtrema
        apply_discrete_extrema!(K)
    end

    return K
end

"""
    materialize(k::FKernel) -> Matrix

Discretize a fecundity kernel using the midpoint rule.
`M[i,j] = h * fecundity(z[i], z[j])`
"""
function materialize(k::FKernel)
    z = meshpoints(k.domain)
    h = step_size(k.domain)
    m = length(z)
    T = typeof(h * k.fecundity(z[1], z[1]))
    K = zeros(T, m, m)

    for j in 1:m
        for i in 1:m
            K[i, j] = k.fecundity(z[i], z[j]) * h
        end
    end

    if k.eviction == DiscreteExtrema
        apply_discrete_extrema!(K)
    end

    return K
end

"""
    materialize(k::CustomKernel) -> Matrix

Discretize a custom kernel. The function is called as `func(z', z)`.
"""
function materialize(k::CustomKernel)
    z = meshpoints(k.domain)
    h = step_size(k.domain)
    m = length(z)
    T = typeof(h * k.func(z[1], z[1]))
    K = zeros(T, m, m)

    for j in 1:m
        for i in 1:m
            K[i, j] = k.func(z[i], z[j]) * h
        end
    end

    if k.eviction == DiscreteExtrema
        apply_discrete_extrema!(K)
    end

    return K
end

"""
    materialize(k::ComposedKernel) -> Matrix

Materialize a composed kernel by summing the materialized sub-kernels.
"""
function materialize(k::ComposedKernel)
    matrices = map(materialize, k.subkernels)
    return sum(matrices)
end

"""
    materialize(k::MatrixKernel) -> Matrix

Return the stored matrix directly.
"""
materialize(k::MatrixKernel) = k.matrix

"""
    materialize(mk::MegaKernel) -> Matrix

Materialize a mega-kernel into a block matrix for general models.

State ordering follows the order of `states`. Each continuous state
contributes `n_meshpoints` rows/columns; each discrete state contributes
`n_states` rows/columns.
"""
function materialize(mk::MegaKernel)
    state_names = keys(mk.states)
    state_sizes = [n_states(mk.states[s]) for s in state_names]
    total = sum(state_sizes)

    # Compute offsets for each state
    offsets = Dict{Symbol, Int}()
    offset = 0
    for (name, sz) in zip(state_names, state_sizes)
        offsets[name] = offset
        offset += sz
    end

    # Infer element type from sub-kernels
    Telm = isempty(mk.kernels) ? Float64 :
           promote_type((eltype(materialize(k)) for (_, k) in mk.kernels)...)
    K = zeros(Telm, total, total)

    for ((from, to), kernel) in mk.kernels
        M = materialize(kernel)
        r_start = offsets[to] + 1
        r_end = offsets[to] + n_states(mk.states[to])
        c_start = offsets[from] + 1
        c_end = offsets[from] + n_states(mk.states[from])
        K[r_start:r_end, c_start:c_end] .= M
    end

    return K
end
