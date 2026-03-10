"""
Time-lagged integral projection models.

Wraps IPM kernels with lag structure information. When materialized,
produces the augmented `(L+1)m × (L+1)m` block matrix.
"""

"""
    LaggedKernel{K, D}

An IPM kernel with time-lagged components.

# Fields
- `immediate::K`: Kernel applied to n(t) (lag 0, typically survival/growth P)
- `lagged::D`: Dict mapping lag index → kernel applied to n(t-lag)
- `lag_structure::TimeLagStructure`: Maximum lag specification
"""
struct LaggedKernel{K, D}
    immediate::K
    lagged::D
    lag_structure::TimeLagStructure
end

"""
    LaggedKernel(immediate, lagged_kernel; lag=1)

Convenience constructor for the common single-lag case.

`n(t+1) = P·n(t) + F·n(t-lag)` where `immediate = P` and `lagged_kernel = F`.
"""
function LaggedKernel(immediate, lagged_kernel; lag::Int=1)
    lag > 0 || throw(ArgumentError("lag must be positive"))
    lagged = Dict{Int, typeof(lagged_kernel)}(lag => lagged_kernel)
    LaggedKernel(immediate, lagged, TimeLagStructure(lag))
end

function Base.show(io::IO, lk::LaggedKernel)
    L = lk.lag_structure.max_lag
    n_lagged = length(lk.lagged)
    print(io, "LaggedKernel(max_lag=$L, $n_lagged lagged component(s))")
end

"""
    materialize(lk::LaggedKernel)

Materialize a `LaggedKernel` into an `(L+1)m × (L+1)m` augmented block matrix.

The immediate kernel is materialized to form the lag-0 block, and each
lagged kernel is materialized to form its respective lag block. Identity
matrices form the sub-diagonal shifts.
"""
function materialize(lk::LaggedKernel)
    L = lk.lag_structure.max_lag

    # Materialize the immediate kernel to determine m
    K0 = materialize(lk.immediate)
    m = size(K0, 1)
    T = eltype(K0)

    # Build lag kernel vector [K_0, K_1, ..., K_L]
    lag_kernels = Vector{Matrix{T}}(undef, L + 1)
    lag_kernels[1] = K0
    for k in 1:L
        if haskey(lk.lagged, k)
            lag_kernels[k + 1] = Matrix{T}(materialize(lk.lagged[k]))
        else
            lag_kernels[k + 1] = zeros(T, m, m)
        end
    end

    return expand_lag_matrix(lag_kernels, lk.lag_structure)
end

"""
    expand_lag_kernels(P_kernel, F_kernel, domain::ContinuousDomain)

Convenience: materialize P and F kernels, then build the `2m × 2m`
augmented matrix for a single-lag model.
"""
function expand_lag_kernels(P_kernel, F_kernel, domain::ContinuousDomain)
    P = materialize(P_kernel)
    F = materialize(F_kernel)
    return expand_lag_matrix(P, F)
end

# --- Solve dispatches for LaggedKernel ---

_is_lagged_kernel(::LaggedKernel) = true

function _solve_lagged_ipm(prob::IPMProblem)
    lk = prob.kernel::LaggedKernel
    K_aug = materialize(lk)
    L = lk.lag_structure.max_lag

    # Extract component matrices
    K0 = materialize(lk.immediate)
    m = size(K0, 1)
    components = extract_lag_components(K_aug, m, lk.lag_structure)
    lag_kernels = components.kernels

    t0, tf = prob.tspan
    n_steps = tf - t0

    # Initialize history: all lag slots start with n0
    history = [copy(float.(prob.n0)) for _ in 0:L]

    u = Vector{typeof(float.(prob.n0))}(undef, n_steps + 1)
    u[1] = copy(float.(prob.n0))
    lambdas = Vector{Float64}(undef, n_steps)

    for t in 1:n_steps
        n_new = lag_kernels[1] * history[1]
        for k in 1:L
            n_new .+= lag_kernels[k + 1] * history[k + 1]
        end

        pop_t = sum(history[1])
        pop_new = sum(n_new)
        lambdas[t] = pop_t > 0 ? pop_new / pop_t : 0.0

        if prob.normalize && pop_new > 0
            n_new ./= pop_new
        end

        # Shift history
        for k in L:-1:1
            history[k + 1] = history[k]
        end
        history[1] = n_new

        u[t + 1] = copy(n_new)
    end

    eigen_result = eigenanalysis_power(K_aug)
    return IPMSolution(collect(t0:tf), u, K_aug, eigen_result, :Success, lambdas)
end

function _solve_lagged_ipm_eigen(prob::IPMProblem)
    lk = prob.kernel::LaggedKernel
    K_aug = materialize(lk)
    ea = eigenanalysis_full(K_aug)
    u = [copy(prob.n0)]
    return IPMSolution([prob.tspan[1]], u, K_aug, ea, :Success, Float64[])
end
