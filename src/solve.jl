"""
    IPMSolution

Result of solving an IPMProblem.

# Fields
- `t`: time steps
- `u`: population state at each time (Vector of Vectors)
- `kernel_matrices`: materialized kernel(s)
- `eigenanalysis`: named tuple `(lambda, stable_dist, repro_value)` or `nothing`
- `retcode`: `:Success` or `:MaxIters`
- `lambdas`: per-step growth rates (for stochastic models)
"""
struct IPMSolution{T, U, K, E} <: AbstractProjectionSolution
    t::T
    u::U
    kernel_matrices::K
    eigenanalysis::E
    retcode::Symbol
    lambdas::Vector{Float64}
end

# --- Direct Iteration ---

"""
    solve(prob::IPMProblem, alg::DirectIteration; kwargs...)

Solve an IPM by iterating `n_{t+1} = K * n_t`.
"""
function CommonSolve.solve(prob::IPMProblem, alg::DirectIteration = DirectIteration();
        kwargs...)
    if _is_lagged_kernel(prob.kernel)
        return _solve_lagged_ipm(prob)
    end
    _solve(prob.structure, prob.density, prob.stochasticity, prob, alg; kwargs...)
end

# Fallback: non-lagged kernels
_is_lagged_kernel(k) = false
# _is_lagged_kernel(::LaggedKernel) = true is defined in time_lag.jl

# Simple, density-independent, deterministic
function _solve(::AbstractIPMStructure, ::DensityIndependent, ::Deterministic,
        prob::IPMProblem, ::DirectIteration; kwargs...)
    K = materialize(prob.kernel)
    t0, tf = prob.tspan
    n_steps = tf - t0

    u = Vector{typeof(prob.n0)}(undef, n_steps + 1)
    u[1] = copy(prob.n0)
    lambdas = Vector{Float64}(undef, n_steps)

    for t in 1:n_steps
        n_new = K * u[t]
        pop_t = sum(u[t])
        pop_new = sum(n_new)
        lambdas[t] = pop_t > 0 ? pop_new / pop_t : 0.0
        if prob.normalize && pop_new > 0
            n_new ./= pop_new
        end
        u[t + 1] = n_new
    end

    eigen_result = eigenanalysis_power(K)
    return IPMSolution(collect(t0:tf), u, K, eigen_result, :Success, lambdas)
end

# Density-dependent (deterministic or stochastic)
function _solve(::AbstractIPMStructure, ::DensityDependent, ::Deterministic,
        prob::IPMProblem, ::DirectIteration; kwargs...)
    t0, tf = prob.tspan
    n_steps = tf - t0

    u = Vector{typeof(prob.n0)}(undef, n_steps + 1)
    u[1] = copy(prob.n0)
    lambdas = Vector{Float64}(undef, n_steps)
    kernel_matrices = Vector{AbstractMatrix}(undef, n_steps)

    for t in 1:n_steps
        # For DD models, kernel is a function of current pop state
        K = _build_dd_kernel(prob, u[t], t + t0 - 1)
        kernel_matrices[t] = K
        n_new = K * u[t]
        pop_t = sum(u[t])
        pop_new = sum(n_new)
        lambdas[t] = pop_t > 0 ? pop_new / pop_t : 0.0
        if prob.normalize && pop_new > 0
            n_new ./= pop_new
        end
        u[t + 1] = n_new
    end

    return IPMSolution(collect(t0:tf), u, kernel_matrices, nothing, :Success, lambdas)
end

# Stochastic kernel-resampled
function _solve(::AbstractIPMStructure, ::Any, ::StochasticKernelResampled,
        prob::IPMProblem, ::DirectIteration;
        kernel_seq = nothing, kwargs...)
    t0, tf = prob.tspan
    n_steps = tf - t0

    # prob.kernel should be a vector/tuple of kernels or a function returning indexed kernels
    kernel_set = _build_kernel_set(prob)
    n_kernels = length(kernel_set)

    if kernel_seq === nothing
        kernel_seq = rand(1:n_kernels, n_steps)
    end

    u = Vector{typeof(prob.n0)}(undef, n_steps + 1)
    u[1] = copy(prob.n0)
    lambdas = Vector{Float64}(undef, n_steps)
    used_kernels = Vector{AbstractMatrix}(undef, n_steps)

    for t in 1:n_steps
        K = kernel_set[kernel_seq[t]]
        used_kernels[t] = K
        n_new = K * u[t]
        pop_t = sum(u[t])
        pop_new = sum(n_new)
        lambdas[t] = pop_t > 0 ? pop_new / pop_t : 0.0
        if prob.normalize && pop_new > 0
            n_new ./= pop_new
        end
        u[t + 1] = n_new
    end

    return IPMSolution(collect(t0:tf), u, used_kernels, nothing, :Success, lambdas)
end

# Stochastic parameter-resampled
function _solve(::AbstractIPMStructure, ::Any, ::StochasticParameterResampled,
        prob::IPMProblem, ::DirectIteration; kwargs...)
    t0, tf = prob.tspan
    n_steps = tf - t0

    u = Vector{typeof(prob.n0)}(undef, n_steps + 1)
    u[1] = copy(prob.n0)
    lambdas = Vector{Float64}(undef, n_steps)
    kernel_matrices = Vector{AbstractMatrix}(undef, n_steps)

    for t in 1:n_steps
        # env_state is a function that returns sampled parameters
        params = prob.env_state(t + t0 - 1)
        K = _build_param_kernel(prob, params)
        kernel_matrices[t] = K
        n_new = K * u[t]
        pop_t = sum(u[t])
        pop_new = sum(n_new)
        lambdas[t] = pop_t > 0 ? pop_new / pop_t : 0.0
        if prob.normalize && pop_new > 0
            n_new ./= pop_new
        end
        u[t + 1] = n_new
    end

    return IPMSolution(collect(t0:tf), u, kernel_matrices, nothing, :Success, lambdas)
end

# --- Eigen Analysis ---

function CommonSolve.solve(prob::IPMProblem, alg::EigenAnalysis; kwargs...)
    if _is_lagged_kernel(prob.kernel)
        return _solve_lagged_ipm_eigen(prob)
    end
    K = materialize(prob.kernel)
    ea = eigenanalysis_full(K)
    u = [copy(prob.n0)]
    return IPMSolution([prob.tspan[1]], u, K, ea, :Success, Float64[])
end

# --- Helper functions ---

function _build_dd_kernel(prob::IPMProblem, n_t, t)
    # prob.kernel should be a function: (n_t, t, p) -> AbstractIPMKernel
    kernel = prob.kernel(n_t, t, prob.p)
    return materialize(kernel)
end

function _build_kernel_set(prob::IPMProblem)
    # If kernel is a vector/tuple of AbstractIPMKernel, materialize each
    if prob.kernel isa AbstractVector || prob.kernel isa Tuple
        return [materialize(k) for k in prob.kernel]
    end
    # If kernel is already a vector of matrices
    return prob.kernel
end

function _build_param_kernel(prob::IPMProblem, params)
    # prob.kernel should be a function: params -> AbstractIPMKernel
    kernel = prob.kernel(params)
    return materialize(kernel)
end
