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
    kernel_set = prob.density isa DensityIndependent ?
                 _build_kernel_set(prob, prob.n0, t0) : nothing
    n_kernels = length(kernel_set === nothing ? prob.kernel : kernel_set)

    if kernel_seq === nothing
        kernel_seq = rand(1:n_kernels, n_steps)
    end

    u = Vector{typeof(prob.n0)}(undef, n_steps + 1)
    u[1] = copy(prob.n0)
    lambdas = Vector{Float64}(undef, n_steps)
    used_kernels = Vector{AbstractMatrix}(undef, n_steps)

    for t in 1:n_steps
        if prob.density isa DensityDependent
            kernel_set = _build_kernel_set(prob, u[t], t + t0 - 1)
        end
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
        params = _sample_env_state(prob, u[t], t + t0 - 1)
        K = _build_param_kernel(prob, u[t], t + t0 - 1, params)
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

# --- Demographic stochasticity (finite-population integer counts on the mesh) ---

# Split a kernel into its survival-growth (sub-stochastic, Multinomial) and
# fecundity (rate, Poisson) discretized matrices. Requires PKernel/FKernel
# sub-kernels so the two can be sampled differently.
function _survival_fecundity_matrices(kernel)
    subs = kernel isa ComposedKernel ? kernel.subkernels : (kernel,)
    P = nothing
    F = nothing
    for sub in subs
        if sub isa PKernel
            M = Matrix{Float64}(materialize(sub))
            P = P === nothing ? M : P .+ M
        elseif sub isa FKernel
            M = Matrix{Float64}(materialize(sub))
            F = F === nothing ? M : F .+ M
        else
            error("Demographic IPM solve requires PKernel (survival/growth) and " *
                  "FKernel (fecundity) sub-kernels so survival can be sampled as a " *
                  "Multinomial draw and fecundity as Poisson; got $(typeof(sub)). " *
                  "Build the kernel as PKernel(...) + FKernel(...).")
        end
    end
    P === nothing && F === nothing && error("kernel has no sub-kernels to materialize")
    P === nothing && (P = zeros(size(F)))
    F === nothing && (F = zeros(size(P)))
    return P, F
end

function _solve(::AbstractIPMStructure, ::DensityIndependent, ::Demographic,
        prob::IPMProblem, ::DirectIteration;
        rng::AbstractRNG = Random.default_rng(), kwargs...)
    P, F = _survival_fecundity_matrices(prob.kernel)
    for j in axes(P, 2)
        s = sum(@view P[:, j])
        s <= 1 + 1e-6 || error(
            "Demographic IPM solve requires sub-stochastic survival-growth columns; " *
            "column $j of the P kernel sums to $s > 1. Check that survival ≤ 1 and that " *
            "the growth kernel is normalized (consider an eviction correction).")
    end

    t0, tf = prob.tspan
    n_steps = tf - t0
    k = size(P, 1)

    u = Vector{Vector{Float64}}(undef, n_steps + 1)
    counts = round.(Int, prob.n0)
    u[1] = Float64.(counts)
    lambdas = Vector{Float64}(undef, n_steps)
    n_next = zeros(Int, k)

    for t in 1:n_steps
        prev = sum(counts)
        demographic_step!(rng, n_next, counts, P, F)
        counts = copy(n_next)
        newtot = sum(counts)
        lambdas[t] = prev > 0 ? newtot / prev : 0.0
        u[t + 1] = Float64.(counts)
    end

    return IPMSolution(collect(t0:tf), u, P .+ F, nothing, :Success, lambdas)
end

"""
    demographic_ensemble(prob::IPMProblem; n_reps=100, rng=Random.default_rng())

Run `n_reps` independent demographic-stochastic realizations of an IPM (binned on
the mesh) and return `(totals, sols)` where `totals` is an `(n_time × n_reps)`
matrix of total population sizes (consumable by `quasi_extinction`). If `prob` is
not already a `Demographic` problem it is `remake`d as one.
"""
function demographic_ensemble(prob::IPMProblem; n_reps::Int = 100,
        rng::AbstractRNG = Random.default_rng())
    dprob = prob.stochasticity isa Demographic ? prob :
            remake(prob; stochasticity = Demographic())
    sols = [solve(dprob, DirectIteration(); rng = rng) for _ in 1:n_reps]
    n_time = length(sols[1].t)
    totals = Matrix{Float64}(undef, n_time, n_reps)
    for (r, s) in enumerate(sols)
        @inbounds for tt in 1:n_time
            totals[tt, r] = sum(s.u[tt])
        end
    end
    return totals, sols
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

function _build_kernel_set(prob::IPMProblem, n_t, t)
    if prob.kernel isa AbstractVector || prob.kernel isa Tuple
        return [_materialize_kernel_entry(prob, k, n_t, t) for k in prob.kernel]
    end
    return prob.kernel
end

function _build_param_kernel(prob::IPMProblem, n_t, t, params)
    kernel = prob.density isa DensityDependent ?
             prob.kernel(n_t, t, params) :
             prob.kernel(params)
    return materialize(kernel)
end

function _materialize_kernel_entry(prob::IPMProblem, kernel, n_t, t)
    if prob.density isa DensityDependent && kernel isa Function
        return materialize(kernel(n_t, t, prob.p))
    end
    return materialize(kernel)
end

function _sample_env_state(prob::IPMProblem, n_t, t)
    applicable(prob.env_state, n_t, t) && return prob.env_state(n_t, t)
    return prob.env_state(t)
end

# Continuous-time solve methods are provided by ContinuousStatePopulationDynamics
# on the shared CommonSolve.solve generic. IntegralProjectionModels re-exports
# those types during the compatibility phase without redefining the methods.
