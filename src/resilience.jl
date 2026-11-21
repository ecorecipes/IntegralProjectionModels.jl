# Resilience operators: `ResilienceCore.linearise` for integral projection models.

# The mesh weights of a kernel's states. A MegaKernel stacks its states in the order of
# `keys(states)`, as `materialize` does.
_kernel_weights(k::MegaKernel, domain) = reduce(vcat, [ResilienceCore.mesh_weights(k.states[s]) for s in keys(k.states)])
_kernel_weights(k, domain::StructuredPopulationCore.AbstractStateDomain) = ResilienceCore.mesh_weights(domain)
_kernel_weights(k, domain) =
    throw(ArgumentError("cannot find the mesh of an IPMProblem with domain $(typeof(domain))"))

"""
    linearise(prob::IPMProblem) -> ResilienceCore.ProjectionMatrix

The projection operator of a density-independent, deterministic integral projection model, for
the stability and transient metrics of [ResilienceCore.jl](https://github.com/ecorecipes/ResilienceCore.jl)
(`growth_rate`, `resilience`, `reactivity`, `max_amplification`, …).

The kernel is discretised with `materialize`, and the operator carries the mesh weights of its
domain (`ResilienceCore.mesh_weights`: the bin width for each mesh point, one for each discrete
state), so its metrics approximate those of the integral operator and converge as the mesh is
refined.

Only such models have a single kernel. A density-dependent `IPMProblem` is linearised at a fixed
point instead, with `linearise(to_discrete_problem(prob), n_star)`; a stochastic one has a
stochastic growth rate but no single spectrum. Both throw an `ArgumentError`.
"""
function ResilienceCore.linearise(prob::IPMProblem)
    prob.density isa StructuredPopulationCore.DensityIndependent ||
        throw(ArgumentError("a density-dependent IPM has no single kernel; linearise its map at a fixed point " *
                            "instead, with `linearise(to_discrete_problem(prob), n_star)`"))
    prob.stochasticity isa StructuredPopulationCore.Deterministic ||
        throw(ArgumentError("a stochastic IPM has no single kernel; use its stochastic growth rate"))
    K = Matrix{Float64}(materialize(prob.kernel))
    w = _kernel_weights(prob.kernel, prob.domain)
    length(w) == size(K, 1) ||
        throw(ArgumentError("the kernel has $(size(K, 1)) states but the domain $(length(w)) mesh points"))
    return ResilienceCore.ProjectionMatrix(K; weights = w)
end
