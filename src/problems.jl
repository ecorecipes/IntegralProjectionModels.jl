"""
    IPMProblem

Central problem type parameterized by trait types. Encapsulates
all information needed to build and iterate an IPM.

# Type parameters
- `S`: AbstractIPMStructure (SimpleIPM or GeneralIPM)
- `D`: AbstractDensityDependence (DensityIndependent or DensityDependent)
- `T`: AbstractStochasticity (Deterministic, StochasticKernelResampled, StochasticParameterResampled)
- `K`: Kernel type (AbstractIPMKernel or function for DD/stoch)
- `Dom`: Domain type
- `U`: Population state type
- `P`: Parameters type
- `E`: Environment state type
"""
struct IPMProblem{S <: AbstractIPMStructure, D <: AbstractDensityDependence,
    T <: AbstractStochasticity, K, Dom, U, P, E}
    structure::S
    density::D
    stochasticity::T
    kernel::K
    domain::Dom
    n0::U
    tspan::Tuple{Int, Int}
    p::P
    env_state::E
    normalize::Bool
    uses_age::Bool
end

# Full constructor
function IPMProblem(structure::AbstractIPMStructure,
        density::AbstractDensityDependence,
        stochasticity::AbstractStochasticity,
        kernel, domain, n0, tspan;
        p = nothing, env_state = nothing,
        normalize = false, uses_age = false)
    IPMProblem(structure, density, stochasticity,
        kernel, domain, n0, tspan, p, env_state, normalize, uses_age)
end

# Convenience: simple_di_det (most common case)
function IPMProblem(kernel, domain, n0, tspan; kwargs...)
    IPMProblem(SimpleIPM(), DensityIndependent(), Deterministic(),
        kernel, domain, n0, tspan; kwargs...)
end

# Convenience: specify density dependence
function IPMProblem(density::DensityDependent, kernel, domain, n0, tspan; kwargs...)
    IPMProblem(SimpleIPM(), density, Deterministic(),
        kernel, domain, n0, tspan; kwargs...)
end

# Convenience: specify stochasticity
function IPMProblem(stoch::AbstractStochasticity, kernel, domain, n0, tspan; kwargs...)
    IPMProblem(SimpleIPM(), DensityIndependent(), stoch,
        kernel, domain, n0, tspan; kwargs...)
end

# Convenience: specify density + stochasticity
function IPMProblem(density::AbstractDensityDependence,
        stoch::AbstractStochasticity,
        kernel, domain, n0, tspan; kwargs...)
    IPMProblem(SimpleIPM(), density, stoch,
        kernel, domain, n0, tspan; kwargs...)
end

# Convenience: general model
function IPMProblem(structure::GeneralIPM, kernel, domain, n0, tspan; kwargs...)
    IPMProblem(structure, DensityIndependent(), Deterministic(),
        kernel, domain, n0, tspan; kwargs...)
end

"""
    remake(prob::IPMProblem; kwargs...)

Create a new IPMProblem with selected fields replaced.
"""
function remake(prob::IPMProblem;
        structure = prob.structure,
        density = prob.density,
        stochasticity = prob.stochasticity,
        kernel = prob.kernel,
        domain = prob.domain,
        n0 = prob.n0,
        tspan = prob.tspan,
        p = prob.p,
        env_state = prob.env_state,
        normalize = prob.normalize,
        uses_age = prob.uses_age)
    IPMProblem(structure, density, stochasticity,
        kernel, domain, n0, tspan, p, env_state, normalize, uses_age)
end
