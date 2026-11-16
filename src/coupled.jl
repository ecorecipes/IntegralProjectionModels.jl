"""
Coupled integral projection models: several populations, each an IPM, stepped together.

A `CoupledIPMProblem` is a `PopulationSystem` (the coupled-system data layer of
StructuredPopulationCore) whose components are stepped by IPM kernels. It runs on the
stepping engine that StructuredPopulationCore shares with MatrixProjectionModels.jl: each
time step the scheduled events fire, the substeps run, every component is stepped by its
kernel for that step, and the rules are applied. What is specific to IPMs is here: a kernel
that is already discretised (`DiscretizedKernel`), how a component's kernel is found and
discretised each step, and a demographic mode in which individuals are whole.
"""

# ---------------------------------------------------------------------------
# A kernel that is already discretised
# ---------------------------------------------------------------------------

"""
    DiscretizedKernel(survival, fecundity, states)

A kernel given by its discretisation: the matrices of survival and growth (`survival`) and
of births (`fecundity`) on a layout of states, as a model that builds them on its own mesh
(or a projection net) would produce them. `states` is a `NamedTuple` of domains in the
order of the rows and columns, as for a `MegaKernel`: a `ContinuousDomain` takes
as many rows as it has mesh points, a `DiscreteDomain` one per level.

`materialize(k)` is `survival .+ fecundity`. Keeping the two apart lets the demographic
mode draw survival as a multinomial and births as Poisson numbers.
"""
struct DiscretizedKernel{M <: AbstractMatrix, S <: NamedTuple} <: AbstractIPMKernel
    survival::M
    fecundity::M
    states::S

    function DiscretizedKernel(survival::M, fecundity::M, states::S) where {M <: AbstractMatrix, S <: NamedTuple}
        k = size(survival, 1)
        size(survival) == (k, k) || throw(DimensionMismatch("the survival matrix must be square"))
        size(fecundity) == (k, k) || throw(DimensionMismatch(
            "the fecundity matrix has size $(size(fecundity)); expected ($k, $k)"))
        all(d -> d isa StructuredPopulationCore.AbstractStateDomain, values(states)) ||
            throw(ArgumentError("states must be a NamedTuple of ContinuousDomains and DiscreteDomains"))
        total = isempty(states) ? 0 : sum(n_states, values(states))
        total == k || throw(DimensionMismatch("the states have $total classes in all; the matrices have $k rows"))
        return new{M, S}(survival, fecundity, states)
    end
end

function DiscretizedKernel(survival::AbstractMatrix, fecundity::AbstractMatrix, states::NamedTuple)
    T = promote_type(eltype(survival), eltype(fecundity))
    return DiscretizedKernel(Matrix{T}(survival), Matrix{T}(fecundity), states)
end

"""
    materialize(k::DiscretizedKernel) -> Matrix

The projection matrix of a discretised kernel: `k.survival .+ k.fecundity`.
"""
materialize(k::DiscretizedKernel) = k.survival .+ k.fecundity

# the names of the classes of a layout: the mesh cells of a continuous state are numbered
_class_names(name::Symbol, d::ContinuousDomain) = [Symbol(name, i) for i in 1:d.n_meshpoints]
_class_names(name::Symbol, d::DiscreteDomain) = copy(d.labels)

function component_stage_names(k::DiscretizedKernel, n::Int)
    validate_component_model(k, n)
    return reduce(vcat, [_class_names(name, d) for (name, d) in pairs(k.states)]; init = Symbol[])
end

function validate_component_model(k::DiscretizedKernel, n::Int)
    size(k.survival) == (n, n) || throw(DimensionMismatch(
        "kernel has size $(size(k.survival)); expected ($n, $n)"))
end

function Base.show(io::IO, k::DiscretizedKernel)
    print(io, "DiscretizedKernel(", size(k.survival, 1), " classes; states = ", join(keys(k.states), ", "), ")")
end

# A single-population demographic solve can split it too.
_survival_fecundity_matrices(k::DiscretizedKernel) = (Matrix{Float64}(k.survival), Matrix{Float64}(k.fecundity))

# ---------------------------------------------------------------------------
# Problem and solution
# ---------------------------------------------------------------------------

"""
    CoupledIPMProblem(system, tspan; p = nothing, substeps = [], rules = [], events = [],
                      observables = [], normalize = false, inflows = Dict(),
                      stochasticity = Deterministic())

Several populations, each an IPM, stepped together. `system` is a `PopulationSystem` whose
components carry their numbers per class in a vector. The model of a component is a kernel
(`PKernel`, `FKernel`, `CustomKernel`, `MatrixKernel`, a sum of them, a `MegaKernel` or a
[`DiscretizedKernel`](@ref)), a matrix, or a function `(system, day, p)` returning one
of these, which is called every time step: this is how a kernel comes to depend on the other
populations, through the system's state.

Each time step, as for a `CoupledMPMProblem` of MatrixProjectionModels.jl: the events that are
due fire; the substeps run in order; every component is stepped by its kernel for the step,
discretised by `materialize`, all from the same state; the expected `inflows` are added; the
rules are applied together; and with `normalize` each component is scaled to sum to one.

`inflows` maps a component's name to a function `(system, day, p)` returning the expected
number of individuals that arrive in each of its classes in a step (immigration, a seed
rain), evaluated with the state the substeps left.

With `stochasticity = Demographic()`, or an `rng` given to `solve`, the step has
demographic stochasticity: see `solve`.
"""
struct CoupledIPMProblem{Sys, P, S, R, E, O, I, T <: AbstractStochasticity}
    system::Sys
    tspan::Tuple{Int, Int}
    p::P
    substeps::S
    rules::R
    events::E
    observables::O
    normalize::Bool
    inflows::I
    stochasticity::T
end

function CoupledIPMProblem(system::PopulationSystem, tspan::Tuple{Int, Int};
        p = nothing,
        substeps::AbstractVector = AbstractDailySubstep[],
        rules::AbstractVector = AbstractTransitionRule[],
        events::AbstractVector = AbstractScheduledEvent[],
        observables::AbstractVector = Observable[],
        normalize::Bool = false,
        inflows::AbstractDict = Dict{Symbol, Any}(),
        stochasticity::AbstractStochasticity = Deterministic())
    tspan[2] >= tspan[1] || throw(ArgumentError("tspan must satisfy tspan[2] >= tspan[1]"))
    stochasticity isa Union{Deterministic, Demographic} || throw(ArgumentError(
        "a CoupledIPMProblem is Deterministic() or Demographic(), got $(typeof(stochasticity))"))
    _ensure_unique_names(substeps, "substep")
    _ensure_unique_names(rules, "rule")
    _ensure_unique_names(observables, "observable")
    _require_vector_carriers(system, "CoupledIPMProblem")
    for name in keys(inflows)
        haskey(system, name) || throw(ArgumentError("inflows into :$name, which is not a component of the system"))
    end
    return CoupledIPMProblem(system, tspan, p, collect(substeps), collect(rules), collect(events),
        collect(observables), normalize, Dict{Symbol, Any}(pairs(inflows)), stochasticity)
end

"""
    remake(prob::CoupledIPMProblem; kwargs...)

A `CoupledIPMProblem` with selected fields replaced.
"""
function remake(prob::CoupledIPMProblem;
        system = prob.system,
        tspan = prob.tspan,
        p = prob.p,
        substeps = prob.substeps,
        rules = prob.rules,
        events = prob.events,
        observables = prob.observables,
        normalize = prob.normalize,
        inflows = prob.inflows,
        stochasticity = prob.stochasticity)
    CoupledIPMProblem(system, tspan; p, substeps, rules, events, observables, normalize, inflows, stochasticity)
end

function Base.show(io::IO, prob::CoupledIPMProblem)
    print(io,
        "CoupledIPMProblem($(length(prob.system)) components, ",
        "substeps=$(length(prob.substeps)), rules=$(length(prob.rules)), ",
        "events=$(length(prob.events)), ",
        prob.stochasticity isa Demographic ? "demographic, " : "",
        "tspan=$(prob.tspan))")
end

"""
    CoupledIPMSolution

Result of solving a [`CoupledIPMProblem`](@ref), with the fields of a coupled solution of
MatrixProjectionModels.jl: `t`, `u` (for each time point, the classes of every component by
name), `component_names`, `component_matrices` (`component_matrices[name][step]` is the
discretised kernel used at a step; empty with `save_matrices = false`), `component_totals`,
`substep_log`, `observables`, `event_log`, `rule_log` and `retcode`. `sol[name]` is the
total of a component at each time point.
"""
struct CoupledIPMSolution{T, U, M, S, O, R} <: AbstractProjectionSolution
    t::Vector{Int}
    u::U
    component_names::Vector{Symbol}
    component_matrices::M
    component_totals::Dict{Symbol, Vector{T}}
    substep_log::S
    observables::O
    event_log::Vector{Tuple{Int, Symbol}}
    rule_log::R
    retcode::Symbol
end

Base.getindex(sol::CoupledIPMSolution, name::Symbol) = sol.component_totals[name]

function Base.show(io::IO, sol::CoupledIPMSolution)
    print(io,
        "CoupledIPMSolution($(length(sol.component_names)) components, ",
        "$(length(sol.t)) time points, retcode=$(sol.retcode))")
end

# ---------------------------------------------------------------------------
# Stepping a component
# ---------------------------------------------------------------------------

const _IPMOperator = Union{AbstractIPMKernel, MegaKernel, AbstractMatrix}

# The kernel of a component for this step: its model, or what its model returns.
function _component_kernel(comp::PopulationComponent, sys::PopulationSystem, day::Int, p)
    model = comp.model
    kernel = if model isa _IPMOperator
        model
    elseif applicable(model, sys, day, p)
        model(sys, day, p)
    else
        throw(ArgumentError(
            "component model $(typeof(model)) is neither a kernel, a matrix nor callable as (system, day, p)"))
    end
    kernel isa _IPMOperator || throw(ArgumentError(
        "a component model must give an IPM kernel, a MegaKernel or a matrix, got $(typeof(kernel))"))
    return kernel
end

_kernel_matrix(k::AbstractMatrix) = k
_kernel_matrix(k) = materialize(k)

function _check_size(matrix::AbstractMatrix, n::Int)
    size(matrix) == (n, n) || throw(DimensionMismatch(
        "component matrix has size $(size(matrix)); expected ($n, $n)"))
    return matrix
end

# Deterministic: the discretised kernel, applied to the stage vector.
_discretised_kernel(comp::PopulationComponent, sys::PopulationSystem, day::Int, p) =
    _check_size(_kernel_matrix(_component_kernel(comp, sys, day, p)), length(comp.population))

function _with_inflows(inflows)
    return function (name::Symbol, A, n::AbstractVector, sys::PopulationSystem, day::Int, p)
        stepped = Vector(A * n)
        haskey(inflows, name) && (stepped .+= inflows[name](sys, day, p))
        return stepped
    end
end

# Demographic: survival and growth drawn as a multinomial down each column, births as Poisson
# numbers, then the inflows as Poisson numbers, for one component after the other.
_survival_fecundity(k::DiscretizedKernel) = (k.survival, k.fecundity)
_survival_fecundity(k::Union{PKernel, FKernel, ComposedKernel}) = _survival_fecundity_matrices(k)
_survival_fecundity(k) = throw(ArgumentError(
    "the demographic mode needs survival and births apart: give a component a DiscretizedKernel, or " *
    "PKernel and FKernel sub-kernels, not a $(typeof(k))"))

function _demographic_step(rng::AbstractRNG, inflows)
    return function (name::Symbol, kernel, n::AbstractVector, sys::PopulationSystem, day::Int, p)
        all(isinteger, n) || throw(ArgumentError(
            "demographic stochasticity needs whole numbers of individuals in :$name"))
        U, F = _survival_fecundity(kernel)
        _check_size(U, length(n))
        _check_size(F, length(n))
        worst = maximum(sum(U; dims = 1))
        worst <= 1 + 1e-6 || throw(ArgumentError(
            "survival matrix of :$name has a column summing to $worst; demographic stochasticity needs " *
            "survival probabilities (column sums at most 1)"))
        counts = demographic_step!(rng, zeros(Int, length(n)), round.(Int, n), U, F)
        stepped = Float64.(counts)
        if haskey(inflows, name)
            for (i, r) in enumerate(inflows[name](sys, day, p))
                r > 0 && (stepped[i] += rand_poisson(rng, r))
            end
        end
        return stepped
    end
end

_saved_kernel_matrix(k) = Matrix(_kernel_matrix(k))

# ---------------------------------------------------------------------------
# Solve
# ---------------------------------------------------------------------------

function CommonSolve.solve(prob::CoupledIPMProblem, ::EigenAnalysis; kwargs...)
    throw(ArgumentError("EigenAnalysis is not defined for CoupledIPMProblem; freeze a component as an IPMProblem"))
end

"""
    solve(prob::CoupledIPMProblem, DirectIteration(); save_matrices = true, rng = nothing)

Iterate a coupled IPM one time step at a time and return a [`CoupledIPMSolution`](@ref).
With `save_matrices = false` the discretised kernel of each component and step is not kept.

The problem is solved with demographic stochasticity if its `stochasticity` is
`Demographic()` or an `rng` is given (with `Demographic()` and no `rng`, the default random
number generator is used). Individuals are then whole: the numbers in every class must be
whole numbers, and each component's kernel must keep survival and births apart (a
`DiscretizedKernel`, or `PKernel` and `FKernel` sub-kernels), with survival probabilities
that sum to at most one down each column. Each step, for one component after the other, in
the order of the system:

1. the individuals of each class survive and move according to a multinomial draw from
   their column of the survival matrix, and give births drawn as Poisson numbers with means
   from the fecundity matrix (`demographic_step!` of StructuredPopulationCore);
2. the component's expected inflows, if any, are drawn as Poisson numbers and added.

The rules are then applied to the drawn state as in the deterministic mode. Given the same
`rng` in the same state, the draws are the same.
"""
function CommonSolve.solve(prob::CoupledIPMProblem, ::DirectIteration = DirectIteration();
        save_matrices::Bool = true, rng = nothing, kwargs...)
    demographic = prob.stochasticity isa Demographic || rng !== nothing
    r = if demographic
        generator = rng === nothing ? Random.default_rng() : rng
        solve_coupled(_component_kernel, prob.system, prob.tspan, prob.p,
            prob.substeps, prob.rules, prob.events, prob.observables, prob.normalize;
            save_matrices, step_population = _demographic_step(generator, prob.inflows),
            saved_matrix = _saved_kernel_matrix)
    else
        step = isempty(prob.inflows) ? _project_population : _with_inflows(prob.inflows)
        solve_coupled(_discretised_kernel, prob.system, prob.tspan, prob.p,
            prob.substeps, prob.rules, prob.events, prob.observables, prob.normalize;
            save_matrices, step_population = step)
    end
    return CoupledIPMSolution(r.t, r.u, r.component_names, r.component_matrices, r.component_totals,
        r.substep_log, r.observables, r.event_log, r.rule_log, r.retcode)
end
