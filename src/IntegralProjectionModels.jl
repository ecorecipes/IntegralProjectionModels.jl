module IntegralProjectionModels

using CommonSolve
using Distributions
using LinearAlgebra
using Random
using StructuredPopulationCore
using RecipesBase
using SciMLBase
using Statistics
using StatsFuns: logistic

# Resilience operators: `linearise` is ResilienceCore's function, extended in resilience.jl.
import ResilienceCore
using ResilienceCore: linearise

# Functions shared across the stack are extended here rather than redefined, so that their
# names stay usable when several packages are loaded together. The categorical operations
# and two verbs belong to StructuredPopulationCore, `remake` to SciMLBase.
import StructuredPopulationCore: left_kan_extension, right_kan_extension, stratify, coarsen,
    compose_from_uwd, demographic_ensemble, to_discrete_problem
import SciMLBase: remake
# The coupled-system data layer and its stepping engine are shared with
# MatrixProjectionModels.jl and live, unexported, in StructuredPopulationCore. They are
# imported here, extended where IPMs need it, and their public names exported below.
import StructuredPopulationCore: AbstractStateVariable, AbstractInteractionRule, AbstractScheduledEvent,
    PopulationComponent, PopulationSystem, Observable, PulseRelease, SingleDayRelease, SprayEvent,
    TransferRule, ReproductionRule, MortalityRule, CustomRule, get_state, set_state!, has_state,
    total_population, component_total, component_totals, inject!, remove_fraction!, apply_event!,
    component_stage_names, validate_component_model
import StructuredPopulationCore: AbstractDailySubstep, apply_substep!, StateUpdateSubstep, CustomSubstep,
    StageFlux, StageTransfer, RuleEffect, AbstractTransitionRule, apply_rule, CustomEvent, solve_coupled
import StructuredPopulationCore: _project_population, _ensure_unique_names, _require_vector_carriers

# Domain types (re-exported shared abstractions from StructuredPopulationCore)
export ContinuousDomain, DiscreteDomain
export meshpoints, step_size, bounds, n_states

# Vital rate abstract types
include("vital_rates/abstract.jl")
export AbstractVitalRate, AbstractSurvivalRate, AbstractGrowthRate
export AbstractFecundityRate, AbstractRecruitmentRate

# Vital rate implementations
include("vital_rates/survival.jl")
export LinearSurvival, QuadraticSurvival, ConstantSurvival

include("vital_rates/growth.jl")
export NormalGrowth, LogNormalGrowth, mean_size

include("vital_rates/fecundity.jl")
export FecundityRate, LogisticFecundityRate, RecruitmentDistribution

include("vital_rates/custom.jl")
export CustomVitalRate

# Individual-level samplers for the vital-rate types (individual-based realizations)
include("individual_samplers.jl")

# Kernel types
include("kernels/abstract.jl")
export AbstractIPMKernel, AbstractSubKernel, KernelFamily, CC, CD, DC, DD

include("kernels/eviction.jl")
export EvictionCorrection, NoCorrection, TruncatedDistributions, DiscreteExtrema
export truncated_growth, apply_discrete_extrema!

include("kernels/subkernels.jl")
export PKernel, FKernel, CustomKernel, MatrixKernel

include("kernels/composed.jl")
export ComposedKernel, MegaKernel

include("kernels/materialize.jl")
export materialize

# Trait types (re-exported shared abstractions from StructuredPopulationCore)
export AbstractIPMStructure, SimpleIPM, GeneralIPM
export AbstractContinuousStateStructure, SimpleContinuousState, GeneralContinuousState
# Re-export shared types from StructuredPopulationCore
export AbstractProjectionStructure
export AbstractDensityDependence, DensityIndependent, DensityDependent
export AbstractStochasticity, Deterministic, StochasticKernelResampled,
       StochasticParameterResampled, Demographic
export DirectIteration, EigenAnalysis

# Problem type
include("problems.jl")
export IPMProblem, remake

# Solve
include("solve.jl")
export IPMSolution, demographic_ensemble
# solve is re-exported via CommonSolve
using CommonSolve: solve
export solve

# Re-export shared analysis and eigenanalysis from StructuredPopulationCore
export AbstractProjectionSolution
export eigenanalysis_power, eigenanalysis_full
export lambda, stable_distribution, reproductive_value
export sensitivity, elasticity, damping_ratio
export stochastic_growth_rate, mean_kernel
export is_irreducible, is_primitive, is_ergodic
export area_under_curve

# Coupled IPMs: several populations stepped together, on the engine shared with
# MatrixProjectionModels.jl
include("coupled.jl")
export DiscretizedKernel, CoupledIPMProblem, CoupledIPMSolution
export PopulationComponent, PopulationSystem, get_state, set_state!, has_state
export total_population, component_total, component_totals, inject!, remove_fraction!
export AbstractDailySubstep, apply_substep!, StateUpdateSubstep, CustomSubstep
export StageFlux, StageTransfer, RuleEffect, apply_rule, apply_event!
export AbstractTransitionRule, AbstractInteractionRule, AbstractScheduledEvent
export TransferRule, ReproductionRule, MortalityRule, CustomRule
export PulseRelease, SingleDayRelease, SprayEvent, CustomEvent, Observable

# Time-lagged models
include("time_lag.jl")
export LaggedKernel, expand_lag_kernels
# Re-export from StructuredPopulationCore
export TimeLagStructure, expand_lag_matrix, extract_lag_components
export augment_population, extract_population
export net_repro_rate_lagged, generation_time_lagged

# Age×size models
include("age_size.jl")
export AgeStructure, ages, n_ages, expand_age_kernels

# Categorical composition
include("categorical.jl")
export left_kan_extension, right_kan_extension, stratify, coarsen
export compose_kernels, compose_from_uwd

# SciML interface
include("sciml_interface.jl")
export to_discrete_problem

# Resilience operators (ResilienceCore)
include("resilience.jl")
export linearise

# Plotting
include("plotting.jl")

# Utilities
include("utils.jl")
export rate_to_proportion, proportion_to_rate
export uniform_population, point_population, normal_population

# PADRINO database integration (stubs; full implementation via extension)
include("padrino.jl")
export PadrinoDB, PadrinoModel
export pdb_download, pdb_load, pdb_save
export pdb_subset, pdb_species, pdb_citations, pdb_metadata, pdb_test_targets
export pdb_make_proto_ipm, pdb_make_ipm, pdb_validate

end # module
