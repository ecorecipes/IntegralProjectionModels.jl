module IntegralProjectionModels

using CommonSolve
using Distributions
using LinearAlgebra
using ProjectionModels
using RecipesBase
using SciMLBase
using Statistics
using StatsFuns: logistic

# Domain types
include("domains.jl")
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

# Trait types
include("types.jl")
export AbstractIPMStructure, SimpleIPM, GeneralIPM
# Re-export shared types from ProjectionModels
export AbstractProjectionStructure
export AbstractDensityDependence, DensityIndependent, DensityDependent
export AbstractStochasticity, Deterministic, StochasticKernelResampled,
       StochasticParameterResampled
export DirectIteration, EigenAnalysis

# Problem type
include("problems.jl")
export IPMProblem, remake

# Solve
include("solve.jl")
export IPMSolution
# solve is re-exported via CommonSolve
using CommonSolve: solve
export solve

# Re-export shared analysis and eigenanalysis from ProjectionModels
export AbstractProjectionSolution
export eigenanalysis_power, eigenanalysis_full
export lambda, stable_distribution, reproductive_value
export sensitivity, elasticity, damping_ratio
export stochastic_growth_rate, mean_kernel
export is_irreducible, is_primitive, is_ergodic
export area_under_curve

# Time-lagged models
include("time_lag.jl")
export LaggedKernel, expand_lag_kernels
# Re-export from ProjectionModels
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
