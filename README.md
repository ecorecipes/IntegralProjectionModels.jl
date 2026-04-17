# IntegralProjectionModels.jl

A Julia package for building and analyzing Integral Projection Models (IPMs) — mathematical models for structured population dynamics based on continuous state variables.

## Features

- **All 12 ipmr model types**: simple/general × density-independent/dependent × deterministic/stochastic (kernel-resampled, parameter-resampled)
- **SciML-compatible**: `IPMProblem`/`solve()` pattern, conversion to `DiscreteProblem`
- **AD-compatible**: all vital rate types work with ForwardDiff.jl
- **Age×size models**: age-structured expansion of size-based kernels
- **Eviction correction**: truncated distribution and discrete extrema methods
- **Analysis**: asymptotic growth rate (λ), sensitivity, elasticity, stochastic growth rate
- **Extensible**: package extensions for Turing.jl, Gen.jl, and ModelingToolkit.jl

## Quick Start

```julia
using IntegralProjectionModels
using Distributions

# Define domain
domain = ContinuousDomain(0.0, 50.0, 100)

# Define vital rates
survival = LinearSurvival(2.2, 0.25)
growth = NormalGrowth(0.2, 1.02, 0.7)
fecundity = FecundityRate(0.003, 0.015, 2.0, 0.3, 1.0)

# Build kernels
P = PKernel(survival, growth, domain)
F = FKernel(fecundity, domain)

# Create and solve problem
n0 = uniform_population(domain)
prob = IPMProblem(P + F, domain, n0, (0, 100))
sol = solve(prob)

# Analysis
println("λ = ", lambda(sol))
E = elasticity(sol)
```

## Model Types

| Structure | Density | Stochasticity | Constructor |
|-----------|---------|---------------|-------------|
| Simple | DI | Deterministic | `IPMProblem(kernel, domain, n0, tspan)` |
| Simple | DD | Deterministic | `IPMProblem(DensityDependent(), kernel_fn, domain, n0, tspan)` |
| Simple | DI | Stoch (kern) | `IPMProblem(StochasticKernelResampled(), kernels, domain, n0, tspan)` |
| Simple | DI | Stoch (param) | `IPMProblem(StochasticParameterResampled(), builder, domain, n0, tspan; env_state=sampler)` |
| General | DI | Deterministic | `IPMProblem(GeneralIPM(), mega_kernel, states, n0, tspan)` |
| ... | ... | ... | All 12 combinations supported |

## Vignettes

| # | Vignette | Description |
|---|----------|-------------|
| 1 | [Introduction to Integral Projection Models](https://github.com/ecorecipes/IntegralProjectionModels.jl/blob/main/vignettes/01_introduction/01_introduction.md) | Core concepts: kernels, vital rates, domain discretization, eigenanalysis |
| 2 | [Ungulate IPM: Soay Sheep](https://github.com/ecorecipes/IntegralProjectionModels.jl/blob/main/vignettes/02_ungulate/02_ungulate.md) | Complete worked example fitting an IPM to Soay sheep data |
| 3 | [Eviction Correction](https://github.com/ecorecipes/IntegralProjectionModels.jl/blob/main/vignettes/03_eviction/03_eviction.md) | Correcting for individuals lost outside domain bounds |
| 4 | [Sensitivity and Elasticity Analysis](https://github.com/ecorecipes/IntegralProjectionModels.jl/blob/main/vignettes/04_sensitivity_elasticity/04_sensitivity_elasticity.md) | Perturbation analysis of continuous-state models |
| 5 | [Density-Dependent IPMs](https://github.com/ecorecipes/IntegralProjectionModels.jl/blob/main/vignettes/05_density_dependent/05_density_dependent.md) | Models with density-dependent vital rates |
| 6 | [Stochastic IPMs](https://github.com/ecorecipes/IntegralProjectionModels.jl/blob/main/vignettes/06_stochastic/06_stochastic.md) | Kernel-resampled and parameter-resampled stochastic models |
| 7 | [General IPM with Seed Bank](https://github.com/ecorecipes/IntegralProjectionModels.jl/blob/main/vignettes/07_general_ipm/07_general_ipm.md) | Multi-state models mixing continuous and discrete states |
| 8 | [Age x Size IPM](https://github.com/ecorecipes/IntegralProjectionModels.jl/blob/main/vignettes/08_age_size/08_age_size.md) | Age-structured expansion of size-based kernels |
| 9 | [Categorical Composition of IPMs](https://github.com/ecorecipes/IntegralProjectionModels.jl/blob/main/vignettes/09_categorical/09_categorical.md) | Composing IPMs using categorical/functorial methods |
| 10 | [PADRINO Database Integration](https://github.com/ecorecipes/IntegralProjectionModels.jl/blob/main/vignettes/10_padrino/10_padrino.md) | Downloading, building, and analyzing models from the PADRINO database |
| 11 | [Evolving IPMs: Evolutionary Demography](https://github.com/ecorecipes/IntegralProjectionModels.jl/blob/main/vignettes/11_evolving_ipm/11_evolving_ipm.md) | IPMs with heritable trait variation and selection |
| 12 | [Time-Lagged Integral Projection Models](https://github.com/ecorecipes/IntegralProjectionModels.jl/blob/main/vignettes/12_time_lag/12_time_lag.md) | State-augmented models with delayed fecundity (Kuss et al. 2008) |

## Installation

This package is not yet registered in the Julia General registry. Install directly from GitHub (the [ProjectionModels.jl](https://github.com/ecorecipes/StructuredPopulationCore.jl) dependency must be installed first):

```julia
using Pkg
Pkg.add(url="https://github.com/ecorecipes/StructuredPopulationCore.jl")
Pkg.add(url="https://github.com/ecorecipes/IntegralProjectionModels.jl")
```

## Related

- [StructuredPopulationCore.jl](https://github.com/ecorecipes/StructuredPopulationCore.jl) — shared abstractions
- [MatrixProjectionModels.jl](https://github.com/ecorecipes/MatrixProjectionModels.jl) — discrete-stage matrix models
- [FiniteStatePopulationDynamics.jl](https://github.com/ecorecipes/FiniteStatePopulationDynamics.jl) — discrete-state continuous-time dynamics
- [ContinuousStatePopulationDynamics.jl](https://github.com/ecorecipes/ContinuousStatePopulationDynamics.jl) — continuous-state continuous-time dynamics
- [CategoricalPopulationDynamics.jl](https://github.com/ecorecipes/CategoricalPopulationDynamics.jl) — categorical/functorial framework
- [PhysiologicallyBasedDemographicModels.jl](https://github.com/ecorecipes/PhysiologicallyBasedDemographicModels.jl) — application-level PBDM reference suite
- [ipmr](https://github.com/levisc8/ipmr) — R package for IPMs (design reference)
- [PADRINO](https://github.com/padrinoDB/Rpadrino) — IPM database
