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

## Installation

```julia
using Pkg
Pkg.add("IntegralProjectionModels")
```

## Related

- [ProjectionModels.jl](https://github.com/ecorecipes/ProjectionModels.jl) — shared abstractions
- [MatrixProjectionModels.jl](https://github.com/ecorecipes/MatrixProjectionModels.jl) — discrete-stage matrix models
- [CategoricalProjectionModels.jl](https://github.com/ecorecipes/CategoricalProjectionModels.jl) — categorical/functorial framework
- [ipmr](https://github.com/levisc8/ipmr) — R package for IPMs (design reference)
- [PADRINO](https://github.com/padrinoDB/Rpadrino) — IPM database
