# IntegralProjectionModels.jl

A Julia package for building and analyzing Integral Projection Models (IPMs).

## Overview

Integral Projection Models are mathematical tools for studying structured population dynamics where individuals are characterized by a continuous state variable (e.g., body size, biomass). This package provides:

- **Typed vital rates** (survival, growth, fecundity) as callable structs
- **Composable kernels** with midpoint rule discretization
- **SciML-compatible** Problem/solve interface
- **All 12 model types** from ipmr (simple/general × DI/DD × det/stoch)
- **AD-compatible** via ForwardDiff.jl
- **Bayesian inference** extensions for Turing.jl and Gen.jl

## Quick Start

```julia
using IntegralProjectionModels, Distributions

# Define domain and vital rates
domain = ContinuousDomain(0.0, 50.0, 100)
survival = LinearSurvival(2.2, 0.25)
growth = NormalGrowth(0.2, 1.02, 0.7)
fecundity = FecundityRate(0.003, 0.015, 2.0, 0.3, 1.0)

# Build and solve
P = PKernel(survival, growth, domain)
F = FKernel(fecundity, domain)
sol = solve(IPMProblem(P + F, domain, uniform_population(domain), (0, 100)))

# Analyze
println("λ = ", lambda(sol))
```

## API Reference

```@autodocs
Modules = [IntegralProjectionModels]
```
