# Getting Started

## What is an IPM?

An Integral Projection Model describes how a population structured by a continuous state variable (e.g., body size) changes over discrete time steps. The core equation is:

```math
n(z', t+1) = \int_L^U K(z', z) \, n(z, t) \, dz
```

where `K(z', z)` is the **kernel** describing transitions from size `z` to `z'`, and `n(z, t)` is the size distribution at time `t`.

## Building Your First IPM

### 1. Define the Domain

```julia
using IntegralProjectionModels, Distributions

domain = ContinuousDomain(0.0, 50.0, 100)  # [0, 50] with 100 meshpoints
```

### 2. Define Vital Rates

```julia
survival = LinearSurvival(2.2, 0.25)        # s(z) = logistic(2.2 + 0.25z)
growth = NormalGrowth(0.2, 1.02, 0.7)       # g(z'|z) = N(0.2 + 1.02z, 0.7)
fecundity = FecundityRate(0.003, 0.015, 2.0, 0.3, 1.0)
```

### 3. Build Kernels

```julia
P = PKernel(survival, growth, domain)   # survival-growth kernel
F = FKernel(fecundity, domain)          # fecundity kernel
kernel = P + F                           # composed kernel
```

### 4. Create and Solve the Problem

```julia
n0 = uniform_population(domain)
prob = IPMProblem(kernel, domain, n0, (0, 100))
sol = solve(prob)
```

### 5. Analyze Results

```julia
lambda(sol)                  # asymptotic growth rate
stable_distribution(sol)     # stable size distribution
sensitivity(sol)             # sensitivity matrix
elasticity(sol)              # elasticity matrix
```
