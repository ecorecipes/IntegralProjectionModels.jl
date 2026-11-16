# Problem & Solution

The `IPMProblem` and `IPMSolution` types follow the SciML problem/solve convention. Construct a problem, solve it, then extract results.

```@docs
IPMProblem
remake
IPMSolution
solve
```

## Coupled IPMs

Several populations, each an IPM, stepped together on the coupled-system engine that
StructuredPopulationCore shares with MatrixProjectionModels.jl: scheduled events, substeps,
rules and observables act on a `PopulationSystem`, and each component is stepped by its
kernel for the step. A deterministic or a demographic run.

```@docs
CoupledIPMProblem
CoupledIPMSolution
```
