# Types & Traits

Trait types classify IPMs along several axes: simple vs. general structure, density dependence, stochasticity, and projection method.

## IPM Structure

```@docs
AbstractIPMStructure
SimpleIPM
GeneralIPM
```

## Projection Structure

```@docs
AbstractProjectionStructure
DirectIteration
EigenAnalysis
```

## Density Dependence

```@docs
AbstractDensityDependence
DensityIndependent
DensityDependent
```

## Stochasticity

```@docs
AbstractStochasticity
Deterministic
StochasticKernelResampled
StochasticParameterResampled
```
