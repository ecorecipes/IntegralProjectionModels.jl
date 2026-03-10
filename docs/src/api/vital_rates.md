# Vital Rates

Vital rates are the building blocks of IPM kernels. Each vital rate is a callable struct that maps individual state (e.g., body size) to a demographic rate or distribution.

## Abstract Types

```@docs
AbstractVitalRate
AbstractSurvivalRate
AbstractGrowthRate
AbstractFecundityRate
AbstractRecruitmentRate
```

## Survival Rates

Survival rates map individual state to a probability of surviving to the next time step.

```@docs
LinearSurvival
QuadraticSurvival
ConstantSurvival
```

## Growth Rates

Growth rates describe the distribution of state transitions for surviving individuals.

```@docs
NormalGrowth
LogNormalGrowth
mean_size
```

## Fecundity Rates

Fecundity rates describe reproduction: the probability of reproducing and the distribution of offspring states.

```@docs
FecundityRate
LogisticFecundityRate
RecruitmentDistribution
```

## Custom Vital Rates

For user-defined vital rate functions that do not fit the built-in types.

```@docs
CustomVitalRate
```
