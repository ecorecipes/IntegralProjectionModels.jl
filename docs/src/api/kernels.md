# Kernels

Kernels are the core of an IPM, encoding how individuals transition between states over one time step. This page covers kernel types, eviction correction, subkernels, and kernel composition.

## Abstract Types

```@docs
AbstractIPMKernel
AbstractSubKernel
```

## Kernel Families

Kernel families classify transitions between continuous (C) and discrete (D) state variables.

```@docs
KernelFamily
CC
CD
DC
DD
```

## Eviction Correction

Eviction correction prevents loss of individuals outside domain boundaries due to truncated distributions.

```@docs
EvictionCorrection
NoCorrection
TruncatedDistributions
DiscreteExtrema
truncated_growth
apply_discrete_extrema!
```

## Subkernels

Subkernels represent individual demographic processes (survival/growth, fecundity, etc.).

```@docs
PKernel
FKernel
CustomKernel
MatrixKernel
```

## Composed Kernels

Composed kernels combine subkernels into a full projection kernel.

```@docs
ComposedKernel
MegaKernel
materialize
```
