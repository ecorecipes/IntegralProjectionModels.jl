"""
Trait types for IPM classification.

Every IPM is classified by 3 traits (matching ipmr's 4-parameter system,
with kern/param folded into stochasticity):

1. Structure: SimpleIPM (single continuous state) vs GeneralIPM (multiple states)
2. Density dependence: DensityIndependent vs DensityDependent (from ProjectionModels)
3. Stochasticity: Deterministic, StochasticKernelResampled, StochasticParameterResampled (from ProjectionModels)
"""

# Structure traits (IPM-specific)
abstract type AbstractIPMStructure <: AbstractProjectionStructure end
struct SimpleIPM <: AbstractIPMStructure end
struct GeneralIPM <: AbstractIPMStructure end

# Density dependence, stochasticity, and algorithm types are
# imported from ProjectionModels via `using ProjectionModels`
