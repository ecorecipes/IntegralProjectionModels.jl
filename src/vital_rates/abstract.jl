"""
Abstract type hierarchy for vital rates.

All vital rate types are callable structs with plain numeric fields,
making them compatible with automatic differentiation.
"""
abstract type AbstractVitalRate end
abstract type AbstractSurvivalRate <: AbstractVitalRate end
abstract type AbstractGrowthRate <: AbstractVitalRate end
abstract type AbstractFecundityRate <: AbstractVitalRate end
abstract type AbstractRecruitmentRate <: AbstractVitalRate end

# logistic(x) imported from StatsFuns in main module
