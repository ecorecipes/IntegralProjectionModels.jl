"""
Abstract type hierarchy for IPM kernels.

Kernels represent transition functions between states. They are
materialized (discretized) into matrices via the midpoint rule.
"""
abstract type AbstractIPMKernel end
abstract type AbstractSubKernel <: AbstractIPMKernel end

"""
Kernel family classification for general (multi-state) models.
"""
@enum KernelFamily begin
    CC  # continuous → continuous
    CD  # continuous → discrete
    DC  # discrete → continuous
    DD  # discrete → discrete
end
