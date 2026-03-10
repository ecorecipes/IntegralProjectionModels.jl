"""
Concrete sub-kernel types for IPM construction.
"""

"""
    PKernel{S,G,D}

Survival-growth (P) kernel: `P(z', z) = s(z) * g(z', z)`.
Represents the probability of surviving and transitioning from size `z` to `z'`.
"""
struct PKernel{S, G, D} <: AbstractSubKernel
    survival::S
    growth::G
    domain::D
    eviction::EvictionCorrection
end

function PKernel(survival, growth, domain; eviction = NoCorrection)
    PKernel(survival, growth, domain, eviction)
end

"""
    FKernel{F,D}

Fecundity (F) kernel: `F(z', z) = f(z', z)`.
Represents the production of new individuals of size `z'` from parents of size `z`.
"""
struct FKernel{F, D} <: AbstractSubKernel
    fecundity::F
    domain::D
    eviction::EvictionCorrection
end

function FKernel(fecundity, domain; eviction = NoCorrection)
    FKernel(fecundity, domain, eviction)
end

"""
    CustomKernel{F,D}

Kernel defined by an arbitrary function `func(z', z)` or `func(z', z, params)`.
"""
struct CustomKernel{F, D} <: AbstractSubKernel
    func::F
    domain::D
    family::KernelFamily
    eviction::EvictionCorrection
end

function CustomKernel(func, domain; family = CC, eviction = NoCorrection)
    CustomKernel(func, domain, family, eviction)
end

"""
    MatrixKernel{M} <: AbstractSubKernel

Wraps a pre-computed matrix as a kernel. Useful for cross-domain transitions
in general models where the matrix dimensions don't match standard kernel
materialization (e.g., discrete → continuous transitions).

`materialize(mk::MatrixKernel)` returns the stored matrix directly.
"""
struct MatrixKernel{M} <: AbstractSubKernel
    matrix::M
end
