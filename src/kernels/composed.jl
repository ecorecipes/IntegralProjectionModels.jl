"""
Composed and mega-kernel types for combining sub-kernels.
"""

"""
    ComposedKernel{T<:Tuple}

Sum of sub-kernels: `K(z', z) = K₁(z', z) + K₂(z', z) + ...`.
Created by adding sub-kernels with `+`.
"""
struct ComposedKernel{T <: Tuple} <: AbstractIPMKernel
    subkernels::T
end

Base.:+(k1::AbstractSubKernel, k2::AbstractSubKernel) = ComposedKernel((k1, k2))
Base.:+(k1::ComposedKernel, k2::AbstractSubKernel) = ComposedKernel((k1.subkernels..., k2))
Base.:+(k1::AbstractSubKernel, k2::ComposedKernel) = ComposedKernel((k1, k2.subkernels...))
function Base.:+(k1::ComposedKernel, k2::ComposedKernel)
    ComposedKernel((k1.subkernels..., k2.subkernels...))
end

"""
    MegaKernel{K,S}

Block-structured kernel for general (multi-state) models.
Maps `(state_from, state_to)` pairs to sub-kernels, with states
described by a collection of continuous and discrete domains.

# Fields
- `kernels`: Dict mapping `(Symbol, Symbol)` pairs to kernels
- `states`: NamedTuple of domains (ContinuousDomain or DiscreteDomain)
"""
struct MegaKernel{K, S}
    kernels::K
    states::S
end

function MegaKernel(; states, kernels...)
    kernel_dict = Dict{Tuple{Symbol, Symbol}, AbstractIPMKernel}()
    for (key, kern) in kernels
        # Parse key like :size_to_size or :seedbank_to_size
        parts = split(String(key), "_to_")
        length(parts) == 2 ||
            throw(ArgumentError("Kernel key must be in format :from_to_to, got :$key"))
        from = Symbol(parts[1])
        to = Symbol(parts[2])
        kernel_dict[(from, to)] = kern
    end
    MegaKernel(kernel_dict, states)
end
