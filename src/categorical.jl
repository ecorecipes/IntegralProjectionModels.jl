"""
Categorical composition tools for Integral Projection Models.

General-purpose "lowering" functions inspired by applied category theory:

- **Left Kan extension**: kernel function → matrix (IPM → MPM discretisation)
- **Right Kan extension**: matrix → piecewise-constant kernel (MPM → IPM refinement)
- **Stratification**: local kernel + dispersal → block-structured spatial matrix
- **Coarsening**: fine matrix → coarse matrix via bin aggregation (pushforward)
- **Kernel composition**: sum sub-kernels additively (UWD evaluation without Catlab)

`left_kan_extension`, `right_kan_extension`, `stratify`, `coarsen` and `compose_from_uwd`
are defined in StructuredPopulationCore.jl, for the domain and matrix types shared by the
stack, and re-exported from here. They are the same functions that
CategoricalPopulationDynamics.jl extends, so the names can be used unqualified with both
packages loaded. `compose_kernels` is defined below.

None of these needs Catlab.jl; the Catlab extension (`IntegralProjectionModelsCatlabExt`)
adds UWD- and FinFunction-aware methods.
"""

"""
    compose_kernels(sub_kernels, domain::ContinuousDomain)

Compose sub-kernel functions additively and discretise onto the given domain.

This is the concrete evaluation of an **undirected wiring diagram** (UWD)
without requiring Catlab: each sub-kernel shares the same trait axes, and
the composition rule sums their contributions.

``K_{ij} = \\sum_k h \\cdot K_k(z_i, z_j)``

# Arguments
- `sub_kernels`: a `NamedTuple`, `Dict`, or iterable of `(name, fn)` pairs
  where each `fn` is a callable `(z_new, z) -> Real`
- `domain`: a `ContinuousDomain`

# Returns
- `K::Matrix{Float64}`: the summed, discretised transition matrix

See also: [`compose_from_uwd`](@ref) for the form that takes a Catlab wiring diagram.
"""
function compose_kernels(sub_kernels, domain::ContinuousDomain)
    z = meshpoints(domain)
    h = step_size(domain)
    m = length(z)
    K = zeros(m, m)
    for (_, kfn) in pairs(sub_kernels)
        for j in 1:m
            for i in 1:m
                K[i, j] += h * kfn(z[i], z[j])
            end
        end
    end
    return K
end
