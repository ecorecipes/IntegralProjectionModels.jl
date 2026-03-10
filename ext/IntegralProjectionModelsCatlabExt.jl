module IntegralProjectionModelsCatlabExt

using IntegralProjectionModels
using Catlab
using Catlab.CategoricalAlgebra  # FinFunction, FinSet
using Catlab.WiringDiagrams      # UWD boxes
using Catlab.Programs             # @relation

"""
    compose_from_uwd(uwd, sub_kernels::Dict, domain::ContinuousDomain)

Evaluate a Catlab undirected wiring diagram with concrete sub-kernel functions.

Each box in the UWD is looked up by name in `sub_kernels`; the composition
rule sums all sub-kernel contributions (additive operadic algebra).

# Example
```julia
using Catlab.Programs: @relation
uwd = @relation (z, z_new) begin
    survive_grow(z, z_new)
    reproduce(z, z_new)
end

sub_kernels = Dict(
    :survive_grow => (z_new, z) -> s(z) * g(z_new, z),
    :reproduce    => (z_new, z) -> f(z) * r(z_new),
)

K = compose_from_uwd(uwd, sub_kernels, domain)
```
"""
function IntegralProjectionModels.compose_from_uwd(
    uwd, sub_kernels::Dict, domain::IntegralProjectionModels.ContinuousDomain)

    z = IntegralProjectionModels.meshpoints(domain)
    h = IntegralProjectionModels.step_size(domain)
    m = length(z)
    K = zeros(m, m)

    # Extract box names from the UWD
    n_boxes = length(boxes(uwd))
    for b in 1:n_boxes
        box_name = Symbol(uwd[:name][b])
        haskey(sub_kernels, box_name) || error(
            "No sub-kernel provided for UWD box :$box_name. " *
            "Available: $(collect(keys(sub_kernels)))")
        kfn = sub_kernels[box_name]
        for j in 1:m
            for i in 1:m
                K[i, j] += h * kfn(z[i], z[j])
            end
        end
    end
    return K
end

"""
    coarsen(A::Matrix, f::FinFunction)

Coarsen a transition matrix via pushforward along a Catlab `FinFunction`.

The `FinFunction` `f: FinSet(n_fine) → FinSet(n_coarse)` maps fine bin
indices to coarse bin indices. The coarse matrix is computed by summing
over fibres and normalising by fibre size:

``A_{\\mathrm{coarse}}[f(i), f(j)] += A[i, j] / |f^{-1}(f(j))|``

# Arguments
- `A`: an `n_fine × n_fine` transition matrix
- `f`: a Catlab `FinFunction` from `FinSet(n_fine)` to `FinSet(n_coarse)`

# Returns
- `A_coarse::Matrix{Float64}`: the coarsened matrix
"""
function IntegralProjectionModels.coarsen(A::Matrix, f::FinFunction)
    n_fine = length(dom(f))
    n_coarse = length(codom(f))
    n_fine == size(A, 1) || throw(DimensionMismatch(
        "Matrix size $(size(A, 1)) does not match FinFunction domain $n_fine"))

    # Compute fibre sizes for normalisation
    fibre_sizes = zeros(Int, n_coarse)
    for i in 1:n_fine
        fibre_sizes[f(i)] += 1
    end

    A_coarse = zeros(n_coarse, n_coarse)
    for j_fine in 1:n_fine
        c_j = f(j_fine)
        for i_fine in 1:n_fine
            c_i = f(i_fine)
            A_coarse[c_i, c_j] += A[i_fine, j_fine] / fibre_sizes[c_j]
        end
    end
    return A_coarse
end

end # module
