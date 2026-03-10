"""
    ContinuousDomain{T<:Real}

Represents a continuous state variable domain discretized via the midpoint rule.

The domain `[lower, upper]` is divided into `n_meshpoints` bins.
Meshpoints are the midpoints of each bin.
"""
struct ContinuousDomain{T <: Real}
    lower::T
    upper::T
    n_meshpoints::Int

    function ContinuousDomain(lower::T, upper::T, n_meshpoints::Int) where {T <: Real}
        lower < upper || throw(ArgumentError("lower must be less than upper"))
        n_meshpoints > 0 || throw(ArgumentError("n_meshpoints must be positive"))
        new{T}(lower, upper, n_meshpoints)
    end
end

function ContinuousDomain(lower::Real, upper::Real, n_meshpoints::Int)
    T = promote_type(typeof(lower), typeof(upper))
    ContinuousDomain(T(lower), T(upper), n_meshpoints)
end

"""
    step_size(d::ContinuousDomain)

Width of each bin: `(upper - lower) / n_meshpoints`.
"""
step_size(d::ContinuousDomain) = (d.upper - d.lower) / d.n_meshpoints

"""
    meshpoints(d::ContinuousDomain)

Return midpoints of the `n_meshpoints` bins spanning `[lower, upper]`.

Matches ipmr convention: `bounds = linspace(L, U, n+1)`, `midpoints = (bounds[1:n] + bounds[2:n+1]) / 2`.
"""
function meshpoints(d::ContinuousDomain{T}) where {T}
    n = d.n_meshpoints
    bounds = range(d.lower, d.upper; length = n + 1)
    return [(bounds[i] + bounds[i + 1]) / 2 for i in 1:n]
end

"""
    bounds(d::ContinuousDomain)

Return the `n_meshpoints + 1` bin edges spanning `[lower, upper]`.
"""
function bounds(d::ContinuousDomain)
    return collect(range(d.lower, d.upper; length = d.n_meshpoints + 1))
end

"""
    DiscreteDomain

Represents a discrete state variable with named levels.
"""
struct DiscreteDomain
    labels::Vector{Symbol}
end

n_states(d::DiscreteDomain) = length(d.labels)
n_states(d::ContinuousDomain) = d.n_meshpoints
