"""
PADRINO eviction correction mapping.

Maps PADRINO eviction types to IntegralProjectionModels.jl eviction corrections.
"""

"""
    map_eviction_type(evict_type::String) -> EvictionCorrection

Map a PADRINO eviction type string to the corresponding
IntegralProjectionModels.jl EvictionCorrection enum value.
"""
function map_eviction_type(evict_type::String)
    if evict_type == "truncated_distributions"
        return IPM.TruncatedDistributions
    elseif evict_type == "discrete_extrema"
        return IPM.DiscreteExtrema
    elseif evict_type == "stretched_domain"
        return IPM.NoCorrection
    elseif evict_type == "rescale_kernel"
        # rescale_kernel is similar to discrete_extrema
        return IPM.DiscreteExtrema
    else
        return IPM.NoCorrection
    end
end

"""
    apply_eviction!(K::Matrix{Float64}, evict_type::String, domain::ContinuousDomain)

Apply eviction correction to a kernel matrix.
"""
function apply_eviction!(K::Matrix{Float64}, evict_type::String, domain::IPM.ContinuousDomain)
    correction = map_eviction_type(evict_type)

    if correction == IPM.DiscreteExtrema
        IPM.apply_discrete_extrema!(K)
    elseif correction == IPM.TruncatedDistributions
        # For truncated distributions, the correction is applied during
        # the distribution evaluation (the distributions are wrapped in
        # Truncated(dist, lower, upper) during translation)
        # Nothing to do here post-hoc
    end

    return K
end
