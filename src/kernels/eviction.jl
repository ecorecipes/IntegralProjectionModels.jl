"""
Eviction correction methods for IPM kernels.

Eviction occurs when individuals are "lost" outside domain bounds because
growth transitions produce sizes beyond the modeled range. Correction
ensures column sums of the discretized kernel are preserved.
"""

"""
    EvictionCorrection

Enum specifying the type of eviction correction to apply.
- `NoCorrection`: no adjustment
- `TruncatedDistributions`: divide growth density by CDF within bounds
- `DiscreteExtrema`: redistribute lost mass to domain boundaries
"""
@enum EvictionCorrection begin
    NoCorrection
    TruncatedDistributions
    DiscreteExtrema
end

"""
    truncated_growth(g::AbstractGrowthRate, z_prime, z, domain::ContinuousDomain)

Evaluate growth kernel with truncated distribution correction.
Divides the density by the probability mass within the domain bounds,
ensuring the growth distribution integrates to 1 over the domain.
"""
function truncated_growth(g::NormalGrowth, z_prime, z, domain::ContinuousDomain)
    mu = g.intercept + g.slope * z
    d = Normal(mu, g.sigma)
    raw = pdf(d, z_prime)
    correction = cdf(d, domain.upper) - cdf(d, domain.lower)
    return correction > 0 ? raw / correction : raw
end

function truncated_growth(g::LogNormalGrowth, z_prime, z, domain::ContinuousDomain)
    mu = g.intercept + g.slope * z
    d = LogNormal(mu, g.sigma)
    raw = pdf(d, z_prime)
    correction = cdf(d, domain.upper) - cdf(d, domain.lower)
    return correction > 0 ? raw / correction : raw
end

function truncated_growth(g::AbstractGrowthRate, z_prime, z, domain::ContinuousDomain)
    g(z_prime, z)
end

"""
    truncated_recruit(r, z_prime, domain::ContinuousDomain)

Apply truncation correction to a recruitment/offspring size distribution.
"""
function truncated_recruit(recruit_mean, recruit_sd, z_prime, domain::ContinuousDomain)
    d = Normal(recruit_mean, recruit_sd)
    raw = pdf(d, z_prime)
    correction = cdf(d, domain.upper) - cdf(d, domain.lower)
    return correction > 0 ? raw / correction : raw
end

"""
    apply_discrete_extrema!(K::AbstractMatrix)

Discrete extrema eviction correction: redistribute lost probability mass
to domain boundaries so that column integrals are preserved.
"""
function apply_discrete_extrema!(K::AbstractMatrix)
    m = size(K, 1)
    mid = m ÷ 2
    for j in axes(K, 2)
        col_sum = sum(view(K, :, j))
        if col_sum > 0
            deficit = one(col_sum) - col_sum
            if j <= mid
                K[1, j] += deficit
            else
                K[m, j] += deficit
            end
        end
    end
    return K
end
