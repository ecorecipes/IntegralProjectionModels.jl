"""
Utility functions for IPM construction and analysis.
"""

"""
    rate_to_proportion(rate)

Convert a rate (potentially > 1) to a proportion (0, 1) via `1 - exp(-rate)`.
"""
rate_to_proportion(rate) = 1 - exp(-rate)

"""
    proportion_to_rate(prop)

Convert a proportion (0, 1) to a rate via `-log(1 - prop)`.
"""
proportion_to_rate(prop) = -log(1 - prop)

"""
    uniform_population(domain::ContinuousDomain)

Create a uniform initial population vector with total mass 1.
"""
function uniform_population(domain::ContinuousDomain)
    n = n_states(domain)
    return fill(1.0 / n, n)
end

"""
    point_population(domain::ContinuousDomain, z0::Real)

Create an initial population concentrated at size `z0`.
"""
function point_population(domain::ContinuousDomain, z0::Real)
    z = meshpoints(domain)
    idx = argmin(abs.(z .- z0))
    n = zeros(n_states(domain))
    n[idx] = 1.0
    return n
end

"""
    normal_population(domain::ContinuousDomain, mean, sd)

Create an initial population following a normal distribution.
"""
function normal_population(domain::ContinuousDomain, mean, sd)
    z = meshpoints(domain)
    h = step_size(domain)
    n = [pdf(Normal(mean, sd), zi) * h for zi in z]
    s = sum(n)
    return s > 0 ? n ./ s : n
end
