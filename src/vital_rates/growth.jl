"""
    NormalGrowth{T}

Growth kernel: `g(z', z) = pdf(Normal(intercept + slope * z, sigma), z')`.
Transition from size `z` to size `z'` follows a normal distribution.
"""
struct NormalGrowth{T} <: AbstractGrowthRate
    intercept::T
    slope::T
    sigma::T
end

function NormalGrowth(intercept, slope, sigma)
    T = promote_type(typeof(intercept), typeof(slope), typeof(sigma))
    NormalGrowth{T}(T(intercept), T(slope), T(sigma))
end

function (g::NormalGrowth)(z_prime, z)
    mu = g.intercept + g.slope * z
    return pdf(Normal(mu, g.sigma), z_prime)
end

"""
    mean_size(g::NormalGrowth, z)

Expected size at next time step given current size `z`.
"""
mean_size(g::NormalGrowth, z) = g.intercept + g.slope * z

"""
    LogNormalGrowth{T}

Growth kernel using a log-normal distribution for size transitions.
"""
struct LogNormalGrowth{T} <: AbstractGrowthRate
    intercept::T
    slope::T
    sigma::T
end

function LogNormalGrowth(intercept, slope, sigma)
    T = promote_type(typeof(intercept), typeof(slope), typeof(sigma))
    LogNormalGrowth{T}(T(intercept), T(slope), T(sigma))
end

function (g::LogNormalGrowth)(z_prime, z)
    mu = g.intercept + g.slope * z
    return pdf(LogNormal(mu, g.sigma), z_prime)
end
