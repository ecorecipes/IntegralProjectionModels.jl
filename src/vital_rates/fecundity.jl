"""
    FecundityRate{T}

Combined fecundity kernel computing the probability of producing an offspring
of size `z'` from a parent of size `z`:

`f(z', z) = establishment_prob * exp(intercept + slope * z) * pdf(Normal(recruit_mean, recruit_sd), z')`
"""
struct FecundityRate{T} <: AbstractFecundityRate
    intercept::T
    slope::T
    recruit_mean::T
    recruit_sd::T
    establishment_prob::T
end

function FecundityRate(intercept, slope, recruit_mean, recruit_sd, establishment_prob)
    T = promote_type(typeof(intercept), typeof(slope), typeof(recruit_mean),
        typeof(recruit_sd), typeof(establishment_prob))
    FecundityRate{T}(T(intercept), T(slope), T(recruit_mean), T(recruit_sd),
        T(establishment_prob))
end

function (f::FecundityRate)(z_prime, z)
    seed_prod = exp(f.intercept + f.slope * z)
    recruit_dist = pdf(Normal(f.recruit_mean, f.recruit_sd), z_prime)
    return f.establishment_prob * seed_prod * recruit_dist
end

"""
    LogisticFecundityRate{T}

Fecundity where reproduction probability is logistic (e.g., flowering probability):

`f(z', z) = logistic(repr_int + repr_slope * z) * exp(seed_int + seed_slope * z) * pdf(Normal(recruit_mean, recruit_sd), z')`
"""
struct LogisticFecundityRate{T} <: AbstractFecundityRate
    repr_int::T
    repr_slope::T
    seed_int::T
    seed_slope::T
    recruit_mean::T
    recruit_sd::T
end

function LogisticFecundityRate(repr_int, repr_slope, seed_int, seed_slope,
        recruit_mean, recruit_sd)
    T = promote_type(typeof(repr_int), typeof(repr_slope), typeof(seed_int),
        typeof(seed_slope), typeof(recruit_mean), typeof(recruit_sd))
    LogisticFecundityRate{T}(T(repr_int), T(repr_slope), T(seed_int), T(seed_slope),
        T(recruit_mean), T(recruit_sd))
end

function (f::LogisticFecundityRate)(z_prime, z)
    f_r = logistic(f.repr_int + f.repr_slope * z)
    f_s = exp(f.seed_int + f.seed_slope * z)
    f_d = pdf(Normal(f.recruit_mean, f.recruit_sd), z_prime)
    return f_r * f_s * f_d
end

"""
    RecruitmentDistribution{T}

Distribution of recruit sizes, independent of parent size.
"""
struct RecruitmentDistribution{T} <: AbstractRecruitmentRate
    mean::T
    sd::T
end

function (r::RecruitmentDistribution)(z_prime)
    return pdf(Normal(r.mean, r.sd), z_prime)
end
(r::RecruitmentDistribution)(z_prime, _) = r(z_prime)
