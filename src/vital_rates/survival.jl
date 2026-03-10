"""
    LinearSurvival{T}

Survival probability as a logistic function of size:
`s(z) = logistic(intercept + slope * z)`.
"""
struct LinearSurvival{T} <: AbstractSurvivalRate
    intercept::T
    slope::T
end

function LinearSurvival(intercept, slope)
    T = promote_type(typeof(intercept), typeof(slope))
    LinearSurvival{T}(T(intercept), T(slope))
end

(s::LinearSurvival)(z) = logistic(s.intercept + s.slope * z)

"""
    QuadraticSurvival{T}

Survival probability with a quadratic term:
`s(z) = logistic(intercept + slope * z + quadratic * z^2)`.
"""
struct QuadraticSurvival{T} <: AbstractSurvivalRate
    intercept::T
    slope::T
    quadratic::T
end

function QuadraticSurvival(intercept, slope, quadratic)
    T = promote_type(typeof(intercept), typeof(slope), typeof(quadratic))
    QuadraticSurvival{T}(T(intercept), T(slope), T(quadratic))
end

(s::QuadraticSurvival)(z) = logistic(s.intercept + s.slope * z + s.quadratic * z^2)

"""
    ConstantSurvival{T}

Constant survival probability independent of size.
"""
struct ConstantSurvival{T} <: AbstractSurvivalRate
    prob::T
end

(s::ConstantSurvival)(_) = s.prob
