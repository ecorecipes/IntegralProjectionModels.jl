module IntegralProjectionModelsTuringExt

using IntegralProjectionModels
using Turing
using Distributions
using LinearAlgebra

"""
    ipm_vital_rates(surv_size, surv_status, growth_size_t, growth_size_t1,
                    fecund_size, fecund_seeds)

Turing model for fitting IPM vital rates from individual-level data.

# Arguments
- `surv_size`: size of individuals observed for survival
- `surv_status`: 1 = survived, 0 = died
- `growth_size_t`: size at time t (for growth data)
- `growth_size_t1`: size at time t+1 (for growth data)
- `fecund_size`: size of reproducing individuals
- `fecund_seeds`: seed/offspring count
"""
Turing.@model function ipm_vital_rates(surv_size, surv_status,
        growth_size_t, growth_size_t1,
        fecund_size, fecund_seeds)
    # Survival priors
    s_int ~ Normal(0, 5)
    s_slope ~ Normal(0, 5)

    # Growth priors
    g_int ~ Normal(0, 5)
    g_slope ~ Normal(0, 5)
    g_sigma ~ truncated(Normal(0, 2); lower = 0.01)

    # Fecundity priors
    f_int ~ Normal(0, 5)
    f_slope ~ Normal(0, 5)

    # Survival likelihood
    for i in eachindex(surv_size)
        p = StatsFuns.logistic(s_int + s_slope * surv_size[i])
        surv_status[i] ~ Bernoulli(p)
    end

    # Growth likelihood
    for i in eachindex(growth_size_t)
        mu = g_int + g_slope * growth_size_t[i]
        growth_size_t1[i] ~ Normal(mu, g_sigma)
    end

    # Fecundity likelihood (Poisson)
    for i in eachindex(fecund_size)
        rate = exp(f_int + f_slope * fecund_size[i])
        fecund_seeds[i] ~ Poisson(rate)
    end
end

"""
    ipm_from_posterior(chain, domain; samples=nothing)

Build IPMProblems from posterior samples.

Returns a vector of IPMProblem instances, one per posterior sample,
suitable for posterior predictive checks or stochastic parameter-resampled IPMs.
"""
function IntegralProjectionModels.ipm_from_posterior(chain, domain::ContinuousDomain;
        samples = nothing, recruit_mean = 0.0, recruit_sd = 1.0,
        establishment_prob = 1.0, n0 = nothing, tspan = (0, 100))
    if n0 === nothing
        n0 = IntegralProjectionModels.uniform_population(domain)
    end

    n_samples = length(chain[:s_int])
    if samples !== nothing
        indices = rand(1:n_samples, samples)
    else
        indices = 1:n_samples
    end

    problems = Vector{IPMProblem}(undef, length(indices))
    for (k, i) in enumerate(indices)
        surv = LinearSurvival(chain[:s_int][i], chain[:s_slope][i])
        growth = NormalGrowth(chain[:g_int][i], chain[:g_slope][i], chain[:g_sigma][i])
        fecund = FecundityRate(chain[:f_int][i], chain[:f_slope][i],
            recruit_mean, recruit_sd, establishment_prob)
        P = PKernel(surv, growth, domain)
        F = FKernel(fecund, domain)
        problems[k] = IPMProblem(P + F, domain, n0, tspan)
    end

    return problems
end

export ipm_vital_rates, ipm_from_posterior

end # module
