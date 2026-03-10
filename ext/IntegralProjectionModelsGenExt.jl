module IntegralProjectionModelsGenExt

using IntegralProjectionModels
using Gen
using Distributions
using LinearAlgebra

"""
    ipm_gen_model(domain, n_steps, obs_sigma)

Gen generative function for IPM inference.
Samples vital rate parameters, builds kernel, iterates population,
and generates observations.
"""
@gen function ipm_gen_model(domain::ContinuousDomain, n_steps::Int, obs_sigma::Float64)
    # Sample vital rate parameters
    s_int = {:s_int} ~ normal(0, 5)
    s_slope = {:s_slope} ~ normal(0, 5)
    g_int = {:g_int} ~ normal(0, 5)
    g_slope = {:g_slope} ~ normal(0, 5)
    g_sigma = {:g_sigma} ~ Gen.gamma(2, 1)
    f_int = {:f_int} ~ normal(0, 5)
    f_slope = {:f_slope} ~ normal(0, 5)

    # Build kernel
    surv = LinearSurvival(s_int, s_slope)
    growth = NormalGrowth(g_int, g_slope, g_sigma)
    fecund = FecundityRate(f_int, f_slope, 0.0, 1.0, 1.0)
    P = PKernel(surv, growth, domain)
    F = FKernel(fecund, domain)
    K = materialize(P + F)

    # Iterate and observe
    m = IntegralProjectionModels.n_states(domain)
    n = IntegralProjectionModels.uniform_population(domain)

    for t in 1:n_steps
        n = K * n
        s = sum(n)
        if s > 0
            n ./= s
        end
        # Observe total population (log scale)
        {:log_N, t} ~ normal(log(s), obs_sigma)
    end

    return n
end

export ipm_gen_model

end # module
