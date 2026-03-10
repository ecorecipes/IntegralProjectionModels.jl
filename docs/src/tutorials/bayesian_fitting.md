# Bayesian IPM Fitting

IntegralProjectionModels.jl integrates with Turing.jl for Bayesian inference of vital rate parameters.

## Fitting Vital Rates

Given individual-level data on survival, growth, and fecundity:

```julia
using IntegralProjectionModels, Turing

# Assume you have:
# surv_size, surv_status - survival data
# growth_size_t, growth_size_t1 - growth data
# fecund_size, fecund_seeds - fecundity data

model = ipm_vital_rates(surv_size, surv_status,
    growth_size_t, growth_size_t1,
    fecund_size, fecund_seeds)

chain = sample(model, NUTS(), 1000)
```

## Posterior IPMs

Convert posterior samples to IPMs:

```julia
domain = ContinuousDomain(0.0, 50.0, 100)
problems = ipm_from_posterior(chain, domain;
    recruit_mean = 2.0, recruit_sd = 0.3,
    samples = 100)

# Posterior distribution of lambda
lambdas = [lambda(solve(p, EigenAnalysis())) for p in problems]
```

## Gen.jl Alternative

For programmatic Bayesian inference with Gen.jl:

```julia
using IntegralProjectionModels, Gen

# Use ipm_gen_model for importance sampling or MCMC
(traces, weights, lml) = importance_resampling(
    ipm_gen_model, (domain, n_steps, obs_sigma),
    constraints, 1000)
```
