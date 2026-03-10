# Stochastic IPMs

IntegralProjectionModels.jl supports two types of environmental stochasticity.

## Kernel-Resampled Stochasticity

Pre-define a set of kernels (e.g., from different years) and randomly sample one each time step:

```julia
using IntegralProjectionModels, Distributions

domain = ContinuousDomain(0.0, 50.0, 100)

# Build kernels for "good" and "bad" years
g = NormalGrowth(0.2, 1.02, 0.7)
f = FecundityRate(0.003, 0.015, 2.0, 0.3, 1.0)

kern_good = PKernel(LinearSurvival(2.5, 0.3), g, domain) + FKernel(f, domain)
kern_bad  = PKernel(LinearSurvival(1.5, 0.2), g, domain) + FKernel(f, domain)

prob = IPMProblem(StochasticKernelResampled(),
    [kern_good, kern_bad], domain,
    uniform_population(domain), (0, 500))

sol = solve(prob)
λ_s = stochastic_growth_rate(sol; burn_in = 100)
```

## Parameter-Resampled Stochasticity

Sample vital rate parameters from distributions each time step:

```julia
function build_kernel(params)
    s = LinearSurvival(params.s_int, 0.25)
    g = NormalGrowth(0.2, 1.02, 0.7)
    f = FecundityRate(0.003, 0.015, 2.0, 0.3, 1.0)
    PKernel(s, g, domain) + FKernel(f, domain)
end

env_sampler(t) = (s_int = 2.0 + 0.5*randn(),)

prob = IPMProblem(StochasticParameterResampled(),
    build_kernel, domain,
    uniform_population(domain), (0, 500);
    env_state = env_sampler, normalize = true)

sol = solve(prob)
λ_s = stochastic_growth_rate(sol; burn_in = 100)
```
