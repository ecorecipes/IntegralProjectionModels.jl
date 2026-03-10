# Bayesian inference of IPM vital rates using Turing.jl
# Fits survival, growth, and fecundity models from individual-level data,
# then constructs posterior IPMs

using IntegralProjectionModels
using Turing
using Distributions
using Random

Random.seed!(42)

# --- Simulate individual-level data ---
n_obs = 200

# True parameters
true_s_int = 2.2
true_s_slope = 0.25
true_g_int = 0.2
true_g_slope = 1.02
true_g_sigma = 0.7
true_f_int = 0.003
true_f_slope = 0.015

# Generate sizes
sizes = rand(Uniform(0, 50), n_obs)

# Survival data
surv_prob = [1 / (1 + exp(-(true_s_int + true_s_slope * z))) for z in sizes]
surv_status = [rand(Bernoulli(p)) for p in surv_prob]

# Growth data (only for survivors)
survivors = findall(surv_status .== 1)
growth_size_t = sizes[survivors]
growth_size_t1 = [true_g_int + true_g_slope * z + true_g_sigma * randn()
                  for z in growth_size_t]

# Fecundity data
fecund_rates = [exp(true_f_int + true_f_slope * z) for z in sizes]
fecund_seeds = [rand(Poisson(r)) for r in fecund_rates]

# --- Fit model ---
# (Uses the @model from IntegralProjectionModelsTuringExt)
model = ipm_vital_rates(sizes, surv_status,
    growth_size_t, growth_size_t1,
    sizes, fecund_seeds)

chain = sample(model, NUTS(), 500)
println(chain)

# --- Build IPMs from posterior ---
domain = ContinuousDomain(0.0, 50.0, 100)
problems = ipm_from_posterior(chain, domain;
    recruit_mean = 2.0, recruit_sd = 0.3,
    establishment_prob = 1.0, samples = 50)

# Compute lambda for each posterior sample
lambdas = [lambda(solve(prob, EigenAnalysis())) for prob in problems]
println("Posterior lambda: mean=$(mean(lambdas)), sd=$(std(lambdas))")
