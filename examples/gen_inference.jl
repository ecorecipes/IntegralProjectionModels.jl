# IPM inference using Gen.jl
# Uses importance sampling to infer vital rate parameters from population data

using IntegralProjectionModels
using Gen
using Distributions

# --- Setup ---
domain = ContinuousDomain(0.0, 10.0, 20)
n_steps = 10
obs_sigma = 0.5

# --- Generate synthetic observed data ---
# True model
true_params = (s_int = 2.2, s_slope = 0.25, g_int = 0.2, g_slope = 1.02,
    g_sigma = 0.7, f_int = 0.003, f_slope = 0.015)

surv = LinearSurvival(true_params.s_int, true_params.s_slope)
growth = NormalGrowth(true_params.g_int, true_params.g_slope, true_params.g_sigma)
fecund = FecundityRate(true_params.f_int, true_params.f_slope, 0.0, 1.0, 1.0)

P = PKernel(surv, growth, domain)
F = FKernel(fecund, domain)
prob = IPMProblem(P + F, domain, uniform_population(domain), (0, n_steps);
    normalize = true)
sol = solve(prob)

# Observed: log total population size
observed_log_N = [log(sum(u)) for u in sol.u[2:end]]

# --- Create constraint choicemap ---
constraints = Gen.choicemap()
for t in 1:n_steps
    constraints[:log_N => t] = observed_log_N[t]
end

# --- Run importance sampling ---
(traces, log_norm_weights, lml_est) = Gen.importance_resampling(
    ipm_gen_model, (domain, n_steps, obs_sigma),
    constraints, 100)

# Report
println("Log marginal likelihood estimate: $lml_est")
println("Inferred s_int: ", traces[1][:s_int])
println("Inferred g_slope: ", traces[1][:g_slope])
