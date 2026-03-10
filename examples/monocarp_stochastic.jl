# Stochastic IPM with parameter-resampled environmental variation

using IntegralProjectionModels
using Distributions
using Statistics

# Domain
domain = ContinuousDomain(0.3, 200.0, 200)

# Fixed parameters
g_int = 8.0;  g_slope = 0.92;  sd_g = 0.9
f_r_int = 0.09; f_r_slope = 0.05
f_s_int = 0.01; f_s_slope = 0.0005
mu_fd = 3.0;  sd_fd = 0.7

# Kernel builder: takes sampled params, returns kernel
function build_kernel(params)
    s_int = params.s_int
    s_slope = params.s_slope

    surv_fn(z) = (1 / (1 + exp(-(s_int + s_slope * z)))) *
                 (1 - 1 / (1 + exp(-(f_r_int + f_r_slope * z))))

    growth_fn(z_prime, z) = pdf(Normal(g_int + g_slope * z, sd_g), z_prime)

    function fecund_fn(z_prime, z)
        f_r = 1 / (1 + exp(-(f_r_int + f_r_slope * z)))
        f_s = exp(f_s_int + f_s_slope * z)
        return f_r * f_s * pdf(Normal(mu_fd, sd_fd), z_prime)
    end

    P = PKernel(CustomVitalRate(surv_fn), CustomVitalRate(growth_fn), domain)
    F = FKernel(CustomVitalRate(fecund_fn), domain)
    return P + F
end

# Environment sampler: survival parameters vary stochastically
env_sampler(t) = (s_int = 1.03 + 0.3 * randn(), s_slope = 0.19 + 0.02 * randn())

# Solve
n0 = uniform_population(domain)
prob = IPMProblem(StochasticParameterResampled(),
    build_kernel, domain, n0, (0, 500);
    env_state = env_sampler, normalize = true)

sol = solve(prob)

# Stochastic growth rate
λ_s = stochastic_growth_rate(sol; burn_in = 100)
println("Stochastic growth rate: ", λ_s)
println("Log stochastic growth rate: ", log(λ_s))
