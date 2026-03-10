# Density-Dependent IPM
# Survival decreases as population density increases

using IntegralProjectionModels
using Distributions

# Domain
domain = ContinuousDomain(0.0, 50.0, 100)

# Fixed parameters
g_int = 0.2;  g_slope = 1.02;  sd_g = 0.7
f_int = 0.003; f_slope = 0.015
mu_fd = 2.0;  sd_fd = 0.3

# Density-dependent kernel builder
function dd_kernel(n_t, t, p)
    total_N = sum(n_t)
    # Survival intercept decreases with population density
    s_int = 2.2 - 0.001 * total_N
    s_slope = 0.25

    survival(z) = 1 / (1 + exp(-(s_int + s_slope * z)))
    growth(z_prime, z) = pdf(Normal(g_int + g_slope * z, sd_g), z_prime)
    function fecundity(z_prime, z)
        exp(f_int + f_slope * z) * pdf(Normal(mu_fd, sd_fd), z_prime)
    end

    P = PKernel(CustomVitalRate(survival), CustomVitalRate(growth), domain)
    F = FKernel(CustomVitalRate(fecundity), domain)
    return P + F
end

# Solve
n0 = ones(100) .* 10.0
prob = IPMProblem(DensityDependent(), dd_kernel, domain, n0, (0, 200))
sol = solve(prob)

println("Final lambda: ", sol.lambdas[end])
println("Initial population: ", sum(sol.u[1]))
println("Final population: ", sum(sol.u[end]))
