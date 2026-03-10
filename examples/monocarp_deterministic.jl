# Monocarp (Oenothera glazioviana) Deterministic IPM
# Reproduces the classic textbook example from Ellner, Childs & Rees

using IntegralProjectionModels
using Distributions

# --- Parameters ---
s_int = 1.03;  s_slope = 0.19      # survival (logit scale)
g_int = 8.0;   g_slope = 0.92;  sd_g = 0.9  # growth
f_r_int = 0.09; f_r_slope = 0.05   # flowering probability (logit)
f_s_int = 0.01; f_s_slope = 0.0005 # seed production (log)
mu_fd = 3.0;    sd_fd = 0.7        # recruit size distribution

# --- Domain ---
domain = ContinuousDomain(0.3, 200.0, 500)

# --- Vital rates ---
# Survival × (1 - flowering probability)
total_survival(z) = (1 / (1 + exp(-(s_int + s_slope * z)))) *
                    (1 - 1 / (1 + exp(-(f_r_int + f_r_slope * z))))

# Growth: Normal(g_int + g_slope * z, sd_g)
growth(z_prime, z) = pdf(Normal(g_int + g_slope * z, sd_g), z_prime)

# Fecundity: f_r(z) * f_s(z) * f_d(z')
function fecundity(z_prime, z)
    f_r = 1 / (1 + exp(-(f_r_int + f_r_slope * z)))
    f_s = exp(f_s_int + f_s_slope * z)
    f_d = pdf(Normal(mu_fd, sd_fd), z_prime)
    return f_r * f_s * f_d
end

# --- Build kernels ---
P = PKernel(CustomVitalRate(total_survival), CustomVitalRate(growth), domain)
F = FKernel(CustomVitalRate(fecundity), domain)
kernel = P + F

# --- Solve ---
n0 = uniform_population(domain)
prob = IPMProblem(kernel, domain, n0, (0, 100))

# Eigenanalysis
sol = solve(prob, EigenAnalysis())
println("Asymptotic growth rate (λ): ", lambda(sol))

# Iteration
sol_iter = solve(prob, DirectIteration())
println("Lambda from iteration: ", sol_iter.lambdas[end])

# Sensitivity and elasticity
S = sensitivity(sol)
E = elasticity(sol)
println("Sum of elasticities: ", sum(E))
