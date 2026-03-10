# Age × Size structured IPM
# Population structured by both continuous size and discrete age

using IntegralProjectionModels
using Distributions
using LinearAlgebra

# Domain and age structure
domain = ContinuousDomain(0.0, 20.0, 50)
age_struct = AgeStructure(5)

# Age-specific P kernels (survival-growth)
function p_func(age)
    # Survival increases with age, then plateaus
    s = LinearSurvival(0.5 + 0.3 * min(age, 3), 0.1)
    g = NormalGrowth(0.5, 0.9, 0.5)
    PKernel(s, g, domain)
end

# Age-specific F kernels (fecundity)
function f_func(age)
    if age <= 1
        # Juveniles don't reproduce
        return CustomKernel((z_prime, z) -> 0.0, domain)
    end
    # Fecundity increases with age
    f = FecundityRate(0.01 * age, 0.005, 2.0, 0.5, 1.0)
    FKernel(f, domain)
end

# Build expanded kernel matrix
K = expand_age_kernels(p_func, f_func, age_struct, domain)

# Initial population: small, young
m = n_states(domain)
total = n_ages(age_struct) * m
n0 = zeros(total)
n0[1:m] .= 1.0 / m  # all in age class 1

# Iterate
u = Vector{Vector{Float64}}(undef, 101)
u[1] = copy(n0)
for t in 1:100
    u[t + 1] = K * u[t]
end

# Compute lambda
lambdas = [sum(u[t + 1]) / sum(u[t]) for t in 1:100]
println("Asymptotic λ: ", lambdas[end])
println("Eigenvalue λ: ", real(eigen(K).values[end]))
