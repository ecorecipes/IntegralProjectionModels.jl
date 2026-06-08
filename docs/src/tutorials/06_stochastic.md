# Stochastic IPMs

## Overview

In stochastic IPMs, the environment varies over time, causing the projection kernel to change at each time step. There are two main approaches:

1. **Kernel-resampled**: pre-compute a set of kernels (e.g., one per observed year) and randomly sample from them each time step
2. **Parameter-resampled**: sample environmental parameters each time step and build a new kernel

The **stochastic growth rate** is the geometric mean of per-step growth rates:

$$
\log \lambda_s = \lim_{T \to \infty} \frac{1}{T} \sum_{t=1}^{T} \log \lambda_t
$$

A key result is that $\lambda_s \leq \lambda_{\text{det}}$ — environmental stochasticity always reduces long-term growth relative to the deterministic case (Tuljapurkar's inequality).

## Setup

```@example ipm
using IntegralProjectionModels
using Distributions
using LinearAlgebra
using Statistics
using Random
using Plots

```

## Part 1: Kernel-Resampled Stochastic IPM

We build 5 year-specific kernels for a monocarpic perennial, each with different random intercepts for survival, growth, and fecundity. This follows the `ipmr` test case from `test-simple_di_stoch_kern.R`.

### Parameters

```@example ipm
# Fixed parameters
s_int = 1.03;    s_slope = 2.2
g_int = 8.0;     g_slope = 0.92;  sd_g = 0.9
f_r_int = 0.09;  f_r_slope = 0.05
f_s_int = 0.1;   f_s_slope = 0.005
mu_fd = 9.0;     sd_fd = 2.0

# Year-specific random intercepts
Random.seed!(50127)
g_r = randn(5) .* 0.3     # growth random intercepts
s_r = randn(5) .* 0.7     # survival random intercepts
f_s_r = randn(5) .* 0.2   # fecundity random intercepts

# Domain
L_stoch = 0.2;  U_stoch = 40.0;  n_stoch = 100
domain_stoch = ContinuousDomain(L_stoch, U_stoch, n_stoch)
z_stoch = meshpoints(domain_stoch)
h_stoch = step_size(domain_stoch)
```

### Vital Rate Functions

```@example ipm
# Flowering probability
f_r(z) = 1.0 / (1.0 + exp(-(f_r_int + f_r_slope * z)))

# Survival × (1 - flowering) for year i
function surv_yr(z, s_r_i)
    s = 1.0 / (1.0 + exp(-(s_int + s_slope * z + s_r_i)))
    return s * (1.0 - f_r(z))
end

# Growth for year i (with truncated distributions correction)
function growth_yr(z_prime, z, g_r_i)
    mu = g_int + g_slope * z + g_r_i
    d = Normal(mu, sd_g)
    ev = cdf(d, U_stoch) - cdf(d, L_stoch)
    return ev > 0 ? pdf(d, z_prime) / ev : pdf(d, z_prime)
end

# Offspring size (with truncated distributions correction)
function f_d(z_prime)
    d = Normal(mu_fd, sd_fd)
    ev = cdf(d, U_stoch) - cdf(d, L_stoch)
    return ev > 0 ? pdf(d, z_prime) / ev : pdf(d, z_prime)
end

# Fecundity for year i
function fec_yr(z_prime, z, f_s_r_i)
    return f_r(z) * exp(f_s_int + f_s_slope * z + f_s_r_i) * f_d(z_prime)
end
```

### Build Year-Specific Kernels

```@example ipm
kernels_stoch = []
for i in 1:5
    P_i = PKernel(
        CustomVitalRate(z -> surv_yr(z, s_r[i])),
        CustomVitalRate((z_prime, z) -> growth_yr(z_prime, z, g_r[i])),
        domain_stoch
    )
    F_i = FKernel(
        CustomVitalRate((z_prime, z) -> fec_yr(z_prime, z, f_s_r[i])),
        domain_stoch
    )
    push!(kernels_stoch, P_i + F_i)
end

# Deterministic lambdas for each year
det_lambdas = Float64[]
for (i, k) in enumerate(kernels_stoch)
    K_i = materialize(k)
    e = eigen(K_i)
    λ_i = real(e.values[argmax(real.(e.values))])
    push!(det_lambdas, λ_i)
    println("Year $i: λ = ", round(λ_i, digits=4))
end
```

### Deterministic Lambda (Average Kernel)

```@example ipm
K_matrices = [materialize(k) for k in kernels_stoch]
K_mean = mean(K_matrices)
e_mean = eigen(K_mean)
λ_det = real(e_mean.values[argmax(real.(e_mean.values))])
```

```@example ipm
println("Deterministic λ (mean kernel): ", round(λ_det, digits=6))
```

### Stochastic Simulation

```@example ipm
Random.seed!(42)
n0_stoch = uniform_population(domain_stoch)

prob_stoch = IPMProblem(
    StochasticKernelResampled(),
    kernels_stoch,
    domain_stoch,
    n0_stoch,
    (0, 1000);
    normalize=true
)

sol_stoch = solve(prob_stoch)
```

### Stochastic Growth Rate

```@example ipm
λ_s = stochastic_growth_rate(sol_stoch; burn_in=100)
log_λ_s = log(λ_s)
```

```@example ipm
println("Stochastic growth rate:     λ_s = ", round(λ_s, digits=6))
println("Log stochastic growth rate: log(λ_s) = ", round(log_λ_s, digits=6))
println("Deterministic growth rate:  λ_det = ", round(λ_det, digits=6))
println()
println("Tuljapurkar's inequality: λ_s ≤ λ_det? ", λ_s ≤ λ_det)
```

### Visualize Stochastic Dynamics

```@example ipm
plot(sol_stoch.lambdas[1:200],
    xlabel="Time step", ylabel="λ(t)",
    title="Per-step growth rates (kernel-resampled)",
    label="λ(t)", linewidth=0.5, alpha=0.7)
hline!([λ_s], label="λ_s (stochastic)", color=:red, linewidth=2)
hline!([λ_det], label="λ_det (deterministic)", color=:blue, linewidth=2, linestyle=:dash)
```

```@example ipm
# Running average of log lambda
log_lambdas = log.(sol_stoch.lambdas)
running_mean = [mean(log_lambdas[101:t]) for t in 101:length(log_lambdas)]

plot(101:length(log_lambdas), running_mean,
    xlabel="Time step", ylabel="Running mean log(λ)",
    title="Convergence of stochastic growth rate",
    label="Running mean", linewidth=1)
hline!([log_λ_s], label="log(λ_s)", color=:red, linewidth=2)
```

### Fixed Kernel Sequence

We can also specify a deterministic sequence of kernels for reproducibility.

```@example ipm
Random.seed!(50127)
kern_seq = rand(1:5, 50)

n0_fixed = uniform_population(domain_stoch)
n0_fixed ./= sum(n0_fixed)

prob_fixed = IPMProblem(
    StochasticKernelResampled(),
    kernels_stoch,
    domain_stoch,
    n0_fixed,
    (0, 50);
    normalize=true
)

sol_fixed = solve(prob_fixed; kernel_seq=kern_seq)
```

```@example ipm
println("Lambda sequence (first 10): ", round.(sol_fixed.lambdas[1:10], digits=4))
```

## Part 2: Parameter-Resampled Stochastic IPM

Instead of pre-computing kernels, we sample parameters from environmental distributions at each time step.

```@example ipm
# Environmental sampler: survival and growth intercepts vary
function env_sampler(t)
    return (
        s_offset = 0.7 * randn(),
        g_offset = 0.3 * randn(),
        f_s_offset = 0.2 * randn()
    )
end

# Kernel builder: takes sampled params, returns kernel
function build_stoch_kernel(params)
    function surv_fn(z)
        s = 1.0 / (1.0 + exp(-(s_int + s_slope * z + params.s_offset)))
        return s * (1.0 - f_r(z))
    end

    function growth_fn(z_prime, z)
        mu = g_int + g_slope * z + params.g_offset
        d = Normal(mu, sd_g)
        ev = cdf(d, U_stoch) - cdf(d, L_stoch)
        return ev > 0 ? pdf(d, z_prime) / ev : pdf(d, z_prime)
    end

    function fec_fn(z_prime, z)
        return f_r(z) * exp(f_s_int + f_s_slope * z + params.f_s_offset) * f_d(z_prime)
    end

    P = PKernel(CustomVitalRate(surv_fn), CustomVitalRate(growth_fn), domain_stoch)
    F = FKernel(CustomVitalRate(fec_fn), domain_stoch)
    return P + F
end
```

```@example ipm
Random.seed!(12345)
n0_param = uniform_population(domain_stoch)

prob_param = IPMProblem(
    StochasticParameterResampled(),
    build_stoch_kernel,
    domain_stoch,
    n0_param,
    (0, 1000);
    env_state=env_sampler,
    normalize=true
)

sol_param = solve(prob_param)
```

```@example ipm
λ_s_param = stochastic_growth_rate(sol_param; burn_in=100)
```

```@example ipm
println("Stochastic growth rate (parameter-resampled): ", round(λ_s_param, digits=6))
```

```@example ipm
plot(sol_param.lambdas[1:200],
    xlabel="Time step", ylabel="λ(t)",
    title="Per-step growth rates (parameter-resampled)",
    label="λ(t)", linewidth=0.5, alpha=0.7)
hline!([λ_s_param], label="λ_s", color=:red, linewidth=2)
```

## Comparing the Two Approaches

| Feature | Kernel-resampled | Parameter-resampled |
|---------|-----------------|---------------------|
| Pre-computation | Kernels built once | Kernel built each step |
| Environmental variation | Discrete set of environments | Continuous distribution |
| Reproducibility | Easy with fixed sequence | Requires seeded RNG |
| Speed | Faster (reuse materialized kernels) | Slower (materialize each step) |
| Flexibility | Limited to observed environments | Can model any distribution |

## Effect of Environmental Variance

Let us examine how increasing environmental variance affects the stochastic growth rate.

```@example ipm
Random.seed!(999)
σ_vals = range(0.0, 1.5, length=8)
λ_s_by_σ = Float64[]

for σ in σ_vals
    function env_σ(t)
        return (s_offset = σ * randn(), g_offset = 0.3 * randn(), f_s_offset = 0.2 * randn())
    end

    prob_σ = IPMProblem(
        StochasticParameterResampled(),
        build_stoch_kernel,
        domain_stoch,
        uniform_population(domain_stoch),
        (0, 500);
        env_state=env_σ,
        normalize=true
    )
    sol_σ = solve(prob_σ)
    push!(λ_s_by_σ, stochastic_growth_rate(sol_σ; burn_in=50))
end

plot(σ_vals, λ_s_by_σ,
    xlabel="Survival SD (σ)",
    ylabel="Stochastic λ_s",
    title="Stochastic growth rate vs environmental variance",
    label="λ_s", linewidth=2, marker=:circle)
hline!([λ_s_by_σ[1]], label="Deterministic λ (σ=0)", linestyle=:dash, color=:red)
```

As environmental variance increases, the stochastic growth rate decreases, consistent with Tuljapurkar's inequality.

## Summary

In this vignette we:

1. Built a **kernel-resampled** stochastic IPM with 5 year-specific kernels
2. Computed the stochastic growth rate $\lambda_s = \exp(\overline{\log \lambda_t})$
3. Verified Tuljapurkar's inequality: $\lambda_s \leq \lambda_{\text{det}}$
4. Built a **parameter-resampled** stochastic IPM with continuous environmental distributions
5. Demonstrated that increasing environmental variance reduces $\lambda_s$
