# Monocarp Plant IPM

This tutorial reproduces the classic monocarp (*Oenothera glazioviana*) example from the textbook "Data-driven Modeling of Structured Populations" (Ellner, Childs & Rees).

## Parameters

```julia
using IntegralProjectionModels, Distributions

# Vital rate parameters (from field data)
s_int = 1.03;  s_slope = 0.19       # survival
g_int = 8.0;   g_slope = 0.92;  sd_g = 0.9   # growth
f_r_int = 0.09; f_r_slope = 0.05    # flowering probability
f_s_int = 0.01; f_s_slope = 0.0005  # seed production
mu_fd = 3.0;    sd_fd = 0.7         # recruit size
```

## Model Construction

The monocarp model uses a modified survival that accounts for flowering probability (monocarpic plants die after flowering):

```julia
domain = ContinuousDomain(0.3, 200.0, 500)

# Custom vital rates
total_survival(z) = (1/(1+exp(-(s_int + s_slope*z)))) *
                    (1 - 1/(1+exp(-(f_r_int + f_r_slope*z))))

growth(z_prime, z) = pdf(Normal(g_int + g_slope*z, sd_g), z_prime)

function fecundity(z_prime, z)
    f_r = 1/(1+exp(-(f_r_int + f_r_slope*z)))
    f_s = exp(f_s_int + f_s_slope*z)
    return f_r * f_s * pdf(Normal(mu_fd, sd_fd), z_prime)
end

P = PKernel(CustomVitalRate(total_survival), CustomVitalRate(growth), domain)
F = FKernel(CustomVitalRate(fecundity), domain)

prob = IPMProblem(P + F, domain, uniform_population(domain), (0, 100))
sol = solve(prob, EigenAnalysis())
println("λ = ", lambda(sol))
```
