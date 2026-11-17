module IntegralProjectionModelsMTKExt

using IntegralProjectionModels
using ModelingToolkit
using Symbolics
using Distributions: Normal, pdf

"""
    IPMModel

Symbolic representation of an IPM, defined via ModelingToolkit parameters
and expressions. Convert to a concrete IPMProblem by supplying parameter values.

# Fields
- `survival_expr`: symbolic survival expression
- `growth_mean_expr`: symbolic growth mean expression
- `growth_sd`: symbolic growth sd parameter
- `fecundity_expr`: symbolic fecundity expression (optional)
- `recruit_mean`: symbolic recruit mean parameter
- `recruit_sd`: symbolic recruit sd parameter
- `parameters`: vector of symbolic parameters
"""
struct IPMModel{S, GM, GS, F, RM, RS, P}
    survival_expr::S
    growth_mean_expr::GM
    growth_sd::GS
    fecundity_expr::F
    recruit_mean::RM
    recruit_sd::RS
    parameters::P
end

"""
    symbolic_ipm(; name=:ipm)

Create a symbolic IPM model with standard parameterization.
Returns an IPMModel with symbolic parameters that can be bound to concrete values.
"""
function symbolic_ipm(; name = :ipm)
    @parameters s_int s_slope g_int g_slope g_sigma f_int f_slope r_mean r_sd p_est
    @variables z

    survival = 1 / (1 + exp(-(s_int + s_slope * z)))
    growth_mean = g_int + g_slope * z
    fecundity = p_est * exp(f_int + f_slope * z)

    return IPMModel(survival, growth_mean, g_sigma, fecundity, r_mean, r_sd,
        [s_int, s_slope, g_int, g_slope, g_sigma, f_int, f_slope, r_mean, r_sd, p_est])
end

"""
    IPMProblem(model::IPMModel, domain, n0, tspan; p)

Convert a symbolic IPMModel to a concrete IPMProblem by substituting parameter values.
`p` should be a Dict mapping symbolic parameters to numeric values.
"""
function IntegralProjectionModels.IPMProblem(model::IPMModel, domain::ContinuousDomain,
        n0, tspan; p::Dict, kwargs...)
    # Compile symbolic expressions to Julia functions. `get_variables` returns the parameters
    # as well as the size variable, so keep the one symbol that is not a parameter.
    z_var = only(filter(v -> !ModelingToolkit.isparameter(v),
        Symbolics.get_variables(model.growth_mean_expr)))

    surv_fn = Symbolics.build_function(
        Symbolics.substitute(model.survival_expr, p),
        z_var; expression = Val{false})

    gmean_fn = Symbolics.build_function(
        Symbolics.substitute(model.growth_mean_expr, p),
        z_var; expression = Val{false})

    # substitution leaves a symbolic constant; unwrap it before converting
    number(x) = Float64(Symbolics.value(Symbolics.substitute(x, p)))
    g_sigma_val = number(model.growth_sd)
    r_mean_val = number(model.recruit_mean)
    r_sd_val = number(model.recruit_sd)

    fecund_fn = Symbolics.build_function(
        Symbolics.substitute(model.fecundity_expr, p),
        z_var; expression = Val{false})

    # Build vital rates
    survival = CustomVitalRate(z -> surv_fn(z))
    growth = CustomVitalRate((z_prime, z) -> begin
        mu = gmean_fn(z)
        pdf(Normal(mu, g_sigma_val), z_prime)
    end)
    fecundity = CustomVitalRate((z_prime, z) -> begin
        fecund_fn(z) * pdf(Normal(r_mean_val, r_sd_val), z_prime)
    end)

    P = PKernel(survival, growth, domain)
    F = FKernel(fecundity, domain)
    return IntegralProjectionModels.IPMProblem(P + F, domain, n0, tspan; kwargs...)
end

export IPMModel, symbolic_ipm

end # module
