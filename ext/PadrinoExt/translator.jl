"""
R → Julia AST translation.

Converts parsed R expressions (RExpr) into Julia Expr objects, handling
PADRINO-specific distribution calls, state variable naming, and R→Julia
function mappings.
"""

# --- Distribution mapping ---

# PADRINO distribution name → (Julia Distributions.jl type, parameter names for PDF)
const PADRINO_DISTRIBUTIONS = Dict{String, Tuple{String, Vector{String}}}(
    "Norm"     => ("Normal", ["mu", "sd"]),
    "Lognorm"  => ("LogNormal", ["meanlog", "sdlog"]),
    "Gamma"    => ("Gamma", ["shape", "rate"]),
    "Beta"     => ("Beta", ["shape1", "shape2"]),
    "Expo"     => ("Exponential", ["rate"]),
    "Unif"     => ("Uniform", ["min", "max"]),
    "Cauchy"   => ("Cauchy", ["location", "scale"]),
    "Weib"     => ("Weibull", ["shape", "scale"]),
    "Pois"     => ("Poisson", ["lambda"]),
    "Binom"    => ("Binomial", ["size", "prob"]),
    "Bernoulli"=> ("Bernoulli", ["prob"]),
    "Negbin"   => ("NegativeBinomial", ["size", "prob"]),
)

# Truncated distribution mapping
const PADRINO_TRUNCATED_DISTRIBUTIONS = Dict{String, String}(
    "TNorm"    => "Normal",
    "TLognorm" => "LogNormal",
    "TGamma"   => "Gamma",
    "TBeta"    => "Beta",
    "TExpo"    => "Exponential",
    "TCauchy"  => "Cauchy",
    "TWeib"    => "Weibull",
)

# R function → Julia function mapping
const R_FUNCTION_MAP = Dict{String, String}(
    "exp"    => "exp",
    "log"    => "log",
    "log10"  => "log10",
    "sqrt"   => "sqrt",
    "abs"    => "abs",
    "sin"    => "sin",
    "cos"    => "cos",
    "tan"    => "tan",
    "asin"   => "asin",
    "acos"   => "acos",
    "atan"   => "atan",
    "ceiling"=> "ceil",
    "floor"  => "floor",
    "round"  => "round",
    "max"    => "_r_max",
    "min"    => "_r_min",
    "sum"    => "sum",
    "prod"   => "prod",
    "mean"   => "mean",
    "pi"     => "pi",
    "Inf"    => "Inf",
    "TRUE"   => "true",
    "FALSE"  => "false",
    "T"      => "true",
    "F"      => "false",
    "NA"     => "NaN",
)

# State variable suffix mapping: PADRINO uses _1 (source) and _2 (target)
const STATE_VAR_SUFFIXES = Dict{String, Symbol}(
    "_1" => :z,        # source state
    "_2" => :z_new,    # target state
)

"""
    translate_rexpr(expr::RExpr, state_vars::Vector{String}) -> Expr

Translate an R expression AST node into a Julia expression.
`state_vars` is a list of state variable base names (e.g., ["size", "leafarea"]).
"""
function translate_rexpr(expr::RNumber, state_vars::Vector{String})
    v = expr.value
    return isinteger(v) ? Int(v) : v
end

function translate_rexpr(expr::RString, state_vars::Vector{String})
    return expr.value
end

function translate_rexpr(expr::RIdent, state_vars::Vector{String})
    name = expr.name

    # Check for state variable references (size_1 → z, size_2 → z_new)
    for sv in state_vars
        if name == "$(sv)_1"
            return :z
        elseif name == "$(sv)_2"
            return :z_new
        elseif name == "d_$(sv)"
            return :h
        end
    end

    # R constants — but be careful not to shadow parameter names.
    # Use :pi (not :π) so that parameter assignment `pi = value` shadows correctly.
    if name == "TRUE" || name == "T"
        return true
    elseif name == "FALSE" || name == "F"
        return false
    elseif name == "NA" || name == "NA_real_"
        return NaN
    elseif name == "pi"
        return :pi  # keep as :pi so parameter shadowing works
    elseif name == "Inf"
        return Inf
    elseif name == "NaN"
        return NaN
    end

    return Symbol(name)
end

function translate_rexpr(expr::RBinOp, state_vars::Vector{String})
    left = translate_rexpr(expr.left, state_vars)
    right = translate_rexpr(expr.right, state_vars)

    op = Symbol(expr.op)
    return Expr(:call, op, left, right)
end

function translate_rexpr(expr::RUnaryOp, state_vars::Vector{String})
    operand = translate_rexpr(expr.operand, state_vars)
    if expr.op == "-"
        return Expr(:call, :(-), operand)
    elseif expr.op == "!"
        return Expr(:call, :(!), operand)
    end
    return operand
end

function translate_rexpr(expr::RCall, state_vars::Vector{String})
    func = expr.func
    args = [translate_rexpr(a, state_vars) for a in expr.args]

    # Handle PADRINO distribution calls (used in "Substituted" vital rates)
    if haskey(PADRINO_DISTRIBUTIONS, func)
        return _translate_distribution(func, expr, state_vars)
    end

    # Handle truncated distributions
    if haskey(PADRINO_TRUNCATED_DISTRIBUTIONS, func)
        return _translate_truncated_distribution(func, expr, state_vars)
    end

    # Handle specific R functions
    if func == "ifelse"
        length(args) >= 3 || error("ifelse requires 3 arguments")
        return Expr(:if, args[1], args[2], args[3])
    end

    if func == "pmin"
        return Expr(:call, :min, args...)
    end

    if func == "pmax"
        return Expr(:call, :max, args...)
    end

    if func == "is.na"
        return Expr(:call, :isnan, args...)
    end

    if func == "pnorm"
        # pnorm(x, mean, sd) → cdf(Normal(mean, sd), x)
        if length(args) >= 3
            return Expr(:call, :cdf, Expr(:call, :Normal, args[2], args[3]), args[1])
        elseif length(args) == 1
            return Expr(:call, :cdf, Expr(:call, :Normal, 0, 1), args[1])
        end
    end

    if func == "qnorm"
        # qnorm(p, mean, sd) → quantile(Normal(mean, sd), p)
        if length(args) >= 3
            return Expr(:call, :quantile, Expr(:call, :Normal, args[2], args[3]), args[1])
        elseif length(args) == 1
            return Expr(:call, :quantile, Expr(:call, :Normal, 0, 1), args[1])
        end
    end

    if func == "dnorm"
        # dnorm(x, mean, sd) → pdf(Normal(mean, sd), x)
        if length(args) >= 3
            return Expr(:call, :pdf, Expr(:call, :Normal, args[2], args[3]), args[1])
        elseif length(args) == 1
            return Expr(:call, :pdf, Expr(:call, :Normal, 0, 1), args[1])
        end
    end

    if func == "dlnorm"
        if length(args) >= 3
            return Expr(:call, :pdf, Expr(:call, :LogNormal, args[2], args[3]), args[1])
        end
    end

    if func == "dgamma"
        if length(args) >= 3
            return Expr(:call, :pdf, Expr(:call, :Gamma, args[2], args[3]), args[1])
        end
    end

    if func == "dbeta"
        if length(args) >= 3
            return Expr(:call, :pdf, Expr(:call, :Beta, args[2], args[3]), args[1])
        end
    end

    if func == "dexp"
        if length(args) >= 2
            return Expr(:call, :pdf, Expr(:call, :Exponential, Expr(:call, :/, 1, args[2])), args[1])
        end
    end

    if func == "dunif"
        if length(args) >= 3
            return Expr(:call, :pdf, Expr(:call, :Uniform, args[2], args[3]), args[1])
        end
    end

    if func == "c"
        # R's c() → Julia vector literal
        return Expr(:vect, args...)
    end

    if func == "seq"
        if length(args) == 2
            return Expr(:call, :(:), args[1], args[2])
        end
    end

    if func == "rep"
        if length(args) >= 2
            return Expr(:call, :fill, args[1], args[2])
        end
    end

    if func == "matrix"
        # R's matrix(x) → just x (no-op for scalars in pointwise context)
        if length(args) >= 1
            return args[1]
        end
    end

    if func == "kronecker"
        # Kronecker delta: kronecker(x, y) → 1.0 when |x - y| ≤ h/2, else 0.0
        # h is available in scope from the kernel function signature
        if length(args) >= 2
            return Expr(:if, Expr(:call, :<=, Expr(:call, :abs, Expr(:call, :-, args[1], args[2])),
                                  Expr(:call, :/, :h, 2)), 1.0, 0.0)
        end
    end

    # Generic function translation
    jl_func = get(R_FUNCTION_MAP, func, func)
    return Expr(:call, Symbol(jl_func), args...)
end

"""
    _translate_distribution(dist_name, call::RCall, state_vars) -> Expr

Translate a PADRINO distribution call (e.g., `Norm(mu_g, sd_g)`)
into a Julia `pdf(Distribution(...), z_new)` expression.
"""
function _translate_distribution(dist_name::String, call::RCall, state_vars::Vector{String})
    julia_dist, _param_names = PADRINO_DISTRIBUTIONS[dist_name]
    args = [translate_rexpr(a, state_vars) for a in call.args]
    kwargs = Dict(k => translate_rexpr(v, state_vars) for (k, v) in call.kwargs)

    # Check for first_arg=TRUE pattern: Lognorm(expr, m, s, first_arg=TRUE)
    has_first_arg = get(kwargs, "first_arg", nothing)
    if has_first_arg !== nothing
        delete!(kwargs, "first_arg")
        # The first positional arg is the evaluation point, rest are distribution params
        eval_point = args[1]
        dist_args = args[2:end]
        dist_expr = _make_dist_expr(julia_dist, dist_name, dist_args)
        return Expr(:call, :pdf, dist_expr, eval_point)
    end

    # Standard case: distribution evaluated at z_new
    dist_expr = _make_dist_expr(julia_dist, dist_name, args)
    return Expr(:call, :pdf, dist_expr, :z_new)
end

"""
    _make_dist_expr(julia_dist, padrino_name, args) -> Expr

Create a Distributions.jl constructor expression with proper parameterization.
"""
function _make_dist_expr(julia_dist::String, padrino_name::String, args)
    dist_sym = Symbol(julia_dist)

    if padrino_name == "Expo"
        # R's Exponential uses rate, Julia uses scale = 1/rate
        length(args) >= 1 || error("Exponential requires at least 1 argument")
        return Expr(:call, dist_sym, Expr(:call, :/, 1, args[1]))
    end

    if padrino_name == "Gamma" && length(args) >= 2
        # R Gamma(shape, rate) → Julia Gamma(shape, 1/rate)
        return Expr(:call, dist_sym, args[1], Expr(:call, :/, 1, args[2]))
    end

    return Expr(:call, dist_sym, args...)
end

"""
    _translate_truncated_distribution(dist_name, call::RCall, state_vars) -> Expr

Translate a truncated PADRINO distribution call.
"""
function _translate_truncated_distribution(dist_name::String, call::RCall, state_vars::Vector{String})
    base_dist_name = PADRINO_TRUNCATED_DISTRIBUTIONS[dist_name]
    args = [translate_rexpr(a, state_vars) for a in call.args]

    # Truncated distributions: params..., a, b (last two are truncation bounds)
    if length(args) < 3
        error("Truncated distribution $dist_name requires at least 3 arguments (params + bounds)")
    end

    # Last two args are truncation bounds
    lower_bound = args[end-1]
    upper_bound = args[end]
    dist_args = args[1:end-2]

    base_expr = _make_dist_expr(base_dist_name, _inv_trunc_name(dist_name), dist_args)
    trunc_expr = Expr(:call, :Truncated, base_expr, lower_bound, upper_bound)
    return Expr(:call, :pdf, trunc_expr, :z_new)
end

function _inv_trunc_name(tname::String)
    # TNorm → Norm, TLognorm → Lognorm, etc.
    return tname[2:end]
end

# --- Random variate generation for environmental variables ---

# PADRINO distribution name → (Julia rand function expression)
const PADRINO_RAND_DISTRIBUTIONS = Dict{String, Tuple{String, Vector{String}}}(
    "Norm"     => ("Normal", ["mean", "sd"]),
    "Lognorm"  => ("LogNormal", ["meanlog", "sdlog"]),
    "Gamma"    => ("Gamma", ["shape", "rate"]),
    "Beta"     => ("Beta", ["shape1", "shape2"]),
    "Expo"     => ("Exponential", ["rate"]),
    "Unif"     => ("Uniform", ["min", "max"]),
    "Cauchy"   => ("Cauchy", ["location", "scale"]),
    "Weib"     => ("Weibull", ["shape", "scale"]),
    "Pois"     => ("Poisson", ["lambda"]),
    "Binom"    => ("Binomial", ["size", "prob"]),
    "sample"   => ("_sample", ["x"]),
    "TNorm"    => ("Normal", ["mean", "sd", "a", "b"]),
    "TLognorm" => ("LogNormal", ["meanlog", "sdlog", "a", "b"]),
)

"""
    translate_env_function(env_func_str, env_vars_data) -> Expr

Translate a PADRINO environmental variable sampling function into a Julia expression
that generates a random variate.
"""
function translate_env_function(env_func_str::String, state_vars::Vector{String})
    rexpr = parse_rexpr(env_func_str)
    return _translate_env_rexpr(rexpr, state_vars)
end

function _translate_env_rexpr(expr::RCall, state_vars::Vector{String})
    func = expr.func
    args = [translate_rexpr(a, state_vars) for a in expr.args]

    # PADRINO sampling distributions → rand(Distribution(...))
    if haskey(PADRINO_RAND_DISTRIBUTIONS, func)
        julia_dist, _pnames = PADRINO_RAND_DISTRIBUTIONS[func]

        if startswith(func, "T") && func != "T_dist"
            # Truncated: last two args are bounds
            if length(args) >= 3
                dist_args = args[1:end-2]
                lower = args[end-1]
                upper = args[end]
                base = _make_dist_expr(julia_dist, _inv_trunc_name(func), dist_args)
                dist_expr = Expr(:call, :Truncated, base, lower, upper)
                return Expr(:call, :rand, dist_expr)
            end
        end

        if func == "Expo"
            dist_expr = Expr(:call, Symbol(julia_dist), Expr(:call, :/, 1, args[1]))
        elseif func == "Gamma" && length(args) >= 2
            dist_expr = Expr(:call, Symbol(julia_dist), args[1], Expr(:call, :/, 1, args[2]))
        elseif func == "sample"
            return Expr(:call, :rand, args[1])
        else
            dist_expr = Expr(:call, Symbol(julia_dist), args...)
        end

        return Expr(:call, :rand, dist_expr)
    end

    # Fall back to generic translation
    return translate_rexpr(expr, state_vars)
end

function _translate_env_rexpr(expr::RExpr, state_vars::Vector{String})
    return translate_rexpr(expr, state_vars)
end
