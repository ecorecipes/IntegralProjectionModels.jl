"""
Vital rate compilation and kernel matrix construction.

Compiles PADRINO vital rate expressions into Julia functions, resolves
dependency order, and builds discretized kernel matrices.
"""

using Distributions
using Statistics: mean

# --- Interpreted vectorized evaluation ---
#
# Instead of compiling unique anonymous functions via Base.eval (which triggers
# expensive JIT compilation per model), we interpret the scalar expression tree
# directly, broadcasting every function call. This matches R's approach: walk
# the AST and dispatch to vectorized C primitives.

"""Lookup table for functions that can appear in translated PADRINO expressions."""
const _EVAL_FUNCS = let d = Dict{Symbol, Any}(
    # Arithmetic (handled as broadcast)
    :+ => +, :- => -, :* => *, :/ => /, :^ => ^,
    # Math
    :exp => exp, :log => log, :log10 => log10, :log2 => log2,
    :sqrt => sqrt, :abs => abs, :ceil => ceil, :floor => floor, :round => round,
    :sin => sin, :cos => cos, :tan => tan, :asin => asin, :acos => acos, :atan => atan,
    # Comparison
    :< => <, :> => >, :(<=) => <=, :(>=) => >=, :(==) => ==, :(!=) => !=,
    # Logic
    :& => &, :| => |, :! => !,
    # Element-wise min/max (from R's pmin/pmax)
    :max => max, :min => min,
    # Control
    :ifelse => ifelse,
    # Distributions
    :pdf => Distributions.pdf, :cdf => Distributions.cdf, :quantile => Distributions.quantile,
    :Normal => Distributions.Normal, :LogNormal => Distributions.LogNormal,
    :Gamma => Distributions.Gamma, :Beta => Distributions.Beta,
    :Exponential => Distributions.Exponential, :Uniform => Distributions.Uniform,
    :Cauchy => Distributions.Cauchy, :Weibull => Distributions.Weibull,
    :Poisson => Distributions.Poisson, :Binomial => Distributions.Binomial,
    :Bernoulli => Distributions.Bernoulli, :NegativeBinomial => Distributions.NegativeBinomial,
    :Truncated => Distributions.Truncated,
    # Type conversion
    :Float64 => Float64,
    # Aggregation (these should be pre-computed, but just in case)
    :sum => sum, :prod => prod, :mean => mean,
    # NaN check
    :isnan => isnan,
); d end

"""
    _eval_vectorized(expr, env::Dict{Symbol, Any})

Interpret a scalar Julia expression with automatic broadcasting.
Every function call is broadcast over array arguments, matching R's
vectorized evaluation semantics. This avoids Base.eval and JIT compilation.
"""
function _eval_vectorized(expr::Expr, env::Dict{Symbol, Any})
    if expr.head == :call
        func_sym = expr.args[1]
        if func_sym isa Symbol && haskey(_EVAL_FUNCS, func_sym)
            f = _EVAL_FUNCS[func_sym]
            args = Any[_eval_vectorized(a, env) for a in @view expr.args[2:end]]
            return broadcast(f, args...)
        else
            error("Unknown function in kernel expression: $func_sym")
        end
    elseif expr.head == :(=)
        lhs = expr.args[1]::Symbol
        rhs = _eval_vectorized(expr.args[2], env)
        env[lhs] = rhs
        return rhs
    elseif expr.head == :block
        result = nothing
        for arg in expr.args
            arg isa LineNumberNode && continue
            result = _eval_vectorized(arg, env)
        end
        return result
    elseif expr.head == :return
        return _eval_vectorized(expr.args[1], env)
    elseif expr.head == :if && length(expr.args) >= 3
        cond = _eval_vectorized(expr.args[1], env)
        t_val = _eval_vectorized(expr.args[2], env)
        f_val = _eval_vectorized(expr.args[3], env)
        return broadcast(ifelse, cond, t_val, f_val)
    elseif expr.head == :comparison
        # Handle chained comparisons like a <= b
        # Julia AST: Expr(:comparison, a, :(<=), b)
        left = _eval_vectorized(expr.args[1], env)
        op = expr.args[2]
        right = _eval_vectorized(expr.args[3], env)
        if haskey(_EVAL_FUNCS, op)
            return broadcast(_EVAL_FUNCS[op], left, right)
        end
        error("Unknown comparison operator: $op")
    else
        error("Unsupported expression head in kernel: $(expr.head)")
    end
end

_eval_vectorized(sym::Symbol, env::Dict{Symbol, Any}) = env[sym]
_eval_vectorized(x::Number, env::Dict{Symbol, Any}) = x
_eval_vectorized(::LineNumberNode, env::Dict{Symbol, Any}) = nothing
_eval_vectorized(x::Bool, env::Dict{Symbol, Any}) = x

"""
    _sanitize_name(name::AbstractString) -> String

Replace dots in R-style names with underscores for valid Julia identifiers.
E.g., "grow.int_1400_current" → "grow_int_1400_current"
"""
_sanitize_name(name::AbstractString) = replace(string(name), "." => "_")

"""
    _sanitize_formula(formula::AbstractString) -> String

Replace dots in R-style identifier names with underscores, but preserve dots
in numeric literals (e.g., 0.15 stays as 0.15, but grow.int becomes grow_int).
Also preserves R dot-functions by converting them to their standard forms first.
"""
function _sanitize_formula(formula::AbstractString)
    result = string(formula)
    # Pre-process known R dot-functions to protect them from sanitization
    # is.na → isna_RFUNC (temporary marker), then after general sanitization → is.na
    # Actually, simpler: replace with the Julia equivalent directly
    result = replace(result, "is.na(" => "is_na_RFUNC(")
    # Replace dots that are between identifier characters (letters/digits/underscores)
    # but NOT dots that are part of numbers (preceded by digit and followed by digit)
    result = replace(result, r"(?<=[A-Za-z_])\.(?=[A-Za-z_])" => "_")
    # Restore R dot-functions
    result = replace(result, "is_na_RFUNC(" => "is.na(")
    return result
end

"""
    compile_vital_rates(vr_exprs, state_vars) -> Dict{String, Any}

Compile a collection of vital rate expressions into Julia Expr objects,
sorted in dependency order.
"""
function compile_vital_rates(vr_exprs::Vector{Dict{Symbol, Any}},
                             state_vars::Vector{String})
    compiled = Dict{String, Any}()

    for vr in vr_exprs
        formula = string(vr[:formula])
        model_type = string(vr[:model_type])
        name = _sanitize_name(vr[:name])

        # Skip self-assignments like "s = s"
        _is_self_assignment(formula) && continue

        # Sanitize dots in the formula (R names like grow.int → grow_int)
        formula = _sanitize_formula(formula)

        # Parse the right-hand side
        parts = split(formula, "="; limit=2)
        length(parts) == 2 || continue
        rhs = strip(parts[2])

        if model_type == "Evaluated"
            rexpr = parse_rexpr(rhs)
            compiled[name] = translate_rexpr(rexpr, state_vars)
        elseif model_type == "Substituted"
            # Distribution call — translate as PDF
            rexpr = parse_rexpr(rhs)
            compiled[name] = translate_rexpr(rexpr, state_vars)
        end
    end

    return compiled
end

"""
    apply_truncation_to_vrs!(compiled_vrs::Dict{String, Any}, lower, upper)

For discrete_extrema eviction: wrap distribution calls in `Truncated(dist, L, U)`
so that distributions integrate to 1 within the domain bounds.
This achieves the same effect as ipmr's discrete_extrema correction.
"""
function apply_truncation_to_vrs!(compiled_vrs::Dict{String, Any}, lower::Float64, upper::Float64)
    for (name, expr) in compiled_vrs
        compiled_vrs[name] = _wrap_distributions_truncated(expr, lower, upper)
    end
end

const _DISTRIBUTION_NAMES = Set([:Normal, :LogNormal, :Gamma, :Beta, :Exponential,
                                  :Uniform, :Cauchy, :Weibull, :Poisson, :Binomial,
                                  :Bernoulli, :NegativeBinomial])

function _wrap_distributions_truncated(expr::Expr, lower::Float64, upper::Float64)
    # Look for pattern: pdf(Distribution(...), eval_point)
    if expr.head == :call && length(expr.args) >= 3 && expr.args[1] == :pdf
        dist_arg = expr.args[2]
        if dist_arg isa Expr && dist_arg.head == :call && dist_arg.args[1] in _DISTRIBUTION_NAMES
            # Already truncated? Skip
            if !(dist_arg.args[1] == :Truncated)
                # Wrap: pdf(Truncated(dist, L, U), eval_point)
                trunc_expr = Expr(:call, :Truncated, dist_arg, lower, upper)
                return Expr(:call, :pdf, trunc_expr, expr.args[3:end]...)
            end
        end
        # Check if already Truncated
        if dist_arg isa Expr && dist_arg.head == :call && dist_arg.args[1] == :Truncated
            return expr  # already truncated
        end
    end

    # Recurse into sub-expressions
    new_args = Any[]
    for arg in expr.args
        push!(new_args, _wrap_distributions_truncated(arg, lower, upper))
    end
    return Expr(expr.head, new_args...)
end

function _wrap_distributions_truncated(x, lower::Float64, upper::Float64)
    return x  # non-Expr nodes pass through
end

"""
    _is_self_assignment(formula::String) -> Bool

Check if a formula is a trivial self-assignment like "s = s".
"""
function _is_self_assignment(formula::AbstractString)
    parts = split(formula, "="; limit=2)
    length(parts) == 2 || return false
    return strip(parts[1]) == strip(parts[2])
end

"""
    resolve_dependency_order(compiled_vrs::Dict{String, Any}) -> Vector{String}

Topologically sort vital rate expressions by dependency.
"""
function resolve_dependency_order(compiled_vrs::Dict{String, Any})
    names = collect(keys(compiled_vrs))
    name_set = Set(names)

    # Build dependency graph
    deps = Dict{String, Set{String}}()
    for name in names
        deps[name] = _find_dependencies(compiled_vrs[name], name_set)
    end

    # Topological sort (Kahn's algorithm)
    sorted = String[]
    remaining = Set(names)
    resolved = Set{String}()

    while !isempty(remaining)
        # Find nodes with no unresolved dependencies
        ready = [n for n in remaining if isempty(setdiff(deps[n], resolved))]

        if isempty(ready)
            # Circular dependency — just add remaining in arbitrary order
            append!(sorted, collect(remaining))
            break
        end

        for n in ready
            push!(sorted, n)
            push!(resolved, n)
            delete!(remaining, n)
        end
    end

    return sorted
end

"""
    _find_dependencies(expr, name_set::Set{String}) -> Set{String}

Find which vital rate names appear in an expression.
"""
function _find_dependencies(expr::Expr, name_set::Set{String})
    deps = Set{String}()
    _walk_expr!(deps, expr, name_set)
    return deps
end

function _find_dependencies(sym::Symbol, name_set::Set{String})
    s = String(sym)
    return s in name_set ? Set{String}([s]) : Set{String}()
end

function _find_dependencies(expr, name_set::Set{String})
    return Set{String}()
end

function _walk_expr!(deps::Set{String}, expr::Expr, name_set::Set{String})
    if expr.head == :call
        for arg in expr.args[2:end]
            _walk_expr!(deps, arg, name_set)
        end
    else
        for arg in expr.args
            _walk_expr!(deps, arg, name_set)
        end
    end
end

function _walk_expr!(deps::Set{String}, sym::Symbol, name_set::Set{String})
    s = String(sym)
    if s in name_set
        push!(deps, s)
    end
end

function _walk_expr!(deps::Set{String}, ::Any, name_set::Set{String})
    # Numbers, strings, etc.
end

"""
    _find_reachable_vrs(kernel_expr, vital_rates) -> Set{String}

Find all VR names transitively referenced by the kernel expression.
This implements dead code elimination — VRs not reachable from the kernel
expression are excluded, mimicking R's lazy evaluation behavior.
"""
function _find_reachable_vrs(kernel_expr, vital_rates::Dict{String, Any})
    name_set = Set(keys(vital_rates))
    reachable = Set{String}()

    # Start with VRs directly referenced by the kernel expression
    frontier = _find_dependencies(kernel_expr, name_set)

    while !isempty(frontier)
        new_frontier = Set{String}()
        for name in frontier
            if name ∉ reachable && haskey(vital_rates, name)
                push!(reachable, name)
                # Find VRs referenced by this VR's expression
                deps = _find_dependencies(vital_rates[name], name_set)
                for d in deps
                    if d ∉ reachable
                        push!(new_frontier, d)
                    end
                end
            end
        end
        frontier = new_frontier
    end

    return reachable
end

"""
    _dotify_expr(expr)

Transform a scalar Julia expression into a broadcasted version.
Converts `f(args...)` to `f.(args...)` and arithmetic operators to their
broadcast equivalents, enabling evaluation over meshpoint arrays.
"""
function _dotify_expr(expr::Expr)
    if expr.head == :call
        # f(args...) → f.(dotified_args...)
        func = expr.args[1]
        args = map(_dotify_expr, expr.args[2:end])
        return Expr(:., func, Expr(:tuple, args...))
    elseif expr.head == :(=)
        # x = rhs → x = dotified_rhs (only dotify RHS)
        return Expr(:(=), expr.args[1], _dotify_expr(expr.args[2]))
    elseif expr.head == :return
        return Expr(:return, _dotify_expr(expr.args[1]))
    elseif expr.head == :block
        return Expr(:block, map(_dotify_expr, expr.args)...)
    elseif expr.head == :if && length(expr.args) == 3
        # if/else → ifelse.() for broadcasting
        return Expr(:., :ifelse, Expr(:tuple,
            _dotify_expr(expr.args[1]),
            _dotify_expr(expr.args[2]),
            _dotify_expr(expr.args[3])))
    else
        return Expr(expr.head, map(_dotify_expr, expr.args)...)
    end
end

_dotify_expr(x) = x  # Numbers, Symbols, LineNumberNodes pass through

"""
    _build_kernel_body(kernel_formula_rhs, vital_rates, state_vars, params) -> Expr

Internal: build the shared function body (parameter assignments + VR chain + kernel return)
from a kernel formula. This is the common work shared by scalar and vectorized compilation.
"""
function _build_kernel_body(kernel_formula_rhs::AbstractString,
                            vital_rates::Dict{String, Any},
                            state_vars::Vector{String},
                            params::Dict{String, Float64})
    kernel_expr = translate_rexpr(parse_rexpr(kernel_formula_rhs), state_vars)

    reachable = _find_reachable_vrs(kernel_expr, vital_rates)
    reachable_vrs = Dict{String, Any}(k => v for (k, v) in vital_rates if k in reachable)
    order = resolve_dependency_order(reachable_vrs)

    body_exprs = Expr[]
    for vr_name in order
        if haskey(vital_rates, vr_name)
            push!(body_exprs, :($(Symbol(vr_name)) = $(vital_rates[vr_name])))
        end
    end
    push!(body_exprs, :(return $kernel_expr))

    param_assigns = Expr[]
    for (pname, pval) in params
        push!(param_assigns, :($(Symbol(_sanitize_name(pname))) = $pval))
    end

    return Expr(:block, param_assigns..., body_exprs...)
end

"""
    build_kernel_function(kernel_formula_rhs, vital_rates, state_vars, params) -> Function

Compile a kernel expression and its vital rate chain into a scalar function
`(z_new, z, h) -> kernel_value`.
"""
function build_kernel_function(kernel_formula_rhs::AbstractString,
                               vital_rates::Dict{String, Any},
                               state_vars::Vector{String},
                               params::Dict{String, Float64})
    func_body = _build_kernel_body(kernel_formula_rhs, vital_rates, state_vars, params)
    func_expr = :((z_new, z, h) -> $func_body)

    try
        return Base.eval(PadrinoEvalModule, func_expr)
    catch e
        @warn "Failed to compile kernel function" exception=e kernel_formula_rhs
        return (z_new, z, h) -> 0.0
    end
end

"""
    build_vectorized_kernel_function(kernel_formula_rhs, vital_rates, state_vars, params) -> Function or nothing

Compile a broadcasted kernel function `(z_new_col, z_row, h) -> Matrix{Float64}`.
Returns `nothing` on failure.
"""
function build_vectorized_kernel_function(kernel_formula_rhs::AbstractString,
                                          vital_rates::Dict{String, Any},
                                          state_vars::Vector{String},
                                          params::Dict{String, Float64})
    func_body = _build_kernel_body(kernel_formula_rhs, vital_rates, state_vars, params)
    dotified_body = _dotify_expr(func_body)
    func_expr = :((z_new, z, h) -> $dotified_body)

    try
        return Base.eval(PadrinoEvalModule, func_expr)
    catch e
        @warn "Failed to compile vectorized kernel function" exception=e
        return nothing
    end
end

"""
    build_kernel_functions(kernel_formula_rhs, vital_rates, state_vars, params) -> (body, nothing)

Build kernel body expression from a single parse pass. Returns `(body::Expr, nothing)`.

No compilation occurs here — the body is evaluated via `_eval_vectorized` (interpreted
broadcasting) in `build_kernel_matrix`, avoiding Base.eval and JIT overhead entirely.
Falls back to compiled evaluation only if interpretation fails.
"""
function build_kernel_functions(kernel_formula_rhs::AbstractString,
                                vital_rates::Dict{String, Any},
                                state_vars::Vector{String},
                                params::Dict{String, Float64})
    func_body = _build_kernel_body(kernel_formula_rhs, vital_rates, state_vars, params)
    return func_body, nothing
end

"""
    _ensure_compiled(fn) -> Function

Compile a deferred kernel function body (Expr) into a callable function.
If `fn` is already a Function, returns it unchanged.
"""
function _ensure_compiled(fn)
    fn isa Expr || return fn
    try
        return Base.eval(PadrinoEvalModule, :((z_new, z, h) -> $fn))
    catch e
        @warn "Failed to compile deferred kernel function" exception=e
        return (z_new, z, h) -> 0.0
    end
end

"""
    build_kernel_function_with_env(kernel_formula_rhs::AbstractString,
                                   vital_rates::Dict{String, Any},
                                   state_vars::Vector{String},
                                   base_params::Dict{String, Float64}) -> Function

Like `build_kernel_function` but returns a function that also accepts
a parameter dict: `(z_new, z, h, params) -> kernel_value`.
"""
function build_kernel_function_with_env(kernel_formula_rhs::AbstractString,
                                        vital_rates::Dict{String, Any},
                                        state_vars::Vector{String},
                                        base_params::Dict{String, Float64})
    kernel_expr = translate_rexpr(parse_rexpr(kernel_formula_rhs), state_vars)
    order = resolve_dependency_order(vital_rates)

    body_exprs = Expr[]
    for vr_name in order
        if haskey(vital_rates, vr_name)
            push!(body_exprs, :($(Symbol(vr_name)) = $(vital_rates[vr_name])))
        end
    end
    push!(body_exprs, :(return $kernel_expr))

    # Parameter assignments from the merged params dict (sanitize dots in names)
    all_param_names = collect(keys(base_params))
    param_assigns = [:($(Symbol(_sanitize_name(pn))) = params[$(pn)]) for pn in all_param_names]

    func_body = Expr(:block, param_assigns..., body_exprs...)
    func_expr = :((z_new, z, h, params) -> $func_body)

    try
        return Base.eval(PadrinoEvalModule, func_expr)
    catch e
        @warn "Failed to compile kernel function with env" exception=e
        return (z_new, z, h, params) -> 0.0
    end
end

"""
    build_kernel_matrix(kernel_fn, domain::ContinuousDomain, family::AbstractString) -> Matrix{Float64}

Discretize a kernel function into an m×m matrix using the midpoint rule.
For CC kernels, includes the step size h factor.
For CD/DC/DD kernels, h is not included.
"""
function build_kernel_matrix(kernel_fn, domain::IPM.ContinuousDomain, family::AbstractString)
    z = IPM.meshpoints(domain)
    h = IPM.step_size(domain)
    m = length(z)

    include_h = (family == "CC")

    # Scalar fallback using element-wise invokelatest
    K = zeros(Float64, m, m)
    @inbounds for j in 1:m
        zj = z[j]
        for i in 1:m
            val = Base.invokelatest(kernel_fn, z[i], zj, h)::Float64
            K[i, j] = include_h ? val * h : val
        end
    end

    return K
end

"""
    build_kernel_matrix(kernel_fn, vectorized_fn, domain, family) -> Matrix{Float64}

Build kernel matrix. Evaluation strategy (in order of preference):
1. Interpreted vectorized evaluation (if kernel_fn is an Expr body) — no JIT
2. Compiled vectorized evaluation (if vectorized_fn is provided)
3. Compiled scalar evaluation (loop with invokelatest)
"""
function build_kernel_matrix(kernel_fn, vectorized_fn, domain::IPM.ContinuousDomain, family::AbstractString)
    z = IPM.meshpoints(domain)
    h = IPM.step_size(domain)
    m = length(z)

    include_h = (family == "CC")

    # Path 1: Interpreted vectorized evaluation (no compilation, no JIT)
    if kernel_fn isa Expr
        try
            z_new = reshape(z, :, 1)  # m×1 column (target states)
            z_src = reshape(z, 1, :)  # 1×m row (source states)
            env = Dict{Symbol, Any}(:z_new => z_new, :z => z_src, :h => h)
            K = _eval_vectorized(kernel_fn, env)
            K = convert(Matrix{Float64}, K)
            if size(K) == (m, m)
                if include_h
                    K .*= h
                end
                return K
            end
        catch
            # Fall through to compiled paths
        end
    end

    # Path 2: Compiled vectorized evaluation
    if vectorized_fn !== nothing
        try
            z_new = reshape(z, :, 1)
            z_src = reshape(z, 1, :)
            K = Base.invokelatest(vectorized_fn, z_new, z_src, h)
            K = convert(Matrix{Float64}, K)
            if size(K) == (m, m)
                if include_h
                    K .*= h
                end
                return K
            end
        catch
        end
    end

    # Path 3: Compiled scalar evaluation
    return build_kernel_matrix(_ensure_compiled(kernel_fn), domain, family)
end

"""
    build_kernel_matrix(kernel_fn, domain::IPM.ContinuousDomain, family::AbstractString, params::Dict) -> Matrix{Float64}

Build kernel matrix with dynamic parameters (for stochastic models).
"""
function build_kernel_matrix(kernel_fn, domain::IPM.ContinuousDomain, family::AbstractString, params::Dict)
    z = IPM.meshpoints(domain)
    h = IPM.step_size(domain)
    m = length(z)
    K = zeros(Float64, m, m)

    include_h = (family == "CC")

    @inbounds for j in 1:m
        zj = z[j]
        for i in 1:m
            val = Base.invokelatest(kernel_fn, z[i], zj, h, params)::Float64
            K[i, j] = include_h ? val * h : val
        end
    end

    return K
end

"""
    precompute_mesh_sums!(compiled_vrs::Dict{String, Any}, state_vars::Vector{String},
                           params::Dict{String, Float64}, domain::IPM.ContinuousDomain)

Detect `sum(vr_name)` calls in vital rate expressions and replace them with
pre-computed constants (sum of vr over all meshpoints).

In ipmr, `sum(r_s)` means: compute r_s at each meshpoint and sum them.
This is a mesh-level aggregation that must be pre-computed before the kernel
function is built.
"""
function precompute_mesh_sums!(compiled_vrs::Dict{String, Any},
                                state_vars::Vector{String},
                                params::Dict{String, Float64},
                                domain::IPM.ContinuousDomain)
    z = IPM.meshpoints(domain)
    m = length(z)
    h = IPM.step_size(domain)

    # Find all sum(X) calls in the compiled expressions
    sum_targets = Dict{Symbol, Float64}()

    for (name, expr) in compiled_vrs
        _find_sum_calls!(sum_targets, expr, compiled_vrs, state_vars, params, z, h, m)
    end

    if isempty(sum_targets)
        return
    end

    # Replace sum(X) calls with the computed constants
    for (name, expr) in compiled_vrs
        compiled_vrs[name] = _replace_sum_calls(expr, sum_targets)
    end
end

function _find_sum_calls!(sum_targets::Dict{Symbol, Float64}, expr::Expr,
                           compiled_vrs::Dict{String, Any},
                           state_vars::Vector{String},
                           params::Dict{String, Float64},
                           z::Vector{Float64}, h::Float64, m::Int)
    if expr.head == :call && expr.args[1] == :sum && length(expr.args) == 2
        target = expr.args[2]
        if target isa Symbol && haskey(compiled_vrs, String(target))
            # Pre-compute sum of this VR over the mesh
            vr_name = String(target)
            if !haskey(sum_targets, target)
                sum_val = _compute_mesh_sum(vr_name, compiled_vrs, state_vars, params, z, h, m)
                sum_targets[target] = sum_val
            end
        elseif target isa Symbol && (target == :z || target == :z_new)
            # sum(state_variable) — sum of all meshpoint values
            if !haskey(sum_targets, target)
                sum_targets[target] = sum(z)
            end
        end
    end

    # Recurse into sub-expressions
    for arg in expr.args
        _find_sum_calls!(sum_targets, arg, compiled_vrs, state_vars, params, z, h, m)
    end
end

function _find_sum_calls!(::Dict{Symbol, Float64}, ::Any, args...)
    # Non-Expr nodes: nothing to find
end

function _compute_mesh_sum(vr_name::String, compiled_vrs::Dict{String, Any},
                            state_vars::Vector{String}, params::Dict{String, Float64},
                            z::Vector{Float64}, h::Float64, m::Int)
    # Build a function that computes just this VR and its dependencies
    deps_vrs = Dict{String, Any}(vr_name => compiled_vrs[vr_name])
    # Include dependencies
    order = resolve_dependency_order(compiled_vrs)
    for dep_name in order
        if haskey(compiled_vrs, dep_name) && dep_name != vr_name
            if _expr_references_vr(compiled_vrs[vr_name], dep_name, compiled_vrs)
                deps_vrs[dep_name] = compiled_vrs[dep_name]
            end
        end
    end

    vr_fn = build_kernel_function(vr_name, deps_vrs, state_vars, params)

    # Sum over all meshpoints (using source state z)
    total = 0.0
    for j in 1:m
        total += Base.invokelatest(vr_fn, z[j], z[j], h)
    end
    return total
end

"""Check if an expression transitively references a VR name"""
function _expr_references_vr(expr, target_name::String, compiled_vrs::Dict{String, Any})
    if _expr_references(expr, Symbol(target_name))
        return true
    end
    # Check transitive: if expr references X and X references target_name
    for (name, vr_expr) in compiled_vrs
        if name != target_name && _expr_references(expr, Symbol(name))
            if _expr_references_vr(vr_expr, target_name, compiled_vrs)
                return true
            end
        end
    end
    return false
end

function _replace_sum_calls(expr::Expr, sum_targets::Dict{Symbol, Float64})
    if expr.head == :call && expr.args[1] == :sum && length(expr.args) == 2
        target = expr.args[2]
        if target isa Symbol && haskey(sum_targets, target)
            return sum_targets[target]
        end
    end

    new_args = Any[]
    for arg in expr.args
        push!(new_args, _replace_sum_calls(arg, sum_targets))
    end
    return Expr(expr.head, new_args...)
end

function _replace_sum_calls(x, sum_targets::Dict{Symbol, Float64})
    return x
end

# --- R-style max()/min() reduction pre-computation ---

"""
    precompute_r_reductions!(compiled_vrs, state_vars, params, domain)

Handle R's `max()` and `min()` as mesh-level reductions.

In R, `max(a, b)` is a global reduction (returns single scalar = maximum of
all vector elements), unlike `pmax(a, b)` which is element-wise. When evaluating
kernel functions pointwise in Julia, this distinction matters:
- R: `max(0.00001, sigma_vec)` → scalar (global max across all meshpoints)
- Julia pointwise: `max(0.00001, sigma_i)` → varies per meshpoint

This function detects `_r_max()`/`_r_min()` calls (translated from R's
`max()`/`min()`), evaluates them as mesh-level reductions when their arguments
are state-dependent, and replaces them with pre-computed constants.
"""
function precompute_r_reductions!(compiled_vrs::Dict{String, Any},
                                   state_vars::Vector{String},
                                   params::Dict{String, Float64},
                                   domain::IPM.ContinuousDomain)
    z = IPM.meshpoints(domain)
    h = IPM.step_size(domain)
    m = length(z)

    for (name, expr) in compiled_vrs
        new_expr = _replace_r_reduction(expr, compiled_vrs, state_vars, params, z, h, m)
        compiled_vrs[name] = new_expr
    end
end

function _replace_r_reduction(expr::Expr, compiled_vrs, state_vars, params, z, h, m)
    if expr.head == :call && expr.args[1] in (:_r_max, :_r_min)
        fn = expr.args[1]

        # First, recursively process arguments
        processed_args = Any[fn]
        for arg in expr.args[2:end]
            push!(processed_args, _replace_r_reduction(arg, compiled_vrs, state_vars, params, z, h, m))
        end
        processed_expr = Expr(:call, processed_args...)

        # Check if any argument is state-dependent
        if _is_state_dependent(processed_expr, compiled_vrs)
            # Pre-compute as mesh-level reduction
            julia_fn = fn == :_r_max ? :max : :min
            eval_expr = Expr(:call, julia_fn, processed_args[2:end]...)
            reduction = fn == :_r_max ? maximum : minimum
            return _eval_reduction_over_mesh(eval_expr, compiled_vrs, state_vars,
                                              params, z, h, m, reduction)
        else
            # Not state-dependent: just use regular max/min
            julia_fn = fn == :_r_max ? :max : :min
            return Expr(:call, julia_fn, processed_args[2:end]...)
        end
    end

    # Recurse into sub-expressions
    new_args = Any[]
    for arg in expr.args
        push!(new_args, _replace_r_reduction(arg, compiled_vrs, state_vars, params, z, h, m))
    end
    return Expr(expr.head, new_args...)
end

function _replace_r_reduction(x, compiled_vrs, state_vars, params, z, h, m)
    return x
end

"""Check if an expression depends on state variables (z, z_new) transitively through VRs."""
function _is_state_dependent(expr, compiled_vrs::Dict{String, Any})
    return _is_sd(expr, compiled_vrs, Set{String}())
end

function _is_sd(expr::Expr, compiled_vrs::Dict{String, Any}, seen::Set{String})
    if _expr_references(expr, :z) || _expr_references(expr, :z_new)
        return true
    end
    vr_names = Set(keys(compiled_vrs))
    deps = _find_dependencies(expr, vr_names)
    for dep in deps
        dep in seen && continue
        push!(seen, dep)
        if haskey(compiled_vrs, dep) && _is_sd(compiled_vrs[dep], compiled_vrs, seen)
            return true
        end
    end
    return false
end

function _is_sd(sym::Symbol, compiled_vrs::Dict{String, Any}, seen::Set{String})
    sym in (:z, :z_new) && return true
    s = String(sym)
    s in seen && return false
    push!(seen, s)
    if haskey(compiled_vrs, s)
        return _is_sd(compiled_vrs[s], compiled_vrs, seen)
    end
    return false
end

function _is_sd(::Any, ::Dict{String, Any}, ::Set{String})
    return false
end

"""Evaluate an expression over all meshpoints and apply a reduction function."""
function _eval_reduction_over_mesh(eval_expr, compiled_vrs, state_vars, params,
                                     z, h, m, reduction_fn)
    # Find reachable VRs
    reachable = _find_reachable_vrs(eval_expr, compiled_vrs)
    reachable_vrs = Dict{String, Any}(k => v for (k, v) in compiled_vrs if k in reachable)
    order = resolve_dependency_order(reachable_vrs)

    body_exprs = Expr[]
    for vr_name in order
        if haskey(compiled_vrs, vr_name)
            # Replace any remaining _r_max/_r_min in VR expressions with max/min
            vr_expr = _normalize_r_fns(compiled_vrs[vr_name])
            push!(body_exprs, :($(Symbol(vr_name)) = $vr_expr))
        end
    end
    push!(body_exprs, :(return $eval_expr))

    param_assigns = [:($(Symbol(_sanitize_name(pn))) = $pval) for (pn, pval) in params]
    func_body = Expr(:block, param_assigns..., body_exprs...)
    func_expr = :((z_new, z, h) -> $func_body)

    fn = Base.eval(PadrinoEvalModule, func_expr)

    # Evaluate at all meshpoints
    values = Float64[]
    for j in 1:m
        val = Base.invokelatest(fn, z[j], z[j], h)
        push!(values, val)
    end

    return reduction_fn(values)
end

"""Replace _r_max/_r_min with max/min in expressions (for evaluation)."""
function _normalize_r_fns(expr::Expr)
    new_args = Any[]
    for (i, arg) in enumerate(expr.args)
        if i == 1 && expr.head == :call && arg in (:_r_max, :_r_min)
            push!(new_args, arg == :_r_max ? :max : :min)
        else
            push!(new_args, _normalize_r_fns(arg))
        end
    end
    return Expr(expr.head, new_args...)
end

function _normalize_r_fns(x)
    return x
end
