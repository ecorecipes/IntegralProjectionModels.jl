"""
Environmental variable sampler construction for stochastic parameter-resampled models.
"""

"""
    _build_env_sampler(pm::PadrinoModel, state_vars::Vector{String}) -> Function

Build an environment sampler function `t -> Dict{String, Float64}` that generates
random parameter values for stochastic parameter-resampled models.
"""
function _build_env_sampler(pm::IPM.PadrinoModel, state_vars::Vector{String})
    env_vars = pm.environmental_variables
    isempty(env_vars) && return t -> Dict{String, Float64}()

    # Parse each environmental variable specification
    samplers = Pair{String, Any}[]
    param_data = Dict{String, Any}()

    for ev in env_vars
        model_type = ev[:model_type]
        vr_name = ev[:vr_expr_name]
        env_range = get(ev, :env_range, "NULL")

        if model_type == "Parameter" || (env_range != "NULL" && !ismissing(env_range))
            # This is a data parameter used by the sampler
            if !ismissing(env_range) && env_range != "NULL"
                try
                    param_data[vr_name] = _parse_env_range(string(env_range))
                catch e
                    @warn "Failed to parse env_range for $vr_name" exception=e
                end
            end
        end

        if model_type in ("Substituted", "Evaluated")
            env_func = ev[:env_function]
            if !ismissing(env_func) && env_func != "NULL"
                # Build a sampler expression
                try
                    sampler_expr = translate_env_function(string(env_func), state_vars)
                    sampler_fn = _compile_env_sampler(sampler_expr, param_data)
                    push!(samplers, vr_name => sampler_fn)
                catch e
                    @warn "Failed to build env sampler for $vr_name" exception=e
                end
            end
        end
    end

    # Return a function that samples all environmental variables
    function env_sampler(t)
        result = Dict{String, Float64}()
        for (name, fn) in samplers
            try
                result[name] = fn(param_data)
            catch e
                @warn "Failed to sample env variable $name at t=$t" exception=e
            end
        end
        return result
    end

    return env_sampler
end

"""
    _compile_env_sampler(expr::Expr, param_data::Dict) -> Function

Compile an environment sampler expression into a function.
"""
function _compile_env_sampler(expr, param_data::Dict)
    # Build parameter assignments
    param_assigns = Expr[]
    for (k, v) in param_data
        if v isa Number
            push!(param_assigns, :($(Symbol(k)) = $v))
        end
    end

    func_body = Expr(:block, param_assigns..., :(return $expr))
    func_expr = :((data) -> $func_body)

    try
        return Base.eval(PadrinoEvalModule, func_expr)
    catch e
        @warn "Failed to compile env sampler" exception=e
        return (data) -> 0.0
    end
end

"""
    _parse_env_range(range_str::String) -> Any

Parse an R-style range/value specification for environmental data.
"""
function _parse_env_range(range_str::String)
    range_str = strip(range_str)

    if range_str == "NULL" || isempty(range_str)
        return nothing
    end

    # Try parsing as a number
    try
        return parse(Float64, range_str)
    catch
    end

    # Try parsing as a range (e.g., "1:10")
    if occursin(":", range_str)
        parts = split(range_str, ":")
        if length(parts) == 2
            try
                return collect(parse(Float64, strip(parts[1])):parse(Float64, strip(parts[2])))
            catch
            end
        end
    end

    # Try parsing as c(...) vector
    if startswith(range_str, "c(") && endswith(range_str, ")")
        inner = range_str[3:end-1]
        try
            return [parse(Float64, strip(v)) for v in split(inner, ",")]
        catch
        end
    end

    return range_str
end
