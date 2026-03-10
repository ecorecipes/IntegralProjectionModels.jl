"""
IPMProblem builders for all PADRINO model types.

Converts PadrinoModel specifications into IntegralProjectionModels.jl
IPMProblem objects for all supported model types.
"""

"""
    build_ipm(pm::PadrinoModel; tspan=(0, 100)) -> IPMProblem

Build an IPMProblem from a PadrinoModel, dispatching on model type.
"""
function build_ipm(pm::IPM.PadrinoModel; tspan=(0, 100))
    mtype = _determine_model_type(pm)

    if mtype.sim_gen == "simple" && mtype.di_dd == "di" && mtype.det_stoch == "det" && mtype.has_par_sets
        return _build_det_par_sets(pm, tspan, mtype)
    elseif mtype.sim_gen == "simple" && mtype.di_dd == "di" && mtype.det_stoch == "det"
        return _build_simple_di_det(pm, tspan)
    elseif mtype.sim_gen == "simple" && mtype.di_dd == "dd" && mtype.det_stoch == "det"
        return _build_simple_dd_det(pm, tspan)
    elseif mtype.sim_gen == "simple" && mtype.det_stoch == "stoch" && mtype.kern_param == "kern"
        return _build_stoch_kern(pm, tspan, mtype)
    elseif mtype.sim_gen == "simple" && mtype.det_stoch == "stoch" && mtype.kern_param == "param"
        return _build_stoch_param(pm, tspan, mtype)
    elseif mtype.sim_gen == "general" && mtype.det_stoch == "det"
        return _build_general_det(pm, tspan, mtype)
    elseif mtype.sim_gen == "general" && mtype.det_stoch == "stoch"
        return _build_general_stoch(pm, tspan, mtype)
    else
        error("Unsupported model type for ipm_id=$(pm.ipm_id): $mtype")
    end
end

"""
    _determine_model_type(pm::PadrinoModel) -> ModelTypeInfo

Determine the model classification from PADRINO metadata.
"""
function _determine_model_type(pm::IPM.PadrinoModel)
    md = pm.metadata

    # Simple vs general: multiple state variables or any discrete → general
    n_svs = length(pm.state_variables)
    has_discrete = any(sv -> get(sv, :discrete, false), pm.state_variables)
    sim_gen = (n_svs > 1 || has_discrete) ? "general" : "simple"

    # Density dependent
    has_dd = get(md, :has_dd, false)
    if ismissing(has_dd) || has_dd === nothing
        has_dd = false
    end
    di_dd = has_dd ? "dd" : "di"

    # Stochasticity: only environmental variables trigger stochastic classification.
    # par_set_indices alone means deterministic expansion (matching R's default behavior).
    has_env = !isempty(pm.environmental_variables)
    has_par_sets = any(psi -> !(get(psi, :vr_expr_name, "") in ("age", "max_age")),
                       pm.par_set_indices)

    if has_env
        det_stoch = "stoch"
        kern_param = "param"
    else
        det_stoch = "det"
        kern_param = nothing
    end

    uses_age = get(md, :has_age, false)
    if ismissing(uses_age) || uses_age === nothing
        uses_age = false
    end

    return ModelTypeInfo(sim_gen, di_dd, det_stoch, kern_param, uses_age, has_par_sets)
end

# --- Simple density-independent deterministic ---

function _build_simple_di_det(pm::IPM.PadrinoModel, tspan)
    state_vars = _get_continuous_state_vars(pm)
    domain = _build_domain(pm, state_vars[1])

    # Compile all vital rates
    vr_exprs = compile_vital_rates(pm.vital_rates, state_vars)

    # Build kernel matrices for each kernel and sum them
    K_total = zeros(Float64, domain.n_meshpoints, domain.n_meshpoints)

    for kern_spec in pm.kernels
        kern_id = string(kern_spec[:kernel_id])
        formula = string(kern_spec[:formula])
        family = string(get(kern_spec, :model_family, "CC"))

        # Skip iteration procedure / IPM-level aggregation kernels
        family in ("iteration_procedure", "IPM") && continue

        # Get the kernel formula RHS
        kern_rhs = _extract_formula_rhs(formula)
        kern_rhs === nothing && continue

        # Sanitize dots in kernel formula
        kern_rhs = _sanitize_formula(kern_rhs)

        # Filter vital rates for this kernel
        kern_vrs = _filter_vital_rates_for_kernel(pm, kern_id, state_vars)

        # Pre-compute mesh-level sum() aggregations (e.g., sum(r_s))
        precompute_mesh_sums!(kern_vrs, state_vars, pm.parameters, domain)

        # Pre-compute R-style max()/min() reductions (R's max() is a global
        # reduction, not element-wise like pmax(); this matters for VRs like
        # sigma_g = max(0.00001, sigma_g_temp) where sigma_g_temp is 0 for
        # some meshpoints — R returns a single scalar, not per-meshpoint values)
        precompute_r_reductions!(kern_vrs, state_vars, pm.parameters, domain)

        # Handle simple models: remove d_z from CC kernel formulas
        # (ipmr does this automatically for simple models)
        if family == "CC"
            kern_rhs = _remove_dz(kern_rhs, state_vars)
        end

        # Apply eviction correction
        eviction = _get_eviction_type(pm, kern_id)

        if eviction == "discrete_extrema" && family == "CC"
            # Proper discrete_extrema: apply correction to the distribution component
            K_mat = _build_kernel_with_discrete_extrema(kern_rhs, kern_vrs, pm, kern_id,
                                                         state_vars, domain)
        elseif eviction == "truncated_distributions"
            # Wrap distribution calls in Truncated(dist, L, U)
            apply_truncation_to_vrs!(kern_vrs, Float64(domain.lower), Float64(domain.upper))
            kern_fn, vec_fn = build_kernel_functions(kern_rhs, kern_vrs, state_vars, pm.parameters)
            K_mat = build_kernel_matrix(kern_fn, vec_fn, domain, family)
        else
            kern_fn, vec_fn = build_kernel_functions(kern_rhs, kern_vrs, state_vars, pm.parameters)
            K_mat = build_kernel_matrix(kern_fn, vec_fn, domain, family)
        end

        K_total .+= K_mat
    end

    # Initial population state
    n0 = IPM.uniform_population(domain)

    kernel = IPM.MatrixKernel(K_total)
    return IPM.IPMProblem(kernel, domain, n0, tspan)
end

# --- Simple density-dependent deterministic ---

function _build_simple_dd_det(pm::IPM.PadrinoModel, tspan)
    state_vars = _get_continuous_state_vars(pm)
    domain = _build_domain(pm, state_vars[1])

    # For DD models, we build a function that creates kernel from current pop state
    vr_exprs = compile_vital_rates(pm.vital_rates, state_vars)

    kernel_specs = [(ks, let r = _extract_formula_rhs(ks[:formula]); r === nothing ? nothing : _sanitize_formula(r) end,
                     get(ks, :model_family, "CC"))
                    for ks in pm.kernels
                    if get(ks, :model_family, "CC") ∉ ("iteration_procedure", "IPM")]

    function dd_kernel_fn(n_t, t, p)
        K_total = zeros(Float64, domain.n_meshpoints, domain.n_meshpoints)
        for (ks, kern_rhs, family) in kernel_specs
            kern_rhs === nothing && continue
            rhs = family == "CC" ? _remove_dz(kern_rhs, state_vars) : kern_rhs
            kern_vrs = _filter_vital_rates_for_kernel(pm, ks[:kernel_id], state_vars)
            kern_fn, vec_fn = build_kernel_functions(rhs, kern_vrs, state_vars, pm.parameters)
            K_mat = build_kernel_matrix(kern_fn, vec_fn, domain, family)
            K_total .+= K_mat
        end
        return IPM.MatrixKernel(K_total)
    end

    n0 = IPM.uniform_population(domain)
    return IPM.IPMProblem(IPM.DensityDependent(), dd_kernel_fn, domain, n0, tspan)
end

# --- Deterministic with parameter set expansion ---

"""
Build separate deterministic kernel matrices for each par_set combination.
Returns the first combination as a deterministic IPMProblem (matching R's
default behavior where lambda() reports per-combination eigenvalues).
"""
function _build_det_par_sets(pm::IPM.PadrinoModel, tspan, mtype::ModelTypeInfo)
    state_vars = _get_continuous_state_vars(pm)
    domain = _build_domain(pm, state_vars[1])

    # Parse par_set_indices to get the set of index values
    par_set_info, drop_levels = _parse_par_set_indices(pm)

    if isempty(par_set_info)
        # Fallback: build as simple deterministic
        return _build_simple_di_det(pm, tspan)
    end

    # Generate all combinations of index levels
    idx_vars = collect(keys(par_set_info))
    idx_values = [par_set_info[v] for v in idx_vars]
    combos = _cartesian_product(idx_values)

    # Filter out dropped level combinations
    if !isempty(drop_levels)
        combos = filter(combos) do combo
            combo_str = join([
                (val isa Number && isinteger(val)) ? string(Int(val)) : string(val)
                for val in combo
            ], "_")
            return combo_str ∉ drop_levels
        end
    end

    # Build kernel matrix for the first combination (matches R's default lambda reporting)
    combo = combos[1]

    # Create substitution map
    subs = Dict{String, String}()
    for (i, var_name) in enumerate(idx_vars)
        val = combo[i]
        subs[var_name] = (val isa Number && isinteger(val)) ? string(Int(val)) : string(val)
    end

    # Build the set of resolved kernel IDs for this combination
    # par_set_indices kernel_id field has semicolon-separated template IDs like "K_sp; P_sp; F_sp"
    valid_kernel_ids = _resolved_kernel_ids_for_combo(pm, subs)

    K_total = zeros(Float64, domain.n_meshpoints, domain.n_meshpoints)

    for kern_spec in pm.kernels
        family = string(get(kern_spec, :model_family, "CC"))
        family in ("iteration_procedure", "IPM") && continue

        kern_rhs = _extract_formula_rhs(string(kern_spec[:formula]))
        kern_rhs === nothing && continue

        kern_id = string(kern_spec[:kernel_id])

        # Only include kernels that belong to this par_set combination
        resolved_kid = _substitute_indices(kern_id, subs)
        if !isempty(valid_kernel_ids) && resolved_kid ∉ valid_kernel_ids
            continue
        end

        kern_rhs = _sanitize_formula(kern_rhs)
        kern_vrs = _filter_vital_rates_for_kernel(pm, kern_id, state_vars)

        # Substitute index variables in vital rate expressions
        resolved_vrs = _resolve_indexed_vrs(kern_vrs, subs)

        # Pre-compute mesh-level aggregations
        precompute_mesh_sums!(resolved_vrs, state_vars, pm.parameters, domain)
        precompute_r_reductions!(resolved_vrs, state_vars, pm.parameters, domain)

        if family == "CC"
            kern_rhs = _remove_dz(kern_rhs, state_vars)
        end

        # Substitute index vars in kernel formula too
        resolved_rhs = _substitute_indices(kern_rhs, subs)

        # Apply eviction correction
        eviction = _get_eviction_type(pm, kern_id)

        if eviction == "discrete_extrema" && family == "CC"
            K_mat = _build_kernel_with_discrete_extrema_resolved(
                resolved_rhs, resolved_vrs, pm, state_vars, domain)
        elseif eviction == "truncated_distributions"
            apply_truncation_to_vrs!(resolved_vrs, Float64(domain.lower), Float64(domain.upper))
            kern_fn, vec_fn = build_kernel_functions(resolved_rhs, resolved_vrs, state_vars, pm.parameters)
            K_mat = build_kernel_matrix(kern_fn, vec_fn, domain, family)
        else
            kern_fn, vec_fn = build_kernel_functions(resolved_rhs, resolved_vrs, state_vars, pm.parameters)
            K_mat = build_kernel_matrix(kern_fn, vec_fn, domain, family)
        end

        K_total .+= K_mat
    end

    n0 = IPM.uniform_population(domain)
    kernel = IPM.MatrixKernel(K_total)
    return IPM.IPMProblem(kernel, domain, n0, tspan)
end

# --- Stochastic kernel-resampled ---

function _build_stoch_kern(pm::IPM.PadrinoModel, tspan, mtype::ModelTypeInfo)
    state_vars = _get_continuous_state_vars(pm)
    domain = _build_domain(pm, state_vars[1])

    # Parse par_set_indices to get the set of index values
    par_set_info, drop_levels = _parse_par_set_indices(pm)

    if isempty(par_set_info)
        # Fallback: build as deterministic
        return _build_simple_di_det(pm, tspan)
    end

    # Generate all combinations of index levels
    idx_vars = collect(keys(par_set_info))
    idx_values = [par_set_info[v] for v in idx_vars]
    combos = _cartesian_product(idx_values)

    # Filter out dropped level combinations
    if !isempty(drop_levels)
        combos = filter(combos) do combo
            # Build the combo string (e.g., "2200_novel")
            combo_str = join([
                (val isa Number && isinteger(val)) ? string(Int(val)) : string(val)
                for val in combo
            ], "_")
            return combo_str ∉ drop_levels
        end
    end

    # Build one kernel matrix per combination
    kernel_matrices = Matrix{Float64}[]

    for combo in combos
        # Create substitution map: index_var_name => index_value
        subs = Dict{String, String}()
        for (i, var_name) in enumerate(idx_vars)
            val = combo[i]
            # Format integer-valued floats without decimal (1400.0 → "1400")
            subs[var_name] = (val isa Number && isinteger(val)) ? string(Int(val)) : string(val)
        end

        # Build the set of resolved kernel IDs for this combination
        valid_kernel_ids = _resolved_kernel_ids_for_combo(pm, subs)

        K_total = zeros(Float64, domain.n_meshpoints, domain.n_meshpoints)

        for kern_spec in pm.kernels
            family = string(get(kern_spec, :model_family, "CC"))
            family in ("iteration_procedure", "IPM") && continue

            kern_rhs = _extract_formula_rhs(string(kern_spec[:formula]))
            kern_rhs === nothing && continue

            kern_id = string(kern_spec[:kernel_id])

            # Only include kernels that belong to this par_set combination
            resolved_kid = _substitute_indices(kern_id, subs)
            if !isempty(valid_kernel_ids) && resolved_kid ∉ valid_kernel_ids
                continue
            end

            # Sanitize dots in kernel formula
            kern_rhs = _sanitize_formula(kern_rhs)

            kern_vrs = _filter_vital_rates_for_kernel(pm, kern_id, state_vars)

            # Substitute index variables in vital rate expressions
            resolved_vrs = _resolve_indexed_vrs(kern_vrs, subs)

            # Pre-compute mesh-level aggregations
            precompute_mesh_sums!(resolved_vrs, state_vars, pm.parameters, domain)
            precompute_r_reductions!(resolved_vrs, state_vars, pm.parameters, domain)

            if family == "CC"
                kern_rhs = _remove_dz(kern_rhs, state_vars)
            end

            # Substitute index vars in kernel formula too
            resolved_rhs = _substitute_indices(kern_rhs, subs)

            kern_fn, vec_fn = build_kernel_functions(resolved_rhs, resolved_vrs, state_vars, pm.parameters)
            K_mat = build_kernel_matrix(kern_fn, vec_fn, domain, family)
            K_total .+= K_mat
        end

        push!(kernel_matrices, K_total)
    end

    kernels = [IPM.MatrixKernel(K) for K in kernel_matrices]
    n0 = IPM.uniform_population(domain)
    return IPM.IPMProblem(IPM.StochasticKernelResampled(), kernels, domain, n0, tspan;
                          normalize=true)
end

# --- Stochastic parameter-resampled ---

function _build_stoch_param(pm::IPM.PadrinoModel, tspan, mtype::ModelTypeInfo)
    state_vars = _get_continuous_state_vars(pm)
    domain = _build_domain(pm, state_vars[1])

    # Build environment sampler
    env_sampler = _build_env_sampler(pm, state_vars)

    # Build kernel builder function that takes merged params
    vr_exprs_compiled = compile_vital_rates(pm.vital_rates, state_vars)

    base_params = copy(pm.parameters)

    # Kernel builder: params dict → AbstractIPMKernel
    function param_kernel_builder(sampled_params)
        # Merge base params with sampled env params
        merged = copy(base_params)
        for (k, v) in sampled_params
            merged[k] = v
        end

        K_total = zeros(Float64, domain.n_meshpoints, domain.n_meshpoints)
        for kern_spec in pm.kernels
            family = get(kern_spec, :model_family, "CC")
            family in ("iteration_procedure", "IPM") && continue

            kern_rhs = _extract_formula_rhs(kern_spec[:formula])
            kern_rhs === nothing && continue

            # Sanitize dots in kernel formula
            kern_rhs = _sanitize_formula(kern_rhs)

            kern_vrs = _filter_vital_rates_for_kernel(pm, kern_spec[:kernel_id], state_vars)
            precompute_mesh_sums!(kern_vrs, state_vars, merged, domain)
            precompute_r_reductions!(kern_vrs, state_vars, merged, domain)
            if family == "CC"
                kern_rhs = _remove_dz(kern_rhs, state_vars)
            end

            kern_fn, vec_fn = build_kernel_functions(kern_rhs, kern_vrs, state_vars, merged)
            K_mat = build_kernel_matrix(kern_fn, vec_fn, domain, family)
            K_total .+= K_mat
        end

        return IPM.MatrixKernel(K_total)
    end

    n0 = IPM.uniform_population(domain)
    return IPM.IPMProblem(IPM.StochasticParameterResampled(),
                          param_kernel_builder, domain, n0, tspan;
                          env_state=env_sampler, normalize=true)
end

# --- General (multi-state) models ---

function _build_general_det(pm::IPM.PadrinoModel, tspan, mtype::ModelTypeInfo)
    state_vars = _get_continuous_state_vars(pm)
    discrete_vars = _get_discrete_state_vars(pm)

    # Build domains for each state variable
    domains = Dict{Symbol, Any}()
    for sv in state_vars
        domains[Symbol(sv)] = _build_domain(pm, sv)
    end
    for dv in discrete_vars
        # Discrete domains have labels based on the model structure
        domains[Symbol(dv)] = IPM.DiscreteDomain([Symbol(dv)])
    end

    states_nt = NamedTuple(domains)

    # Build kernel matrices for each kernel block
    kernel_dict = Dict{Tuple{Symbol, Symbol}, IPM.AbstractIPMKernel}()

    for kern_spec in pm.kernels
        family = get(kern_spec, :model_family, "CC")
        family in ("iteration_procedure", "IPM") && continue

        kern_id = kern_spec[:kernel_id]
        kern_rhs = _extract_formula_rhs(kern_spec[:formula])
        kern_rhs === nothing && continue

        # Sanitize dots in kernel formula
        kern_rhs = _sanitize_formula(kern_rhs)

        domain_start = get(kern_spec, :domain_start, "")
        domain_end = get(kern_spec, :domain_end, "")

        from_state = Symbol(domain_start)
        to_state = Symbol(domain_end)

        all_state_names = vcat(state_vars, discrete_vars)
        kern_vrs = _filter_vital_rates_for_kernel(pm, kern_id, all_state_names)

        # For general models, don't remove d_z from CC kernels
        # (ipmr keeps d_z for general models)
        if family == "CC" && haskey(domains, from_state) && domains[from_state] isa IPM.ContinuousDomain
            dom = domains[from_state]
            kern_fn, vec_fn = build_kernel_functions(kern_rhs, kern_vrs, all_state_names, pm.parameters)
            K_mat = build_kernel_matrix(kern_fn, vec_fn, dom, family)
            kernel_dict[(from_state, to_state)] = IPM.MatrixKernel(K_mat)
        elseif family in ("CD", "DC", "DD")
            # For cross-domain kernels, determine dimensions and build matrix
            from_n = _state_size(domains, from_state)
            to_n = _state_size(domains, to_state)
            # Build with the appropriate domain (use target domain for row dim)
            if haskey(domains, to_state) && domains[to_state] isa IPM.ContinuousDomain
                dom = domains[to_state]
                kern_fn, vec_fn = build_kernel_functions(kern_rhs, kern_vrs, all_state_names, pm.parameters)
                K_mat = build_kernel_matrix(kern_fn, vec_fn, dom, family)
                kernel_dict[(from_state, to_state)] = IPM.MatrixKernel(K_mat)
            else
                # Discrete target — build a simple matrix
                K_mat = zeros(Float64, to_n, from_n)
                kernel_dict[(from_state, to_state)] = IPM.MatrixKernel(K_mat)
            end
        end
    end

    mega = IPM.MegaKernel(kernel_dict, states_nt)

    # Build initial population state
    total_n = sum(IPM.n_states(d) for d in values(domains))
    n0 = fill(1.0 / total_n, total_n)

    return IPM.IPMProblem(IPM.GeneralIPM(), mega, domains, n0, tspan)
end

function _build_general_stoch(pm::IPM.PadrinoModel, tspan, mtype::ModelTypeInfo)
    # For now, build as deterministic general model
    # Full stochastic general model support would require combining
    # multi-state with environmental or kernel resampling
    return _build_general_det(pm, tspan, mtype)
end

# --- Helper functions ---

function _get_continuous_state_vars(pm::IPM.PadrinoModel)
    svs = String[]
    for sv in pm.state_variables
        if !get(sv, :discrete, false)
            push!(svs, string(sv[:state_variable]))
        end
    end
    return svs
end

function _get_discrete_state_vars(pm::IPM.PadrinoModel)
    svs = String[]
    for sv in pm.state_variables
        if get(sv, :discrete, false)
            push!(svs, string(sv[:state_variable]))
        end
    end
    return svs
end

function _build_domain(pm::IPM.PadrinoModel, state_var::String)
    cd = nothing
    for d in pm.continuous_domains
        if string(d[:state_variable]) == state_var
            cd = d
            break
        end
    end
    cd === nothing && error("No continuous domain found for state variable: $state_var")

    # Find n_bins from state vectors
    n_bins = 100  # default
    for sv in pm.state_vectors
        expr_name = string(get(sv, :expression, ""))
        if occursin(state_var, expr_name)
            n_bins = round(Int, sv[:n_bins])
            break
        end
    end

    return IPM.ContinuousDomain(Float64(cd[:lower]), Float64(cd[:upper]), n_bins)
end

function _state_size(domains, state::Symbol)
    if haskey(domains, state)
        return IPM.n_states(domains[state])
    end
    return 1
end

function _extract_formula_rhs(formula::AbstractString)
    parts = split(formula, "="; limit=2)
    length(parts) == 2 || return nothing
    rhs = strip(parts[2])
    isempty(rhs) && return nothing
    return rhs
end

function _filter_vital_rates_for_kernel(pm::IPM.PadrinoModel, kernel_id::AbstractString,
                                         state_vars::Vector{String})
    compiled = Dict{String, Any}()
    for vr in pm.vital_rates
        # Check if this vital rate belongs to this kernel
        kr_ids = [string(k) for k in vr[:kernel_ids]]
        belongs = any(kid -> _kernel_id_matches(string(kernel_id), kid), kr_ids)
        belongs || continue

        formula = string(vr[:formula])
        model_type = string(vr[:model_type])
        name = _sanitize_name(vr[:name])

        _is_self_assignment(formula) && continue

        # Sanitize dots in the formula (R names like grow.int → grow_int)
        formula = _sanitize_formula(formula)

        parts = split(formula, "="; limit=2)
        length(parts) == 2 || continue
        rhs = strip(parts[2])

        if model_type == "Evaluated"
            rexpr = parse_rexpr(rhs)
            compiled[name] = translate_rexpr(rexpr, state_vars)
        elseif model_type == "Substituted"
            rexpr = parse_rexpr(rhs)
            compiled[name] = translate_rexpr(rexpr, state_vars)
        end
    end

    return compiled
end

function _kernel_id_matches(target::AbstractString, kid::AbstractString)
    # Kernel IDs in VitalRateExpr can be semicolon-separated
    parts = strip.(split(kid, ";"))
    return target in parts
end

function _remove_dz(formula::AbstractString, state_vars::Vector{String})
    result = formula
    for sv in state_vars
        # Remove " * d_size" patterns
        result = replace(result, r"\s*\*\s*d_" * sv => "")
    end
    return result
end

function _get_eviction_type(pm::IPM.PadrinoModel, kernel_id::AbstractString)
    md = pm.metadata
    eviction_used = get(md, :eviction_used, false)
    if ismissing(eviction_used) || eviction_used === nothing || !eviction_used
        return "none"
    end
    evict_type = get(md, :evict_type, "none")
    if ismissing(evict_type)
        return "none"
    end
    return string(evict_type)
end

function _parse_par_set_indices(pm::IPM.PadrinoModel)
    result = Dict{String, Vector{Any}}()
    drop_levels = Set{String}()

    for psi in pm.par_set_indices
        vr_name = string(psi[:vr_expr_name])
        vr_name in ("age", "max_age") && continue

        range_str = string(psi[:range_str])

        # Collect drop_levels
        dl = get(psi, :drop_levels, nothing)
        if dl !== nothing && !ismissing(dl) && string(dl) != ""
            for dlv in split(string(dl), ";")
                push!(drop_levels, strip(dlv))
            end
        end

        # Parse R range expressions like "2004:2012" or "c(1, 2, 3)"
        try
            range_str = replace(range_str, "'" => "")
            range_str = replace(range_str, "\"" => "")

            if occursin(":", range_str)
                # Handle range like "2004:2012"
                parts = split(range_str, ":")
                if length(parts) == 2
                    start_val = parse(Int, strip(parts[1]))
                    end_val = parse(Int, strip(parts[2]))
                    result[vr_name] = collect(start_val:end_val)
                end
            elseif startswith(range_str, "c(")
                # Handle c(1, 2, 3) or c('Open', 'Dense') style
                inner = range_str[3:end-1]
                vals = Any[]
                for v in split(inner, ",")
                    sv = strip(v)
                    try
                        push!(vals, parse(Float64, sv))
                    catch
                        push!(vals, sv)  # Keep as string
                    end
                end
                result[vr_name] = vals
            else
                # Try parsing as a single value
                try
                    result[vr_name] = [parse(Float64, strip(range_str))]
                catch
                    result[vr_name] = [strip(range_str)]
                end
            end
        catch e
            @warn "Failed to parse par_set_index range" range_str exception=e
        end
    end

    return result, drop_levels
end

"""
    _resolved_kernel_ids_for_combo(pm, subs) -> Set{String}

Build the set of resolved kernel IDs that belong to a given par_set combination.
Extracts template kernel IDs from par_set_indices (e.g., "K_sp; P_sp; F_sp"),
resolves them with the substitution map (sp→1 gives "K_1; P_1; F_1"),
and returns the set of valid IDs.
"""
function _resolved_kernel_ids_for_combo(pm::IPM.PadrinoModel, subs::Dict{String, String})
    valid = Set{String}()
    for psi in pm.par_set_indices
        kid_str = get(psi, :kernel_id, "")
        (ismissing(kid_str) || kid_str == "") && continue
        for tmpl in split(string(kid_str), ";")
            resolved = _substitute_indices(strip(tmpl), subs)
            push!(valid, resolved)
        end
    end
    return valid
end

"""
    _build_kernel_with_discrete_extrema_resolved(...)

Like _build_kernel_with_discrete_extrema but works with already-resolved VRs
(index variables already substituted). Finds the Substituted distribution
component from the resolved vital rates dict rather than from pm.vital_rates.
"""
function _build_kernel_with_discrete_extrema_resolved(
        kern_rhs::AbstractString,
        resolved_vrs::Dict{String, Any},
        pm::IPM.PadrinoModel,
        state_vars::Vector{String},
        domain::IPM.ContinuousDomain)
    z = IPM.meshpoints(domain)
    h = IPM.step_size(domain)
    m = length(z)

    # Find the Substituted vital rate (the distribution) among resolved VRs.
    # A distribution VR is one whose compiled expression is a pdf() call.
    dist_vr_name = nothing
    for (name, expr) in resolved_vrs
        if _is_pdf_expr(expr)
            dist_vr_name = name
            break
        end
    end

    if dist_vr_name === nothing
        # No distribution found — fall back to standard build
        kern_fn = build_kernel_function(kern_rhs, resolved_vrs, state_vars, pm.parameters)
        return build_kernel_matrix(kern_fn, domain, "CC")
    end

    # Build scalar function: everything except the distribution
    scalar_rhs = replace(kern_rhs, Regex("\\b" * dist_vr_name * "\\b") => "1.0")
    scalar_vrs = Dict{String, Any}(k => v for (k, v) in resolved_vrs if k != dist_vr_name)
    scalar_fn = build_kernel_function(scalar_rhs, scalar_vrs, state_vars, pm.parameters)

    # Build distribution function
    dist_vrs = Dict{String, Any}()
    dist_vrs[dist_vr_name] = resolved_vrs[dist_vr_name]
    for (name, expr) in resolved_vrs
        if name != dist_vr_name && _expr_references(resolved_vrs[dist_vr_name], Symbol(name))
            dist_vrs[name] = expr
        end
    end
    dist_fn, dist_vec_fn = build_kernel_functions(dist_vr_name, dist_vrs, state_vars, pm.parameters)

    # Build distribution matrix using vectorized evaluation when possible
    D = _build_matrix_vectorized_or_scalar(dist_fn, dist_vec_fn, z, h, m)
    D .*= h

    # Correction: adjust first/last row so columns sum to 1
    mid = m ÷ 2
    @inbounds for j in 1:m
        col_sum = sum(view(D, :, j))
        if col_sum > 0
            deficit = 1.0 - col_sum
            if j <= mid
                D[1, j] += deficit
            else
                D[m, j] += deficit
            end
        end
    end

    # Build kernel: K[i,j] = scalar(z_j) * D_corrected[i,j]
    K = similar(D)
    @inbounds for j in 1:m
        s_val = Base.invokelatest(scalar_fn, z[j], z[j], h)::Float64
        for i in 1:m
            K[i, j] = s_val * D[i, j]
        end
    end

    return K
end

"""
    _build_matrix_vectorized_or_scalar(scalar_fn, vec_fn, z, h, m) -> Matrix{Float64}

Try interpreted vectorized evaluation first, then compiled vectorized, then scalar loop.
"""
function _build_matrix_vectorized_or_scalar(scalar_fn, vec_fn, z::Vector{Float64}, h::Float64, m::Int)
    # Path 1: Interpreted vectorized evaluation
    if scalar_fn isa Expr
        try
            z_new = reshape(z, :, 1)
            z_src = reshape(z, 1, :)
            env = Dict{Symbol, Any}(:z_new => z_new, :z => z_src, :h => h)
            K = _eval_vectorized(scalar_fn, env)
            K = convert(Matrix{Float64}, K)
            if size(K) == (m, m)
                return K
            end
        catch
        end
    end
    # Path 2: Compiled vectorized evaluation
    if vec_fn !== nothing
        try
            z_new = reshape(z, :, 1)
            z_src = reshape(z, 1, :)
            K = Base.invokelatest(vec_fn, z_new, z_src, h)
            K = convert(Matrix{Float64}, K)
            if size(K) == (m, m)
                return K
            end
        catch
        end
    end
    # Path 3: Compiled scalar fallback
    compiled_fn = _ensure_compiled(scalar_fn)
    K = zeros(Float64, m, m)
    @inbounds for j in 1:m
        zj = z[j]
        for i in 1:m
            K[i, j] = Base.invokelatest(compiled_fn, z[i], zj, h)::Float64
        end
    end
    return K
end

"""Check if an expression is a pdf() call."""
function _is_pdf_expr(expr::Expr)
    return expr.head == :call && length(expr.args) >= 1 && expr.args[1] == :pdf
end

_is_pdf_expr(::Any) = false

"""
    _cartesian_product(arrays) -> Vector{Vector}

Generate all combinations of elements from multiple arrays.
"""
function _cartesian_product(arrays::Vector)
    if isempty(arrays)
        return [Any[]]
    end
    if length(arrays) == 1
        return [[x] for x in arrays[1]]
    end

    rest = _cartesian_product(arrays[2:end])
    result = Vector{Any}[]
    for x in arrays[1]
        for r in rest
            push!(result, vcat([x], r))
        end
    end
    return result
end

"""
    _resolve_indexed_params(params, subs) -> Dict{String, Float64}

Resolve parameter names by substituting index variables.
E.g., with subs = {"yr" => "2014", "pl" => "1"}:
  - "s_yearan_yr" is not a parameter, but the expression references it
  - We add "s_yearan_yr" => params["s_yearan_2014"]
  - "s_ptran_pl" => params["s_ptran_1"]
"""
function _resolve_indexed_params(params::Dict{String, Float64}, subs::Dict{String, String})
    resolved = copy(params)

    # For each index variable, find parameters that match the pattern
    # param_name_indexval and create entries for param_name_indexvar
    for (idx_var, idx_val) in subs
        suffix_var = "_$(idx_var)"
        suffix_val = "_$(idx_val)"

        for (name, val) in params
            if endswith(name, suffix_val)
                base = name[1:end-length(suffix_val)]
                var_name = base * suffix_var
                # Map the variable-suffixed name to the resolved value
                resolved[var_name] = val
            end
        end
    end

    return resolved
end

"""
    _resolve_indexed_vrs(vrs::Dict{String, Any}, subs::Dict{String, String}) -> Dict{String, Any}

Resolve index variable references in compiled vital rate expressions.
Substitutes Symbol names like :s_yearan_yr → :s_yearan_2014.
"""
function _resolve_indexed_vrs(vrs::Dict{String, Any}, subs::Dict{String, String})
    result = Dict{String, Any}()
    for (name, expr) in vrs
        # Substitute index variables in the expression
        new_expr = _substitute_expr_indices(expr, subs)
        # Also resolve the vital rate name itself
        new_name = _substitute_indices(name, subs)
        result[new_name] = new_expr
    end
    return result
end

"""
    _substitute_indices(s::AbstractString, subs::Dict{String, String}) -> String

Replace index variable names with their values in a string.
E.g., "s_yearan_yr" with subs={"yr"=>"2014"} → "s_yearan_2014"
"""
function _substitute_indices(s::AbstractString, subs::Dict{String, String})
    result = string(s)
    for (idx_var, idx_val) in subs
        result = replace(result, "_$(idx_var)" => "_$(idx_val)")
    end
    return result
end

"""
    _substitute_expr_indices(expr, subs) -> Expr

Recursively substitute index variables in a Julia expression.
"""
function _substitute_expr_indices(expr::Expr, subs::Dict{String, String})
    new_args = Any[]
    for arg in expr.args
        push!(new_args, _substitute_expr_indices(arg, subs))
    end
    return Expr(expr.head, new_args...)
end

function _substitute_expr_indices(sym::Symbol, subs::Dict{String, String})
    s = String(sym)
    new_s = _substitute_indices(s, subs)
    return new_s == s ? sym : Symbol(new_s)
end

function _substitute_expr_indices(x, subs::Dict{String, String})
    return x  # numbers, strings, etc.
end

"""
    _build_kernel_with_discrete_extrema(kern_rhs, kern_vrs, pm, kern_id, state_vars, domain)

Build a kernel matrix with proper discrete_extrema eviction correction.

The approach:
1. Identify the Substituted vital rate (the distribution component)
2. Build the distribution matrix and apply discrete_extrema correction
   (redistribute lost mass to boundary rows so columns sum to 1)
3. Build the scalar part from remaining vital rates + parameters
4. Combine: K[i,j] = scalar(z_j) * corrected_dist[i,j] * h
"""
function _build_kernel_with_discrete_extrema(kern_rhs::AbstractString,
                                              kern_vrs::Dict{String, Any},
                                              pm::IPM.PadrinoModel,
                                              kern_id::AbstractString,
                                              state_vars::Vector{String},
                                              domain::IPM.ContinuousDomain)
    z = IPM.meshpoints(domain)
    h = IPM.step_size(domain)
    m = length(z)

    # Find the Substituted vital rate (the distribution) for this kernel
    dist_vr_name = nothing
    for vr in pm.vital_rates
        kr_ids = [string(k) for k in vr[:kernel_ids]]
        belongs = any(kid -> _kernel_id_matches(string(kern_id), kid), kr_ids)
        belongs || continue
        if string(vr[:model_type]) == "Substituted"
            dist_vr_name = _sanitize_name(vr[:name])
            break
        end
    end

    if dist_vr_name === nothing
        # No distribution found — fall back to standard build
        kern_fn = build_kernel_function(kern_rhs, kern_vrs, state_vars, pm.parameters)
        return build_kernel_matrix(kern_fn, domain, "CC")
    end

    # Build scalar function: everything in the kernel except the distribution
    # Replace the distribution name with 1.0 in the kernel formula
    scalar_rhs = replace(kern_rhs, Regex("\\b" * dist_vr_name * "\\b") => "1.0")
    scalar_vrs = Dict{String, Any}(k => v for (k, v) in kern_vrs if k != dist_vr_name)
    # Also remove vital rates that are dependencies of the distribution only
    # (but keep vital rates that the scalar part also needs)

    # Build the scalar function: (z_new, z, h) -> scalar_value
    scalar_fn = build_kernel_function(scalar_rhs, scalar_vrs, state_vars, pm.parameters)

    # Build the distribution function
    dist_vrs = Dict{String, Any}()
    dist_vrs[dist_vr_name] = kern_vrs[dist_vr_name]
    # Include dependencies of the distribution vital rate
    for (name, expr) in kern_vrs
        if name != dist_vr_name && _expr_references(kern_vrs[dist_vr_name], Symbol(name))
            dist_vrs[name] = expr
        end
    end
    dist_fn, dist_vec_fn = build_kernel_functions(dist_vr_name, dist_vrs, state_vars, pm.parameters)

    # Build distribution matrix using vectorized evaluation when possible
    D = _build_matrix_vectorized_or_scalar(dist_fn, dist_vec_fn, z, h, m)
    D .*= h

    # Apply discrete_extrema correction: adjust first/last row so columns sum to 1
    mid = m ÷ 2
    @inbounds for j in 1:m
        col_sum = sum(view(D, :, j))
        if col_sum > 0
            deficit = 1.0 - col_sum
            if j <= mid
                D[1, j] += deficit
            else
                D[m, j] += deficit
            end
        end
    end

    # Build the kernel matrix: K[i,j] = scalar(z_j) * D_corrected[i,j]
    K = similar(D)
    @inbounds for j in 1:m
        s_val = Base.invokelatest(scalar_fn, z[j], z[j], h)::Float64
        for i in 1:m
            K[i, j] = s_val * D[i, j]
        end
    end

    return K
end

"""
    _expr_references(expr, sym::Symbol) -> Bool

Check if a Julia expression references a given symbol.
"""
function _expr_references(expr::Expr, sym::Symbol)
    for arg in expr.args
        if _expr_references(arg, sym)
            return true
        end
    end
    return false
end

function _expr_references(s::Symbol, sym::Symbol)
    return s == sym
end

function _expr_references(::Any, ::Symbol)
    return false
end


