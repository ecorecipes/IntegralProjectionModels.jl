"""
PADRINO database integration extension for IntegralProjectionModels.jl.

Activated by `using IntegralProjectionModels, CSV, DataFrames, Downloads`.
Provides functions to download, parse, and build IPMs from the PADRINO database.
"""
module IntegralProjectionModelsPadrinoExt

import IntegralProjectionModels as IPM
import CSV
import DataFrames
import Downloads
using Distributions

# Evaluation module for compiled kernel functions
# This module provides access to Distributions functions in eval'd code
module PadrinoEvalModule
    using Distributions
    using Statistics: mean, std
    # R-style reduction functions: R's max()/min() are global reductions
    # (unlike pmax/pmin which are element-wise). These fallbacks handle
    # cases not pre-computed by precompute_r_reductions!.
    _r_max(args...) = max(args...)
    _r_min(args...) = min(args...)
end

# Include extension components
include("PadrinoExt/types.jl")
include("PadrinoExt/parser.jl")
include("PadrinoExt/translator.jl")
include("PadrinoExt/compiler.jl")
include("PadrinoExt/eviction.jl")
include("PadrinoExt/env_state.jl")
include("PadrinoExt/builder.jl")
include("PadrinoExt/io.jl")
include("PadrinoExt/query.jl")
include("PadrinoExt/validation.jl")

# --- Override stub methods ---

function IPM.pdb_download(; save::Bool=false, destination=nothing)
    return _pdb_download(; save=save, destination=destination)
end

function IPM.pdb_load(path::String)
    return _pdb_load(path)
end

function IPM.pdb_save(pdb::IPM.PadrinoDB, destination::String)
    return _pdb_save(pdb, destination)
end

function IPM.pdb_subset(pdb::IPM.PadrinoDB, ipm_ids::Vector{String})
    return _pdb_subset(pdb, ipm_ids)
end

function IPM.pdb_species(pdb::IPM.PadrinoDB; ipm_id=nothing)
    return _pdb_species(pdb; ipm_id=ipm_id)
end

function IPM.pdb_citations(pdb::IPM.PadrinoDB; ipm_id=nothing)
    return _pdb_citations(pdb; ipm_id=ipm_id)
end

function IPM.pdb_metadata(pdb::IPM.PadrinoDB, column::Symbol; ipm_id=nothing)
    return _pdb_metadata(pdb, column; ipm_id=ipm_id)
end

function IPM.pdb_test_targets(pdb::IPM.PadrinoDB; ipm_id=nothing)
    return _pdb_test_targets(pdb; ipm_id=ipm_id)
end

function IPM.pdb_make_proto_ipm(pdb::IPM.PadrinoDB; ipm_id=nothing,
                                 det_stoch="det", kern_param="kern")
    return _pdb_make_proto_ipm(pdb; ipm_id=ipm_id, det_stoch=det_stoch,
                                kern_param=kern_param)
end

function IPM.pdb_make_ipm(protos::Dict{String, IPM.PadrinoModel}; tspan=(0, 100))
    return _pdb_make_ipm_from_protos(protos; tspan=tspan)
end

function IPM.pdb_make_ipm(pdb::IPM.PadrinoDB; ipm_id=nothing, tspan=(0, 100),
                           det_stoch="det", kern_param="kern")
    protos = _pdb_make_proto_ipm(pdb; ipm_id=ipm_id, det_stoch=det_stoch,
                                  kern_param=kern_param)
    return _pdb_make_ipm_from_protos(protos; tspan=tspan)
end

function IPM.pdb_validate(ipms::Dict, pdb::IPM.PadrinoDB; rtol=0.01)
    return _pdb_validate(ipms, pdb; rtol=rtol)
end

# --- Proto IPM construction ---

"""
    _pdb_make_proto_ipm(pdb; ipm_id, det_stoch, kern_param) -> Dict{String, PadrinoModel}

Parse PADRINO database entries into PadrinoModel objects.
"""
function _pdb_make_proto_ipm(pdb::IPM.PadrinoDB; ipm_id=nothing,
                              det_stoch="det", kern_param="kern")
    # Get IPM IDs to process
    ids = if ipm_id === nothing
        _get_ipm_ids(pdb)
    elseif ipm_id isa AbstractString
        [string(ipm_id)]
    else
        [string(x) for x in ipm_id]
    end

    result = Dict{String, IPM.PadrinoModel}()

    for id in ids
        try
            pm = _parse_single_model(pdb, string(id))
            result[string(id)] = pm
        catch e
            @warn "Failed to parse model $id" exception=e
        end
    end

    return result
end

"""
    _parse_single_model(pdb, ipm_id) -> PadrinoModel

Parse a single PADRINO model entry into a PadrinoModel.
"""
function _parse_single_model(pdb::IPM.PadrinoDB, ipm_id::AbstractString)
    # Extract relevant rows from each table
    metadata = _extract_metadata(pdb, ipm_id)
    state_variables = _extract_state_variables(pdb, ipm_id)
    continuous_domains = _extract_continuous_domains(pdb, ipm_id)
    integration_rules = _extract_integration_rules(pdb, ipm_id)
    state_vectors = _extract_state_vectors(pdb, ipm_id)
    kernels = _extract_kernels(pdb, ipm_id)
    vital_rates = _extract_vital_rates(pdb, ipm_id)
    parameters = _extract_parameters(pdb, ipm_id)
    env_vars = _extract_env_variables(pdb, ipm_id)
    par_set_indices = _extract_par_set_indices(pdb, ipm_id)

    return IPM.PadrinoModel(
        ipm_id,
        metadata,
        state_variables,
        continuous_domains,
        integration_rules,
        state_vectors,
        kernels,
        vital_rates,
        parameters,
        env_vars,
        par_set_indices
    )
end

# --- Table extraction helpers ---

function _extract_metadata(pdb::IPM.PadrinoDB, ipm_id::String)
    md = _get_table(pdb, "Metadata")
    DataFrames.nrow(md) == 0 && return Dict{Symbol, Any}()

    rows = filter(r -> string(r.ipm_id) == ipm_id, md)
    DataFrames.nrow(rows) == 0 && return Dict{Symbol, Any}()

    row = rows[1, :]
    result = Dict{Symbol, Any}()
    for col in DataFrames.propertynames(md)
        result[col] = row[col]
    end
    return result
end

function _extract_state_variables(pdb::IPM.PadrinoDB, ipm_id::String)
    sv = _get_table(pdb, "StateVariables")
    DataFrames.nrow(sv) == 0 && return Dict{Symbol, Any}[]

    rows = filter(r -> string(r.ipm_id) == ipm_id, sv)
    result = Dict{Symbol, Any}[]
    for row in DataFrames.eachrow(rows)
        d = Dict{Symbol, Any}()
        d[:state_variable] = string(row.state_variable)
        disc = get(row, :discrete, false)
        d[:discrete] = !ismissing(disc) && (disc == true || disc == 1 || disc == "TRUE")
        push!(result, d)
    end
    return result
end

function _extract_continuous_domains(pdb::IPM.PadrinoDB, ipm_id::String)
    cd = _get_table(pdb, "ContinuousDomains")
    DataFrames.nrow(cd) == 0 && return Dict{Symbol, Any}[]

    rows = filter(r -> string(r.ipm_id) == ipm_id, cd)
    result = Dict{Symbol, Any}[]
    for row in DataFrames.eachrow(rows)
        d = Dict{Symbol, Any}()
        d[:state_variable] = string(row.state_variable)
        d[:lower] = Float64(row.lower)
        d[:upper] = Float64(row.upper)
        push!(result, d)
    end
    return result
end

function _extract_integration_rules(pdb::IPM.PadrinoDB, ipm_id::String)
    ir = _get_table(pdb, "IntegrationRules")
    DataFrames.nrow(ir) == 0 && return Dict{Symbol, Any}[]

    rows = filter(r -> string(r.ipm_id) == ipm_id, ir)
    result = Dict{Symbol, Any}[]
    for row in DataFrames.eachrow(rows)
        d = Dict{Symbol, Any}()
        d[:kernel_id] = string(row.kernel_id)
        d[:integration_rule] = string(row.integration_rule)
        push!(result, d)
    end
    return result
end

function _extract_state_vectors(pdb::IPM.PadrinoDB, ipm_id::String)
    sv = _get_table(pdb, "StateVectors")
    DataFrames.nrow(sv) == 0 && return Dict{Symbol, Any}[]

    rows = filter(r -> string(r.ipm_id) == ipm_id, sv)
    result = Dict{Symbol, Any}[]
    for row in DataFrames.eachrow(rows)
        d = Dict{Symbol, Any}()
        d[:expression] = string(row.expression)
        d[:n_bins] = ismissing(row.n_bins) ? 100 : Int(row.n_bins)
        push!(result, d)
    end
    return result
end

function _extract_kernels(pdb::IPM.PadrinoDB, ipm_id::String)
    ik = _get_table(pdb, "IpmKernels")
    DataFrames.nrow(ik) == 0 && return Dict{Symbol, Any}[]

    rows = filter(r -> string(r.ipm_id) == ipm_id, ik)
    result = Dict{Symbol, Any}[]
    for row in DataFrames.eachrow(rows)
        d = Dict{Symbol, Any}()
        d[:kernel_id] = string(row.kernel_id)
        d[:model_family] = ismissing(row.model_family) ? "CC" : string(row.model_family)
        d[:formula] = string(row.formula)
        d[:domain_start] = ismissing(row.domain_start) ? "" : string(row.domain_start)
        d[:domain_end] = ismissing(row.domain_end) ? "" : string(row.domain_end)
        push!(result, d)
    end
    return result
end

function _extract_vital_rates(pdb::IPM.PadrinoDB, ipm_id::String)
    vr = _get_table(pdb, "VitalRateExpr")
    DataFrames.nrow(vr) == 0 && return Dict{Symbol, Any}[]

    rows = filter(r -> string(r.ipm_id) == ipm_id, vr)
    result = Dict{Symbol, Any}[]
    for row in DataFrames.eachrow(rows)
        d = Dict{Symbol, Any}()
        formula = string(row.formula)

        # Extract name from LHS of formula
        parts = split(formula, "="; limit=2)
        name = length(parts) >= 1 ? strip(parts[1]) : ""

        d[:name] = name
        d[:formula] = formula
        d[:model_type] = ismissing(row.model_type) ? "Evaluated" : string(row.model_type)
        d[:kernel_ids] = ismissing(row.kernel_id) ? String[] : [strip(string(row.kernel_id))]
        push!(result, d)
    end
    return result
end

function _extract_parameters(pdb::IPM.PadrinoDB, ipm_id::String)
    pv = _get_table(pdb, "ParameterValues")
    DataFrames.nrow(pv) == 0 && return Dict{String, Float64}()

    rows = filter(r -> string(r.ipm_id) == ipm_id, pv)
    result = Dict{String, Float64}()
    for row in DataFrames.eachrow(rows)
        name = string(row.parameter_name)
        val = row.parameter_value
        if !ismissing(val)
            try
                result[name] = Float64(val)
            catch
                @warn "Non-numeric parameter value for $name: $val"
            end
        end
    end
    return result
end

function _extract_env_variables(pdb::IPM.PadrinoDB, ipm_id::String)
    ev = _get_table(pdb, "EnvironmentalVariables")
    DataFrames.nrow(ev) == 0 && return Dict{Symbol, Any}[]

    rows = filter(r -> string(r.ipm_id) == ipm_id, ev)
    result = Dict{Symbol, Any}[]
    for row in DataFrames.eachrow(rows)
        d = Dict{Symbol, Any}()
        d[:env_variable] = ismissing(row.env_variable) ? "" : string(row.env_variable)
        d[:vr_expr_name] = ismissing(row.vr_expr_name) ? "" : string(row.vr_expr_name)
        d[:env_function] = ismissing(row.env_function) ? "" : string(row.env_function)
        d[:env_range] = ismissing(row.env_range) ? "NULL" : string(row.env_range)
        d[:model_type] = ismissing(row.model_type) ? "Evaluated" : string(row.model_type)
        push!(result, d)
    end
    return result
end

function _extract_par_set_indices(pdb::IPM.PadrinoDB, ipm_id::String)
    psi = _get_table(pdb, "ParSetIndices")
    DataFrames.nrow(psi) == 0 && return Dict{Symbol, Any}[]

    rows = filter(r -> string(r.ipm_id) == ipm_id, psi)
    result = Dict{Symbol, Any}[]
    for row in DataFrames.eachrow(rows)
        d = Dict{Symbol, Any}()
        d[:kernel_id] = ismissing(row.kernel_id) ? "" : string(row.kernel_id)
        d[:vr_expr_name] = ismissing(row.vr_expr_name) ? "" : string(row.vr_expr_name)
        d[:range_str] = ismissing(row.range) ? "" : string(row.range)
        d[:drop_levels] = hasproperty(row, :drop_levels) ?
            (ismissing(row.drop_levels) ? "" : string(row.drop_levels)) : ""
        push!(result, d)
    end
    return result
end

# --- IPM building from protos ---

function _pdb_make_ipm_from_protos(protos::Dict{String, IPM.PadrinoModel}; tspan=(0, 100))
    result = Dict{String, IPM.IPMProblem}()

    for (id, pm) in protos
        try
            ipm = build_ipm(pm; tspan=tspan)
            result[id] = ipm
        catch e
            @warn "Failed to build IPM for $id" exception=(e, catch_backtrace())
        end
    end

    return result
end

end # module
