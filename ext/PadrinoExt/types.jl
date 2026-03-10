"""
Types for the PADRINO extension.
"""

# Table names in the PADRINO database
const PADRINO_TABLES = [
    "Metadata", "StateVariables", "ContinuousDomains",
    "IntegrationRules", "StateVectors", "IpmKernels", "VitalRateExpr",
    "ParameterValues", "EnvironmentalVariables", "ParSetIndices"
]

const PADRINO_BASE_URL = "https://raw.githubusercontent.com/padrinoDB/Padrino/main/padrino-database/clean/"

"""
    PadrinoKernelSpec

Parsed specification for a single PADRINO kernel.
"""
struct PadrinoKernelSpec
    kernel_id::String
    formula::String          # e.g., "K = P + F"
    model_family::String     # "CC", "CD", "DC", "DD"
    domain_start::String     # state variable name
    domain_end::String       # state variable name
    integration_rule::String # "midpoint"
end

"""
    PadrinoVitalRate

Parsed vital rate expression from PADRINO.
"""
struct PadrinoVitalRate
    name::String
    formula::String
    model_type::String  # "Evaluated" or "Substituted"
    kernel_ids::Vector{String}
end

"""
    PadrinoEnvVar

Environmental variable specification for stochastic models.
"""
struct PadrinoEnvVar
    env_variable::String
    vr_expr_name::String
    env_function::String
    env_range::String
    model_type::String  # "Substituted", "Evaluated", or "Parameter"
end

"""
    PadrinoParSetIndex

Parameter set index for kernel-resampled stochastic models.
"""
struct PadrinoParSetIndex
    kernel_id::String
    vr_expr_name::String
    range_str::String
    drop_levels::String
end

"""
    ModelTypeInfo

Determined model type classification for a PADRINO model.
"""
struct ModelTypeInfo
    sim_gen::String      # "simple" or "general"
    di_dd::String        # "di" or "dd"
    det_stoch::String    # "det" or "stoch"
    kern_param::Union{String, Nothing}  # "kern" or "param" (only for stoch)
    uses_age::Bool
    has_par_sets::Bool   # par_set_indices present (deterministic expansion)
end
