"""
PADRINO database integration stubs.

The PADRINO database contains ~280 published IPMs stored as tabular data.
Full functionality requires `using CSV, DataFrames, Downloads` to load the
`IntegralProjectionModelsPadrinoExt` extension.
"""

"""
    PadrinoModel

Parsed representation of a single PADRINO model, including vital rate expressions,
kernel formulas, parameter values, domain specifications, and metadata.
Created by `pdb_make_proto_ipm`.
"""
struct PadrinoModel
    ipm_id::String
    metadata::Dict{Symbol, Any}
    state_variables::Vector{Dict{Symbol, Any}}
    continuous_domains::Vector{Dict{Symbol, Any}}
    integration_rules::Vector{Dict{Symbol, Any}}
    state_vectors::Vector{Dict{Symbol, Any}}
    kernels::Vector{Dict{Symbol, Any}}
    vital_rates::Vector{Dict{Symbol, Any}}
    parameters::Dict{String, Float64}
    environmental_variables::Vector{Dict{Symbol, Any}}
    par_set_indices::Vector{Dict{Symbol, Any}}
end

"""
    PadrinoDB

Container for the PADRINO database tables. Holds 10 DataFrames corresponding
to the PADRINO tables: Metadata, StateVariables, ContinuousDomains,
IntegrationRules, StateVectors, IpmKernels, VitalRateExpr, ParameterValues,
EnvironmentalVariables, and ParSetIndices.

Requires `using CSV, DataFrames, Downloads` to access full functionality.
"""
struct PadrinoDB
    tables::Dict{String, Any}
end

"""
    pdb_download(; save=false, destination=nothing) -> PadrinoDB

Download the PADRINO database from GitHub.
Requires the PADRINO extension: `using CSV, DataFrames, Downloads`.
"""
function pdb_download end

"""
    pdb_load(path::String) -> PadrinoDB

Load PADRINO database tables from a directory of tab-delimited text files.
Requires the PADRINO extension: `using CSV, DataFrames, Downloads`.
"""
function pdb_load end

"""
    pdb_save(pdb::PadrinoDB, destination::String)

Save a PadrinoDB to a directory of tab-delimited text files.
Requires the PADRINO extension: `using CSV, DataFrames, Downloads`.
"""
function pdb_save end

"""
    pdb_subset(pdb::PadrinoDB, ipm_ids::Vector{String}) -> PadrinoDB

Subset a PadrinoDB to specific IPM IDs.
Requires the PADRINO extension: `using CSV, DataFrames, Downloads`.
"""
function pdb_subset end

"""
    pdb_species(pdb::PadrinoDB; ipm_id=nothing) -> Vector{String}

Return species names in the database, optionally filtered by IPM ID.
Requires the PADRINO extension: `using CSV, DataFrames, Downloads`.
"""
function pdb_species end

"""
    pdb_citations(pdb::PadrinoDB; ipm_id=nothing)

Return citation information, optionally filtered by IPM ID.
Requires the PADRINO extension: `using CSV, DataFrames, Downloads`.
"""
function pdb_citations end

"""
    pdb_metadata(pdb::PadrinoDB, column::Symbol; ipm_id=nothing)

Return metadata column values, optionally filtered by IPM ID.
Requires the PADRINO extension: `using CSV, DataFrames, Downloads`.
"""
function pdb_metadata end

"""
    pdb_test_targets(pdb::PadrinoDB; ipm_id=nothing)

Return test target values (expected lambda, etc.) for validation.
Requires the PADRINO extension: `using CSV, DataFrames, Downloads`.
"""
function pdb_test_targets end

"""
    pdb_make_proto_ipm(pdb::PadrinoDB; ipm_id=nothing, det_stoch="det", kern_param="kern") -> Dict{String, PadrinoModel}

Parse PADRINO database entries into PadrinoModel objects suitable for building IPMs.
Requires the PADRINO extension: `using CSV, DataFrames, Downloads`.
"""
function pdb_make_proto_ipm end

"""
    pdb_make_ipm(protos::Dict{String, PadrinoModel}; tspan=(0, 100)) -> Dict{String, IPMProblem}
    pdb_make_ipm(pdb::PadrinoDB; ipm_id=nothing, tspan=(0, 100), det_stoch="det", kern_param="kern") -> Dict{String, IPMProblem}

Build IPMProblem objects from PadrinoModel prototypes or directly from a PadrinoDB.
Requires the PADRINO extension: `using CSV, DataFrames, Downloads`.
"""
function pdb_make_ipm end

"""
    pdb_validate(ipms::Dict, pdb::PadrinoDB; rtol=0.01)

Validate built IPMs against PADRINO test targets.
Requires the PADRINO extension: `using CSV, DataFrames, Downloads`.
"""
function pdb_validate end
