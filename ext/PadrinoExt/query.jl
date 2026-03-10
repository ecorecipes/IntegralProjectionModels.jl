"""
PADRINO database querying functions.
"""

"""
    _pdb_subset(pdb::PadrinoDB, ipm_ids::Vector{String}) -> PadrinoDB

Subset a PadrinoDB to specific IPM IDs.
"""
function _pdb_subset(pdb::IPM.PadrinoDB, ipm_ids::Vector{String})
    tables = Dict{String, Any}()

    for (tab_name, df) in pdb.tables
        if DataFrames.ncol(df) > 0 && :ipm_id in DataFrames.propertynames(df)
            tables[tab_name] = DataFrames.subset(df, :ipm_id => x -> x .∈ Ref(ipm_ids))
        else
            tables[tab_name] = df
        end
    end

    return IPM.PadrinoDB(tables)
end

"""
    _pdb_species(pdb::PadrinoDB; ipm_id=nothing) -> Vector{String}

Return species names from the database.
"""
function _pdb_species(pdb::IPM.PadrinoDB; ipm_id=nothing)
    md = _get_metadata(pdb)
    DataFrames.nrow(md) == 0 && return String[]

    if ipm_id !== nothing
        ids = ipm_id isa String ? [ipm_id] : ipm_id
        md = DataFrames.subset(md, :ipm_id => x -> x .∈ Ref(ids))
    end

    species = String[]
    for row in DataFrames.eachrow(md)
        sp = get(row, :species_accepted, missing)
        if !ismissing(sp)
            push!(species, string(sp))
        end
    end

    return unique(species)
end

"""
    _pdb_citations(pdb::PadrinoDB; ipm_id=nothing)

Return citation information.
"""
function _pdb_citations(pdb::IPM.PadrinoDB; ipm_id=nothing)
    md = _get_metadata(pdb)
    DataFrames.nrow(md) == 0 && return DataFrames.DataFrame()

    if ipm_id !== nothing
        ids = ipm_id isa String ? [ipm_id] : ipm_id
        md = DataFrames.subset(md, :ipm_id => x -> x .∈ Ref(ids))
    end

    cols = intersect(DataFrames.propertynames(md),
                     [:ipm_id, :species_accepted, :authors, :journal,
                      :pub_year, :apa_citation, :doi])
    return DataFrames.select(md, cols)
end

"""
    _pdb_metadata(pdb::PadrinoDB, column::Symbol; ipm_id=nothing)

Return values of a metadata column.
"""
function _pdb_metadata(pdb::IPM.PadrinoDB, column::Symbol; ipm_id=nothing)
    md = _get_metadata(pdb)
    DataFrames.nrow(md) == 0 && return []

    if ipm_id !== nothing
        ids = ipm_id isa String ? [ipm_id] : ipm_id
        md = DataFrames.subset(md, :ipm_id => x -> x .∈ Ref(ids))
    end

    if column in DataFrames.propertynames(md)
        return md[!, column]
    else
        error("Column :$column not found in Metadata. Available: $(DataFrames.propertynames(md))")
    end
end

"""
    _pdb_test_targets(pdb::PadrinoDB; ipm_id=nothing)

Return test target values for validation.
PADRINO stores lambda and other targets in the Metadata table.
"""
function _pdb_test_targets(pdb::IPM.PadrinoDB; ipm_id=nothing)
    md = _get_metadata(pdb)
    DataFrames.nrow(md) == 0 && return DataFrames.DataFrame()

    if ipm_id !== nothing
        ids = ipm_id isa String ? [ipm_id] : ipm_id
        md = DataFrames.subset(md, :ipm_id => x -> x .∈ Ref(ids))
    end

    # Test targets are typically stored as lambda values
    target_cols = intersect(DataFrames.propertynames(md),
                            [:ipm_id, :species_accepted, :lambda])
    if !isempty(target_cols)
        return DataFrames.select(md, target_cols)
    end

    return DataFrames.select(md, [:ipm_id])
end

# --- Internal helpers ---

function _get_metadata(pdb::IPM.PadrinoDB)
    return get(pdb.tables, "Metadata", DataFrames.DataFrame())
end

function _get_table(pdb::IPM.PadrinoDB, name::String)
    return get(pdb.tables, name, DataFrames.DataFrame())
end

"""
    _get_ipm_ids(pdb::PadrinoDB) -> Vector{String}

Return all unique IPM IDs in the database.
"""
function _get_ipm_ids(pdb::IPM.PadrinoDB)
    md = _get_metadata(pdb)
    DataFrames.nrow(md) == 0 && return String[]
    return String[string(x) for x in unique(md.ipm_id)]
end
