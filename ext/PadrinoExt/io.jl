"""
PADRINO database I/O functions: download, load, save.
"""

"""
    _pdb_download(; save=false, destination=nothing) -> PadrinoDB

Download the PADRINO database tables from GitHub.
"""
function _pdb_download(; save::Bool=false, destination=nothing)
    if save && destination === nothing
        error("'destination' must be specified when save=true")
    end

    tables = Dict{String, Any}()

    for tab_name in PADRINO_TABLES
        url = PADRINO_BASE_URL * tab_name * ".txt"
        @info "Downloading $tab_name..."
        try
            data = Downloads.download(url)
            df = CSV.read(data, DataFrames.DataFrame;
                          delim='\t',
                          missingstring=["NA", ""],
                          quotechar='"',
                          types=_column_types(tab_name))
            tables[tab_name] = df
        catch e
            @warn "Failed to download $tab_name" exception=e
            tables[tab_name] = DataFrames.DataFrame()
        end
    end

    pdb = IPM.PadrinoDB(tables)

    if save
        _pdb_save(pdb, destination)
    end

    return pdb
end

"""
    _pdb_load(path::String) -> PadrinoDB

Load PADRINO database tables from a directory of tab-delimited text files.
"""
function _pdb_load(path::String)
    # Ensure trailing separator
    if !endswith(path, "/") && !endswith(path, "\\")
        path = path * "/"
    end

    tables = Dict{String, Any}()

    for tab_name in PADRINO_TABLES
        filepath = path * tab_name * ".txt"
        if isfile(filepath)
            df = CSV.read(filepath, DataFrames.DataFrame;
                          delim='\t',
                          missingstring=["NA", ""],
                          quotechar='"',
                          types=_column_types(tab_name))
            tables[tab_name] = df
        else
            @warn "Table file not found: $filepath"
            tables[tab_name] = DataFrames.DataFrame()
        end
    end

    return IPM.PadrinoDB(tables)
end

"""
    _pdb_save(pdb::PadrinoDB, destination::String)

Save PADRINO database tables as tab-delimited text files.
"""
function _pdb_save(pdb::IPM.PadrinoDB, destination::String)
    if !endswith(destination, "/") && !endswith(destination, "\\")
        destination = destination * "/"
    end

    mkpath(destination)

    for tab_name in PADRINO_TABLES
        filepath = destination * tab_name * ".txt"
        if haskey(pdb.tables, tab_name)
            CSV.write(filepath, pdb.tables[tab_name];
                      delim='\t',
                      missingstring="NA",
                      quotechar='"')
        end
    end
end

"""
    _column_types(tab_name::String) -> Dict

Return column type hints for CSV parsing of specific PADRINO tables.
"""
function _column_types(tab_name::String)
    if tab_name == "ParameterValues"
        return Dict(:parameter_value => Float64)
    elseif tab_name == "ContinuousDomains"
        return Dict(:lower => Float64, :upper => Float64)
    elseif tab_name == "StateVectors"
        return Dict(:n_bins => Int)
    end
    return Dict{Symbol, Type}()
end
