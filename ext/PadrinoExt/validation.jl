"""
Validation of built IPMs against PADRINO test targets.
"""

"""
    _pdb_validate(ipms::Dict, pdb::PadrinoDB; rtol=0.01) -> DataFrame

Validate IPM results against PADRINO test targets (expected lambda values).
Returns a DataFrame with columns: ipm_id, species, expected_lambda, computed_lambda, matches.
"""
function _pdb_validate(ipms::Dict, pdb::IPM.PadrinoDB; rtol=0.01)
    md = _get_metadata(pdb)
    DataFrames.nrow(md) == 0 && return DataFrames.DataFrame()

    results = DataFrames.DataFrame(
        ipm_id = String[],
        species = String[],
        expected_lambda = Union{Float64, Missing}[],
        computed_lambda = Union{Float64, Missing}[],
        abs_error = Union{Float64, Missing}[],
        matches = Bool[]
    )

    for (ipm_id, ipm) in ipms
        # Get expected lambda from metadata
        rows = DataFrames.subset(md, :ipm_id => x -> x .== ipm_id)
        DataFrames.nrow(rows) == 0 && continue

        species = get(rows[1, :], :species_accepted, missing)
        species = ismissing(species) ? "" : string(species)

        # Look for lambda in metadata (column name varies)
        expected = missing
        for col in [:lambda, :target_lambda, :lam]
            if col in DataFrames.propertynames(rows)
                val = rows[1, col]
                if !ismissing(val)
                    expected = Float64(val)
                    break
                end
            end
        end

        # Compute lambda from the IPM
        computed = missing
        try
            if ipm isa IPM.IPMProblem
                sol = IPM.solve(ipm)
                computed = IPM.lambda(sol)
            elseif ipm isa IPM.IPMSolution
                computed = IPM.lambda(ipm)
            end
        catch e
            @warn "Failed to compute lambda for $ipm_id" exception=e
        end

        # Check match
        abs_err = missing
        matches = false
        if !ismissing(expected) && !ismissing(computed)
            abs_err = abs(expected - computed)
            matches = isapprox(expected, computed; rtol=rtol)
        end

        push!(results, (ipm_id=ipm_id, species=species,
                        expected_lambda=expected, computed_lambda=computed,
                        abs_error=abs_err, matches=matches))
    end

    return results
end
