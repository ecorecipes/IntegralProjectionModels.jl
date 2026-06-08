# PADRINO Database Integration

## Introduction

The [PADRINO](https://padrinodb.github.io/) database contains ~280 published Integral Projection Models stored as tabular data on GitHub. IntegralProjectionModels.jl provides an extension that downloads, parses, and builds IPMs directly from PADRINO — analogous to the R package [RPadrino](https://github.com/padrinoDB/RPadrino).

The extension is activated by loading `CSV`, `DataFrames`, and `Downloads`:

```@example ipm
using IntegralProjectionModels
using CSV, DataFrames, Downloads

```

## Downloading the Database

```@example ipm
pdb = pdb_download()
md = pdb.tables["Metadata"]
println("Models in database: ", nrow(md))
println("Unique species: ", length(unique(md.species_accepted)))
```

## Exploring the Database

Query species names:

```@example ipm
species = pdb_species(pdb)
first(species, 10)
```

Get citations for a specific model:

```@example ipm
cites = pdb_citations(pdb; ipm_id="aaa310")
cites
```

Get metadata columns:

```@example ipm
# Countries represented
countries = pdb_metadata(pdb, :country)
unique(skipmissing(countries))
```

## Building a Simple Deterministic Model

Model `aaa310` is a simple density-independent deterministic IPM. We can build it in two steps (proto → IPM) or one:

```@example ipm
# Two-step approach
protos = pdb_make_proto_ipm(pdb; ipm_id="aaa310")
pm = protos["aaa310"]
println("Species: ", pm.metadata[:species_accepted])
println("State variable: ", pm.state_variables[1][:state_variable])
println("Domain: [", pm.continuous_domains[1][:lower], ", ", pm.continuous_domains[1][:upper], "]")
println("Kernels: ", [k[:kernel_id] for k in pm.kernels])
println("Parameters: ", length(pm.parameters))
```

```@example ipm
# Build and solve
ipms = pdb_make_ipm(protos; tspan=(0, 100))
sol = solve(ipms["aaa310"])
lam = lambda(sol)
println("λ = ", round(lam, digits=4))
```

The one-step convenience method:

```@example ipm
ipms = pdb_make_ipm(pdb; ipm_id="aaa310", tspan=(0, 100))
sol = solve(ipms["aaa310"])
println("λ = ", round(lambda(sol), digits=4))
```

## Batch Analysis

Build and analyze multiple models at once:

```@example ipm
# Get all model IDs
all_ids = String.(pdb.tables["Metadata"].ipm_id)

# Build first 20 models
protos = pdb_make_proto_ipm(pdb; ipm_id=all_ids[1:20])
ipms = pdb_make_ipm(protos; tspan=(0, 50))
println("Successfully built: ", length(ipms), " / 20")
```

```@example ipm
# Compute lambda for each
results = Dict{String, Float64}()
for (id, ipm) in ipms
    try
        sol = solve(ipm)
        lam = lambda(sol)
        if isfinite(lam)
            results[id] = lam
        end
    catch
    end
end

for (id, lam) in sort(collect(results))
    println("  ", id, ": λ = ", round(lam, digits=4))
end
```

## Subsetting the Database

You can subset the database to work with specific models:

```@example ipm
sub_pdb = pdb_subset(pdb, ["aaa310", "aaa341"])
println("Models in subset: ", nrow(sub_pdb.tables["Metadata"]))
species = pdb_species(sub_pdb)
println("Species: ", species)
```

## Saving and Loading

Save the database locally for offline use:

```@example ipm
save_dir = joinpath(pwd(), "build", "padrino_data")
mkpath(save_dir)
pdb_save(pdb, save_dir)
pdb2 = pdb_load(save_dir)
rm(save_dir; recursive=true, force=true)
```

## Notes

- **Expression parsing**: PADRINO stores vital rates as R expressions. The extension includes a complete R expression parser that translates these to Julia functions.
- **Model types**: The extension handles simple deterministic, stochastic kernel-resampled, stochastic parameter-resampled, density-dependent, and general (multi-state) models.
- **Performance**: Kernel matrices are computed using the midpoint rule, matching the ipmr R package approach.
- **Eviction**: Truncated distribution and discrete extrema eviction corrections are supported.
