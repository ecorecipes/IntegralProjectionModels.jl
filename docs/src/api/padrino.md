# PADRINO Integration

Interface to the PADRINO IPM database, providing access to a large collection of published IPMs. Requires the CSV, DataFrames, and Downloads packages to activate the package extension.

## Database Types

```@docs
PadrinoDB
PadrinoModel
```

## Loading & Saving

```@docs
pdb_download
pdb_load
pdb_save
```

## Querying

```@docs
pdb_subset
pdb_species
pdb_citations
pdb_metadata
pdb_test_targets
```

## Building Models

```@docs
pdb_make_proto_ipm
pdb_make_ipm
pdb_validate
```
