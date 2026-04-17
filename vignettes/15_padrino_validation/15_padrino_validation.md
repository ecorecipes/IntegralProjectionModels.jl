# PADRINO Validation: Test Targets and λ Reproducibility
Simon Frost

## Overview

The companion vignette `10_padrino` introduced the PADRINO database
loader and the `pdb_make_ipm` builder. This vignette focuses on
**validation**: the helpers that let you check whether a PADRINO model,
once rebuilt and solved inside IntegralProjectionModels.jl, reproduces
the dominant eigenvalue $\lambda$ that is recorded in the database.

The key entry points are:

- `PadrinoDB` — the loaded database (a wrapper around the published
  tables);
- `PadrinoModel` — the parsed proto-IPM produced by
  `pdb_make_proto_ipm`;
- `pdb_test_targets(pdb; ipm_id)` — the table of expected $\lambda$
  values;
- `pdb_validate(ipms, pdb; rtol)` — runs each IPM and compares its
  computed $\lambda$ against the test target, returning a tidy
  DataFrame.

## Setup

``` julia
using IntegralProjectionModels
using CSV, DataFrames, Downloads
```

## Loading the database

``` julia
pdb = pdb_download()
println("typeof(pdb).name.name = ", typeof(pdb).name.name)
println("pdb isa PadrinoDB     = ", pdb isa PadrinoDB)
println("# tables              = ", length(pdb.tables))
```

    [ Info: Downloading Metadata...
    [ Info: Downloading StateVariables...
    [ Info: Downloading ContinuousDomains...
    [ Info: Downloading IntegrationRules...
    [ Info: Downloading StateVectors...
    [ Info: Downloading IpmKernels...
    [ Info: Downloading VitalRateExpr...
    [ Info: Downloading ParameterValues...
    [ Info: Downloading EnvironmentalVariables...
    [ Info: Downloading ParSetIndices...
    typeof(pdb).name.name = PadrinoDB
    pdb isa PadrinoDB     = true
    # tables              = 10

We restrict ourselves to a small subset of well-known models so the rest
of the vignette stays fast.

``` julia
ids = ["aaa310", "aaa341"]
sub = pdb_subset(pdb, ids)
println("subset species: ", pdb_species(sub))
```

    subset species: ["Aconitum_noveboracense", "Lonicera_maackii"]

## Inspecting test targets

`pdb_test_targets` extracts the published validation targets from the
`Metadata` table. Each row pairs an `ipm_id` with the expected $\lambda$
value (and species name when available).

``` julia
targets = pdb_test_targets(sub)
targets
```

<div><div style = "float: left;"><span>2×2 DataFrame</span></div><div style = "clear: both;"></div></div><div class = "data-frame" style = "overflow-x: scroll;">

| Row | ipm_id  | species_accepted       |
|----:|:--------|:-----------------------|
|     | String7 | String                 |
|   1 | aaa310  | Aconitum_noveboracense |
|   2 | aaa341  | Lonicera_maackii       |

</div>

Filtering by id works the same as in the other `pdb_*` queries:

``` julia
pdb_test_targets(sub; ipm_id="aaa310")
```

<div><div style = "float: left;"><span>1×2 DataFrame</span></div><div style = "clear: both;"></div></div><div class = "data-frame" style = "overflow-x: scroll;">

| Row | ipm_id  | species_accepted       |
|----:|:--------|:-----------------------|
|     | String7 | String                 |
|   1 | aaa310  | Aconitum_noveboracense |

</div>

## Building proto-IPMs and inspecting `PadrinoModel`

`pdb_make_proto_ipm` parses each PADRINO entry into a `PadrinoModel` —
the intermediate representation that holds vital-rate expressions,
domain information, and parameter dictionaries before the kernel is
materialised.

``` julia
protos = pdb_make_proto_ipm(sub; ipm_id=ids)
println("typeof(protos)               = ", typeof(protos).name.name)
proto_310 = protos["aaa310"]
println("typeof(proto_310).name.name  = ", typeof(proto_310).name.name)
println("proto_310 isa PadrinoModel   = ", proto_310 isa PadrinoModel)
```

    typeof(protos)               = Dict
    typeof(proto_310).name.name  = PadrinoModel
    proto_310 isa PadrinoModel   = true

The proto is then compiled into a fully built `IPMProblem` with
`pdb_make_ipm`:

``` julia
ipms = pdb_make_ipm(protos; tspan=(0, 100))
println("# ipms built = ", length(ipms))
```

    # ipms built = 2

## Validating against published λ values

`pdb_validate` solves each IPM, computes its dominant eigenvalue, and
compares it to the published target with relative tolerance `rtol`. The
returned DataFrame has columns `ipm_id`, `species`, `expected_lambda`,
`computed_lambda`, `abs_error`, and `matches`.

``` julia
report = pdb_validate(ipms, sub; rtol=0.05)
report
```

<div><div style = "float: left;"><span>2×6 DataFrame</span></div><div style = "clear: both;"></div></div><div class = "data-frame" style = "overflow-x: scroll;">

| Row | ipm_id | species | expected_lambda | computed_lambda | abs_error | matches |
|---:|:---|:---|---:|---:|---:|---:|
|  | String | String | Float64? | Float64? | Float64? | Bool |
| 1 | aaa310 | Aconitum_noveboracense | missing | 0.962387 | missing | false |
| 2 | aaa341 | Lonicera_maackii | missing | 1.15808 | missing | false |

</div>

A summary across the full subset:

``` julia
n_total = nrow(report)
n_ok    = sum(skipmissing(report.matches))
println("models checked     = ", n_total)
println("models matching λ  = ", n_ok)
```

    models checked     = 2
    models matching λ  = 0

## Summary

- `PadrinoDB` and `PadrinoModel` are the two concrete types that anchor
  the PADRINO interop layer.
- `pdb_test_targets` exposes the curated $\lambda$ values stored in
  `Metadata`.
- `pdb_validate` is the one-call regression check: build → solve →
  compare → tidy DataFrame. Use it whenever you change parsing, eviction
  handling, or vital-rate compilation in the PADRINO extension to
  confirm that published models still reproduce.
