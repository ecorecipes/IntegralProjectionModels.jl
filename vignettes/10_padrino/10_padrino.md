# PADRINO Database Integration


## Introduction

The [PADRINO](https://padrinodb.github.io/) database contains ~280
published Integral Projection Models stored as tabular data on GitHub.
IntegralProjectionModels.jl provides an extension that downloads,
parses, and builds IPMs directly from PADRINO — analogous to the R
package [RPadrino](https://github.com/padrinoDB/RPadrino).

The extension is activated by loading `CSV`, `DataFrames`, and
`Downloads`:

``` julia
using IntegralProjectionModels
using CSV, DataFrames, Downloads
```

## Downloading the Database

``` julia
pdb = pdb_download()
md = pdb.tables["Metadata"]
println("Models in database: ", nrow(md))
println("Unique species: ", length(unique(md.species_accepted)))
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
    Models in database: 280
    Unique species: 56

## Exploring the Database

Query species names:

``` julia
species = pdb_species(pdb)
first(species, 10)
```

    10-element Vector{String}:
     "Calathea_crotalifera"
     "Heliconia_tortuosa"
     "Geum_radiatum"
     "Tilandsia_macdougallii"
     "Dendrophylax_lindenii"
     "Cecropia_obtusifolia"
     "Simarouba_amara"
     "Minquartia_guianensis"
     "Balizia_elegans"
     "Hymenolobium_mesoamericanum"

Get citations for a specific model:

``` julia
cites = pdb_citations(pdb; ipm_id="aaa310")
cites
```

<div><div style = "float: left;"><span>1×7 DataFrame</span></div><div style = "clear: both;"></div></div><div class = "data-frame" style = "overflow-x: scroll;">

| Row | ipm_id | species_accepted | authors | journal | pub_year | doi | apa_citation |
|---:|:---|:---|:---|:---|---:|:---|:---|
|  | String7 | String | String | String | Int64 | String? | String |
| 1 | aaa310 | Aconitum_noveboracense | Easterling; Ellner; Dixon | Ecology | 2000 | 10.1890/0012-9658(2000)081\[0694:SSSAAN\]2.0.CO;2 | Easterling, M. R., Ellner, S. P., & Dixon, P. M. (2000). Size-specific sensitivity: applying a new structured population model. Ecology, 81(3), 694-708. |

</div>

Get metadata columns:

``` julia
# Countries represented
countries = pdb_metadata(pdb, :country)
unique(skipmissing(countries))
```

    16-element Vector{String7}:
     "CRI"
     "USA"
     "MEX"
     "CUB"
     "GBR"
     "NZL  "
     "CHN"
     "FRA"
     "DEU"
     "CHE"
     "NOR"
     "ISR"
     "CAN"
     "PYF"
     "ESP"
     "JPN"

## Building a Simple Deterministic Model

Model `aaa310` is a simple density-independent deterministic IPM. We can
build it in two steps (proto → IPM) or one:

``` julia
# Two-step approach
protos = pdb_make_proto_ipm(pdb; ipm_id="aaa310")
pm = protos["aaa310"]
println("Species: ", pm.metadata[:species_accepted])
println("State variable: ", pm.state_variables[1][:state_variable])
println("Domain: [", pm.continuous_domains[1][:lower], ", ", pm.continuous_domains[1][:upper], "]")
println("Kernels: ", [k[:kernel_id] for k in pm.kernels])
println("Parameters: ", length(pm.parameters))
```

    Species: Aconitum_noveboracense
    State variable: size
    Domain: [0.0, 5.83]
    Kernels: String31["P", "F"]
    Parameters: 12

``` julia
# Build and solve
ipms = pdb_make_ipm(protos; tspan=(0, 100))
sol = solve(ipms["aaa310"])
lam = lambda(sol)
println("λ = ", round(lam, digits=4))
```

    λ = 0.9624

The one-step convenience method:

``` julia
ipms = pdb_make_ipm(pdb; ipm_id="aaa310", tspan=(0, 100))
sol = solve(ipms["aaa310"])
println("λ = ", round(lambda(sol), digits=4))
```

    λ = 0.9624

## Batch Analysis

Build and analyze multiple models at once:

``` julia
# Get all model IDs
all_ids = String.(pdb.tables["Metadata"].ipm_id)

# Build first 20 models
protos = pdb_make_proto_ipm(pdb; ipm_id=all_ids[1:20])
ipms = pdb_make_ipm(protos; tspan=(0, 50))
println("Successfully built: ", length(ipms), " / 20")
```

    ┌ Warning: Failed to build IPM for aaaa32
    │   exception =
    │    UndefVarError: `gi_le` not defined in `IntegralProjectionModelsPadrinoExt.PadrinoEvalModule`
    │    Suggestion: check for spelling errors or missing imports.
    │    Stacktrace:
    │      [1] (::IntegralProjectionModelsPadrinoExt.PadrinoEvalModule.var"#2#3")(z_new::Float64, z::Float64, h::Float64)
    │        @ IntegralProjectionModelsPadrinoExt.PadrinoEvalModule ~/Projects/ipm/IntegralProjectionModels.jl/ext/PadrinoExt/compiler.jl:484
    │      [2] build_kernel_matrix(kernel_fn::Function, domain::ContinuousDomain{Float64}, family::String3)
    │        @ IntegralProjectionModelsPadrinoExt ~/Projects/ipm/IntegralProjectionModels.jl/ext/PadrinoExt/compiler.jl:549
    │      [3] build_kernel_matrix(kernel_fn::Expr, vectorized_fn::Nothing, domain::ContinuousDomain{Float64}, family::String3)
    │        @ IntegralProjectionModelsPadrinoExt ~/Projects/ipm/IntegralProjectionModels.jl/ext/PadrinoExt/compiler.jl:609
    │      [4] _build_general_det(pm::PadrinoModel, tspan::Tuple{Int64, Int64}, mtype::IntegralProjectionModelsPadrinoExt.ModelTypeInfo)
    │        @ IntegralProjectionModelsPadrinoExt ~/Projects/ipm/IntegralProjectionModels.jl/ext/PadrinoExt/builder.jl:482
    │      [5] build_ipm(pm::PadrinoModel; tspan::Tuple{Int64, Int64})
    │        @ IntegralProjectionModelsPadrinoExt ~/Projects/ipm/IntegralProjectionModels.jl/ext/PadrinoExt/builder.jl:27
    │      [6] _pdb_make_ipm_from_protos(protos::Dict{String, PadrinoModel}; tspan::Tuple{Int64, Int64})
    │        @ IntegralProjectionModelsPadrinoExt ~/Projects/ipm/IntegralProjectionModels.jl/ext/IntegralProjectionModelsPadrinoExt.jl:337
    │      [7] _pdb_make_ipm_from_protos
    │        @ ~/Projects/ipm/IntegralProjectionModels.jl/ext/IntegralProjectionModelsPadrinoExt.jl:337 [inlined]
    │      [8] #pdb_make_ipm#69
    │        @ ~/Projects/ipm/IntegralProjectionModels.jl/ext/IntegralProjectionModelsPadrinoExt.jl:80 [inlined]
    │      [9] top-level scope
    │        @ ~/Projects/ipm/IntegralProjectionModels.jl/vignettes/10_padrino/10_padrino.qmd:95
    │     [10] eval(m::Module, e::Any)
    │        @ Core ./boot.jl:489
    │     [11] (::QuartoNotebookWorker.var"#21#22"{Module, Expr})()
    │        @ QuartoNotebookWorker ~/.julia/packages/QuartoNotebookRunner/evCNi/src/QuartoNotebookWorker/src/render.jl:222
    │     [12] (::QuartoNotebookWorker.Packages.IOCapture.var"#12#13"{Type{InterruptException}, QuartoNotebookWorker.var"#21#22"{Module, Expr}, IOContext{Base.PipeEndpoint}, IOContext{Base.PipeEndpoint}, IOContext{Base.PipeEndpoint}, IOContext{Base.PipeEndpoint}})()
    │        @ QuartoNotebookWorker.Packages.IOCapture ~/.julia/packages/QuartoNotebookRunner/evCNi/src/QuartoNotebookWorker/src/vendor/IOCapture/src/IOCapture.jl:170
    │     [13] with_logstate(f::QuartoNotebookWorker.Packages.IOCapture.var"#12#13"{Type{InterruptException}, QuartoNotebookWorker.var"#21#22"{Module, Expr}, IOContext{Base.PipeEndpoint}, IOContext{Base.PipeEndpoint}, IOContext{Base.PipeEndpoint}, IOContext{Base.PipeEndpoint}}, logstate::Base.CoreLogging.LogState)
    │        @ Base.CoreLogging ./logging/logging.jl:542
    │     [14] with_logger(f::Function, logger::Base.CoreLogging.ConsoleLogger)
    │        @ Base.CoreLogging ./logging/logging.jl:653
    │     [15] capture(f::QuartoNotebookWorker.var"#21#22"{Module, Expr}; rethrow::Type, color::Bool, passthrough::Bool, capture_buffer::IOBuffer, io_context::Vector{Pair{Symbol, Any}})
    │        @ QuartoNotebookWorker.Packages.IOCapture ~/.julia/packages/QuartoNotebookRunner/evCNi/src/QuartoNotebookWorker/src/vendor/IOCapture/src/IOCapture.jl:167
    │     [16] capture
    │        @ ~/.julia/packages/QuartoNotebookRunner/evCNi/src/QuartoNotebookWorker/src/render.jl:248 [inlined]
    │     [17] io_capture(f::Function; cell_options::Dict{String, Any}, kws::@Kwargs{rethrow::DataType, color::Bool, io_context::Vector{Pair{Symbol, Any}}})
    │        @ QuartoNotebookWorker ~/.julia/packages/QuartoNotebookRunner/evCNi/src/QuartoNotebookWorker/src/render.jl:250
    │     [18] io_capture
    │        @ ~/.julia/packages/QuartoNotebookRunner/evCNi/src/QuartoNotebookWorker/src/render.jl:246 [inlined]
    │     [19] include_str(mod::Module, code::String; file::String, line::Int64, cell_options::Dict{String, Any})
    │        @ QuartoNotebookWorker ~/.julia/packages/QuartoNotebookRunner/evCNi/src/QuartoNotebookWorker/src/render.jl:201
    │     [20] #invokelatest_gr#232
    │        @ ./reflection.jl:1297 [inlined]
    │     [21] invokelatest_gr
    │        @ ./reflection.jl:1289 [inlined]
    │     [22] #7
    │        @ ~/.julia/packages/QuartoNotebookRunner/evCNi/src/QuartoNotebookWorker/src/render.jl:18 [inlined]
    │     [23] with_inline_display(f::QuartoNotebookWorker.var"#7#8"{String, String, Int64, Dict{String, Any}}, cell_options::Dict{String, Any})
    │        @ QuartoNotebookWorker ~/.julia/packages/QuartoNotebookRunner/evCNi/src/QuartoNotebookWorker/src/InlineDisplay.jl:31
    │     [24] _render_thunk(thunk::Function, code::String, cell_options::Dict{String, Any}, is_expansion_ref::Base.RefValue{Bool}; inline::Bool)
    │        @ QuartoNotebookWorker ~/.julia/packages/QuartoNotebookRunner/evCNi/src/QuartoNotebookWorker/src/render.jl:43
    │     [25] _render_thunk
    │        @ ~/.julia/packages/QuartoNotebookRunner/evCNi/src/QuartoNotebookWorker/src/render.jl:35 [inlined]
    │     [26] #render#4
    │        @ ~/.julia/packages/QuartoNotebookRunner/evCNi/src/QuartoNotebookWorker/src/render.jl:15 [inlined]
    │     [27] render(code::String, file::String, line::Int64, cell_options::Dict{String, Any})
    │        @ QuartoNotebookWorker ~/.julia/packages/QuartoNotebookRunner/evCNi/src/QuartoNotebookWorker/src/render.jl:1
    │     [28] render(::String, ::Vararg{Any}; kwargs::@Kwargs{})
    │        @ Main ~/.julia/packages/QuartoNotebookRunner/evCNi/src/startup.jl:145
    │     [29] top-level scope
    │        @ none:1
    │     [30] eval(m::Module, e::Any)
    │        @ Core ./boot.jl:489
    │     [31] (::Main.var"#22#23")(chan::Channel{Any})
    │        @ Main ~/.julia/packages/QuartoNotebookRunner/evCNi/src/startup.jl:158
    │     [32] (::Base.var"#562#563"{Main.var"#22#23", Channel{Any}})()
    │        @ Base ./channels.jl:141
    └ @ IntegralProjectionModelsPadrinoExt ~/Projects/ipm/IntegralProjectionModels.jl/ext/IntegralProjectionModelsPadrinoExt.jl:345
    ┌ Warning: Failed to build IPM for aaaa27
    │   exception =
    │    UndefVarError: `gi_le` not defined in `IntegralProjectionModelsPadrinoExt.PadrinoEvalModule`
    │    Suggestion: check for spelling errors or missing imports.
    │    Stacktrace:
    │      [1] (::IntegralProjectionModelsPadrinoExt.PadrinoEvalModule.var"#5#6")(z_new::Float64, z::Float64, h::Float64)
    │        @ IntegralProjectionModelsPadrinoExt.PadrinoEvalModule ~/Projects/ipm/IntegralProjectionModels.jl/ext/PadrinoExt/compiler.jl:484
    │      [2] build_kernel_matrix(kernel_fn::Function, domain::ContinuousDomain{Float64}, family::String3)
    │        @ IntegralProjectionModelsPadrinoExt ~/Projects/ipm/IntegralProjectionModels.jl/ext/PadrinoExt/compiler.jl:549
    │      [3] build_kernel_matrix(kernel_fn::Expr, vectorized_fn::Nothing, domain::ContinuousDomain{Float64}, family::String3)
    │        @ IntegralProjectionModelsPadrinoExt ~/Projects/ipm/IntegralProjectionModels.jl/ext/PadrinoExt/compiler.jl:609
    │      [4] _build_general_det(pm::PadrinoModel, tspan::Tuple{Int64, Int64}, mtype::IntegralProjectionModelsPadrinoExt.ModelTypeInfo)
    │        @ IntegralProjectionModelsPadrinoExt ~/Projects/ipm/IntegralProjectionModels.jl/ext/PadrinoExt/builder.jl:482
    │      [5] build_ipm(pm::PadrinoModel; tspan::Tuple{Int64, Int64})
    │        @ IntegralProjectionModelsPadrinoExt ~/Projects/ipm/IntegralProjectionModels.jl/ext/PadrinoExt/builder.jl:27
    │      [6] _pdb_make_ipm_from_protos(protos::Dict{String, PadrinoModel}; tspan::Tuple{Int64, Int64})
    │        @ IntegralProjectionModelsPadrinoExt ~/Projects/ipm/IntegralProjectionModels.jl/ext/IntegralProjectionModelsPadrinoExt.jl:337
    │      [7] _pdb_make_ipm_from_protos
    │        @ ~/Projects/ipm/IntegralProjectionModels.jl/ext/IntegralProjectionModelsPadrinoExt.jl:337 [inlined]
    │      [8] #pdb_make_ipm#69
    │        @ ~/Projects/ipm/IntegralProjectionModels.jl/ext/IntegralProjectionModelsPadrinoExt.jl:80 [inlined]
    │      [9] top-level scope
    │        @ ~/Projects/ipm/IntegralProjectionModels.jl/vignettes/10_padrino/10_padrino.qmd:95
    │     [10] eval(m::Module, e::Any)
    │        @ Core ./boot.jl:489
    │     [11] (::QuartoNotebookWorker.var"#21#22"{Module, Expr})()
    │        @ QuartoNotebookWorker ~/.julia/packages/QuartoNotebookRunner/evCNi/src/QuartoNotebookWorker/src/render.jl:222
    │     [12] (::QuartoNotebookWorker.Packages.IOCapture.var"#12#13"{Type{InterruptException}, QuartoNotebookWorker.var"#21#22"{Module, Expr}, IOContext{Base.PipeEndpoint}, IOContext{Base.PipeEndpoint}, IOContext{Base.PipeEndpoint}, IOContext{Base.PipeEndpoint}})()
    │        @ QuartoNotebookWorker.Packages.IOCapture ~/.julia/packages/QuartoNotebookRunner/evCNi/src/QuartoNotebookWorker/src/vendor/IOCapture/src/IOCapture.jl:170
    │     [13] with_logstate(f::QuartoNotebookWorker.Packages.IOCapture.var"#12#13"{Type{InterruptException}, QuartoNotebookWorker.var"#21#22"{Module, Expr}, IOContext{Base.PipeEndpoint}, IOContext{Base.PipeEndpoint}, IOContext{Base.PipeEndpoint}, IOContext{Base.PipeEndpoint}}, logstate::Base.CoreLogging.LogState)
    │        @ Base.CoreLogging ./logging/logging.jl:542
    │     [14] with_logger(f::Function, logger::Base.CoreLogging.ConsoleLogger)
    │        @ Base.CoreLogging ./logging/logging.jl:653
    │     [15] capture(f::QuartoNotebookWorker.var"#21#22"{Module, Expr}; rethrow::Type, color::Bool, passthrough::Bool, capture_buffer::IOBuffer, io_context::Vector{Pair{Symbol, Any}})
    │        @ QuartoNotebookWorker.Packages.IOCapture ~/.julia/packages/QuartoNotebookRunner/evCNi/src/QuartoNotebookWorker/src/vendor/IOCapture/src/IOCapture.jl:167
    │     [16] capture
    │        @ ~/.julia/packages/QuartoNotebookRunner/evCNi/src/QuartoNotebookWorker/src/render.jl:248 [inlined]
    │     [17] io_capture(f::Function; cell_options::Dict{String, Any}, kws::@Kwargs{rethrow::DataType, color::Bool, io_context::Vector{Pair{Symbol, Any}}})
    │        @ QuartoNotebookWorker ~/.julia/packages/QuartoNotebookRunner/evCNi/src/QuartoNotebookWorker/src/render.jl:250
    │     [18] io_capture
    │        @ ~/.julia/packages/QuartoNotebookRunner/evCNi/src/QuartoNotebookWorker/src/render.jl:246 [inlined]
    │     [19] include_str(mod::Module, code::String; file::String, line::Int64, cell_options::Dict{String, Any})
    │        @ QuartoNotebookWorker ~/.julia/packages/QuartoNotebookRunner/evCNi/src/QuartoNotebookWorker/src/render.jl:201
    │     [20] #invokelatest_gr#232
    │        @ ./reflection.jl:1297 [inlined]
    │     [21] invokelatest_gr
    │        @ ./reflection.jl:1289 [inlined]
    │     [22] #7
    │        @ ~/.julia/packages/QuartoNotebookRunner/evCNi/src/QuartoNotebookWorker/src/render.jl:18 [inlined]
    │     [23] with_inline_display(f::QuartoNotebookWorker.var"#7#8"{String, String, Int64, Dict{String, Any}}, cell_options::Dict{String, Any})
    │        @ QuartoNotebookWorker ~/.julia/packages/QuartoNotebookRunner/evCNi/src/QuartoNotebookWorker/src/InlineDisplay.jl:31
    │     [24] _render_thunk(thunk::Function, code::String, cell_options::Dict{String, Any}, is_expansion_ref::Base.RefValue{Bool}; inline::Bool)
    │        @ QuartoNotebookWorker ~/.julia/packages/QuartoNotebookRunner/evCNi/src/QuartoNotebookWorker/src/render.jl:43
    │     [25] _render_thunk
    │        @ ~/.julia/packages/QuartoNotebookRunner/evCNi/src/QuartoNotebookWorker/src/render.jl:35 [inlined]
    │     [26] #render#4
    │        @ ~/.julia/packages/QuartoNotebookRunner/evCNi/src/QuartoNotebookWorker/src/render.jl:15 [inlined]
    │     [27] render(code::String, file::String, line::Int64, cell_options::Dict{String, Any})
    │        @ QuartoNotebookWorker ~/.julia/packages/QuartoNotebookRunner/evCNi/src/QuartoNotebookWorker/src/render.jl:1
    │     [28] render(::String, ::Vararg{Any}; kwargs::@Kwargs{})
    │        @ Main ~/.julia/packages/QuartoNotebookRunner/evCNi/src/startup.jl:145
    │     [29] top-level scope
    │        @ none:1
    │     [30] eval(m::Module, e::Any)
    │        @ Core ./boot.jl:489
    │     [31] (::Main.var"#22#23")(chan::Channel{Any})
    │        @ Main ~/.julia/packages/QuartoNotebookRunner/evCNi/src/startup.jl:158
    │     [32] (::Base.var"#562#563"{Main.var"#22#23", Channel{Any}})()
    │        @ Base ./channels.jl:141
    └ @ IntegralProjectionModelsPadrinoExt ~/Projects/ipm/IntegralProjectionModels.jl/ext/IntegralProjectionModelsPadrinoExt.jl:345
    ┌ Warning: Failed to build IPM for aaaa26
    │   exception =
    │    UndefVarError: `gi_le` not defined in `IntegralProjectionModelsPadrinoExt.PadrinoEvalModule`
    │    Suggestion: check for spelling errors or missing imports.
    │    Stacktrace:
    │      [1] (::IntegralProjectionModelsPadrinoExt.PadrinoEvalModule.var"#17#18")(z_new::Float64, z::Float64, h::Float64)
    │        @ IntegralProjectionModelsPadrinoExt.PadrinoEvalModule ~/Projects/ipm/IntegralProjectionModels.jl/ext/PadrinoExt/compiler.jl:484
    │      [2] build_kernel_matrix(kernel_fn::Function, domain::ContinuousDomain{Float64}, family::String3)
    │        @ IntegralProjectionModelsPadrinoExt ~/Projects/ipm/IntegralProjectionModels.jl/ext/PadrinoExt/compiler.jl:549
    │      [3] build_kernel_matrix(kernel_fn::Expr, vectorized_fn::Nothing, domain::ContinuousDomain{Float64}, family::String3)
    │        @ IntegralProjectionModelsPadrinoExt ~/Projects/ipm/IntegralProjectionModels.jl/ext/PadrinoExt/compiler.jl:609
    │      [4] _build_general_det(pm::PadrinoModel, tspan::Tuple{Int64, Int64}, mtype::IntegralProjectionModelsPadrinoExt.ModelTypeInfo)
    │        @ IntegralProjectionModelsPadrinoExt ~/Projects/ipm/IntegralProjectionModels.jl/ext/PadrinoExt/builder.jl:482
    │      [5] build_ipm(pm::PadrinoModel; tspan::Tuple{Int64, Int64})
    │        @ IntegralProjectionModelsPadrinoExt ~/Projects/ipm/IntegralProjectionModels.jl/ext/PadrinoExt/builder.jl:27
    │      [6] _pdb_make_ipm_from_protos(protos::Dict{String, PadrinoModel}; tspan::Tuple{Int64, Int64})
    │        @ IntegralProjectionModelsPadrinoExt ~/Projects/ipm/IntegralProjectionModels.jl/ext/IntegralProjectionModelsPadrinoExt.jl:337
    │      [7] _pdb_make_ipm_from_protos
    │        @ ~/Projects/ipm/IntegralProjectionModels.jl/ext/IntegralProjectionModelsPadrinoExt.jl:337 [inlined]
    │      [8] #pdb_make_ipm#69
    │        @ ~/Projects/ipm/IntegralProjectionModels.jl/ext/IntegralProjectionModelsPadrinoExt.jl:80 [inlined]
    │      [9] top-level scope
    │        @ ~/Projects/ipm/IntegralProjectionModels.jl/vignettes/10_padrino/10_padrino.qmd:95
    │     [10] eval(m::Module, e::Any)
    │        @ Core ./boot.jl:489
    │     [11] (::QuartoNotebookWorker.var"#21#22"{Module, Expr})()
    │        @ QuartoNotebookWorker ~/.julia/packages/QuartoNotebookRunner/evCNi/src/QuartoNotebookWorker/src/render.jl:222
    │     [12] (::QuartoNotebookWorker.Packages.IOCapture.var"#12#13"{Type{InterruptException}, QuartoNotebookWorker.var"#21#22"{Module, Expr}, IOContext{Base.PipeEndpoint}, IOContext{Base.PipeEndpoint}, IOContext{Base.PipeEndpoint}, IOContext{Base.PipeEndpoint}})()
    │        @ QuartoNotebookWorker.Packages.IOCapture ~/.julia/packages/QuartoNotebookRunner/evCNi/src/QuartoNotebookWorker/src/vendor/IOCapture/src/IOCapture.jl:170
    │     [13] with_logstate(f::QuartoNotebookWorker.Packages.IOCapture.var"#12#13"{Type{InterruptException}, QuartoNotebookWorker.var"#21#22"{Module, Expr}, IOContext{Base.PipeEndpoint}, IOContext{Base.PipeEndpoint}, IOContext{Base.PipeEndpoint}, IOContext{Base.PipeEndpoint}}, logstate::Base.CoreLogging.LogState)
    │        @ Base.CoreLogging ./logging/logging.jl:542
    │     [14] with_logger(f::Function, logger::Base.CoreLogging.ConsoleLogger)
    │        @ Base.CoreLogging ./logging/logging.jl:653
    │     [15] capture(f::QuartoNotebookWorker.var"#21#22"{Module, Expr}; rethrow::Type, color::Bool, passthrough::Bool, capture_buffer::IOBuffer, io_context::Vector{Pair{Symbol, Any}})
    │        @ QuartoNotebookWorker.Packages.IOCapture ~/.julia/packages/QuartoNotebookRunner/evCNi/src/QuartoNotebookWorker/src/vendor/IOCapture/src/IOCapture.jl:167
    │     [16] capture
    │        @ ~/.julia/packages/QuartoNotebookRunner/evCNi/src/QuartoNotebookWorker/src/render.jl:248 [inlined]
    │     [17] io_capture(f::Function; cell_options::Dict{String, Any}, kws::@Kwargs{rethrow::DataType, color::Bool, io_context::Vector{Pair{Symbol, Any}}})
    │        @ QuartoNotebookWorker ~/.julia/packages/QuartoNotebookRunner/evCNi/src/QuartoNotebookWorker/src/render.jl:250
    │     [18] io_capture
    │        @ ~/.julia/packages/QuartoNotebookRunner/evCNi/src/QuartoNotebookWorker/src/render.jl:246 [inlined]
    │     [19] include_str(mod::Module, code::String; file::String, line::Int64, cell_options::Dict{String, Any})
    │        @ QuartoNotebookWorker ~/.julia/packages/QuartoNotebookRunner/evCNi/src/QuartoNotebookWorker/src/render.jl:201
    │     [20] #invokelatest_gr#232
    │        @ ./reflection.jl:1297 [inlined]
    │     [21] invokelatest_gr
    │        @ ./reflection.jl:1289 [inlined]
    │     [22] #7
    │        @ ~/.julia/packages/QuartoNotebookRunner/evCNi/src/QuartoNotebookWorker/src/render.jl:18 [inlined]
    │     [23] with_inline_display(f::QuartoNotebookWorker.var"#7#8"{String, String, Int64, Dict{String, Any}}, cell_options::Dict{String, Any})
    │        @ QuartoNotebookWorker ~/.julia/packages/QuartoNotebookRunner/evCNi/src/QuartoNotebookWorker/src/InlineDisplay.jl:31
    │     [24] _render_thunk(thunk::Function, code::String, cell_options::Dict{String, Any}, is_expansion_ref::Base.RefValue{Bool}; inline::Bool)
    │        @ QuartoNotebookWorker ~/.julia/packages/QuartoNotebookRunner/evCNi/src/QuartoNotebookWorker/src/render.jl:43
    │     [25] _render_thunk
    │        @ ~/.julia/packages/QuartoNotebookRunner/evCNi/src/QuartoNotebookWorker/src/render.jl:35 [inlined]
    │     [26] #render#4
    │        @ ~/.julia/packages/QuartoNotebookRunner/evCNi/src/QuartoNotebookWorker/src/render.jl:15 [inlined]
    │     [27] render(code::String, file::String, line::Int64, cell_options::Dict{String, Any})
    │        @ QuartoNotebookWorker ~/.julia/packages/QuartoNotebookRunner/evCNi/src/QuartoNotebookWorker/src/render.jl:1
    │     [28] render(::String, ::Vararg{Any}; kwargs::@Kwargs{})
    │        @ Main ~/.julia/packages/QuartoNotebookRunner/evCNi/src/startup.jl:145
    │     [29] top-level scope
    │        @ none:1
    │     [30] eval(m::Module, e::Any)
    │        @ Core ./boot.jl:489
    │     [31] (::Main.var"#22#23")(chan::Channel{Any})
    │        @ Main ~/.julia/packages/QuartoNotebookRunner/evCNi/src/startup.jl:158
    │     [32] (::Base.var"#562#563"{Main.var"#22#23", Channel{Any}})()
    │        @ Base ./channels.jl:141
    └ @ IntegralProjectionModelsPadrinoExt ~/Projects/ipm/IntegralProjectionModels.jl/ext/IntegralProjectionModelsPadrinoExt.jl:345
    ┌ Warning: Failed to build env sampler for j
    │   exception =
    │    MethodError: no method matching translate_env_function(::String31, ::Vector{String})
    │    The function `translate_env_function` exists, but no method is defined for this combination of argument types.
    │    
    │    Closest candidates are:
    │      translate_env_function(::String, ::Vector{String})
    │       @ IntegralProjectionModelsPadrinoExt ~/Projects/ipm/IntegralProjectionModels.jl/ext/PadrinoExt/translator.jl:370
    │    
    └ @ IntegralProjectionModelsPadrinoExt ~/Projects/ipm/IntegralProjectionModels.jl/ext/PadrinoExt/env_state.jl:44
    ┌ Warning: Failed to build env sampler for A_max
    │   exception =
    │    MethodError: no method matching translate_env_function(::String31, ::Vector{String})
    │    The function `translate_env_function` exists, but no method is defined for this combination of argument types.
    │    
    │    Closest candidates are:
    │      translate_env_function(::String, ::Vector{String})
    │       @ IntegralProjectionModelsPadrinoExt ~/Projects/ipm/IntegralProjectionModels.jl/ext/PadrinoExt/translator.jl:370
    │    
    └ @ IntegralProjectionModelsPadrinoExt ~/Projects/ipm/IntegralProjectionModels.jl/ext/PadrinoExt/env_state.jl:44
    ┌ Warning: Failed to build env sampler for j
    │   exception =
    │    MethodError: no method matching translate_env_function(::String31, ::Vector{String})
    │    The function `translate_env_function` exists, but no method is defined for this combination of argument types.
    │    
    │    Closest candidates are:
    │      translate_env_function(::String, ::Vector{String})
    │       @ IntegralProjectionModelsPadrinoExt ~/Projects/ipm/IntegralProjectionModels.jl/ext/PadrinoExt/translator.jl:370
    │    
    └ @ IntegralProjectionModelsPadrinoExt ~/Projects/ipm/IntegralProjectionModels.jl/ext/PadrinoExt/env_state.jl:44
    ┌ Warning: Failed to build env sampler for A_max
    │   exception =
    │    MethodError: no method matching translate_env_function(::String31, ::Vector{String})
    │    The function `translate_env_function` exists, but no method is defined for this combination of argument types.
    │    
    │    Closest candidates are:
    │      translate_env_function(::String, ::Vector{String})
    │       @ IntegralProjectionModelsPadrinoExt ~/Projects/ipm/IntegralProjectionModels.jl/ext/PadrinoExt/translator.jl:370
    │    
    └ @ IntegralProjectionModelsPadrinoExt ~/Projects/ipm/IntegralProjectionModels.jl/ext/PadrinoExt/env_state.jl:44
    ┌ Warning: Failed to build IPM for aaaa24
    │   exception =
    │    UndefVarError: `gi_le` not defined in `IntegralProjectionModelsPadrinoExt.PadrinoEvalModule`
    │    Suggestion: check for spelling errors or missing imports.
    │    Stacktrace:
    │      [1] (::IntegralProjectionModelsPadrinoExt.PadrinoEvalModule.var"#29#30")(z_new::Float64, z::Float64, h::Float64)
    │        @ IntegralProjectionModelsPadrinoExt.PadrinoEvalModule ~/Projects/ipm/IntegralProjectionModels.jl/ext/PadrinoExt/compiler.jl:484
    │      [2] build_kernel_matrix(kernel_fn::Function, domain::ContinuousDomain{Float64}, family::String3)
    │        @ IntegralProjectionModelsPadrinoExt ~/Projects/ipm/IntegralProjectionModels.jl/ext/PadrinoExt/compiler.jl:549
    │      [3] build_kernel_matrix(kernel_fn::Expr, vectorized_fn::Nothing, domain::ContinuousDomain{Float64}, family::String3)
    │        @ IntegralProjectionModelsPadrinoExt ~/Projects/ipm/IntegralProjectionModels.jl/ext/PadrinoExt/compiler.jl:609
    │      [4] _build_general_det(pm::PadrinoModel, tspan::Tuple{Int64, Int64}, mtype::IntegralProjectionModelsPadrinoExt.ModelTypeInfo)
    │        @ IntegralProjectionModelsPadrinoExt ~/Projects/ipm/IntegralProjectionModels.jl/ext/PadrinoExt/builder.jl:482
    │      [5] build_ipm(pm::PadrinoModel; tspan::Tuple{Int64, Int64})
    │        @ IntegralProjectionModelsPadrinoExt ~/Projects/ipm/IntegralProjectionModels.jl/ext/PadrinoExt/builder.jl:27
    │      [6] _pdb_make_ipm_from_protos(protos::Dict{String, PadrinoModel}; tspan::Tuple{Int64, Int64})
    │        @ IntegralProjectionModelsPadrinoExt ~/Projects/ipm/IntegralProjectionModels.jl/ext/IntegralProjectionModelsPadrinoExt.jl:337
    │      [7] _pdb_make_ipm_from_protos
    │        @ ~/Projects/ipm/IntegralProjectionModels.jl/ext/IntegralProjectionModelsPadrinoExt.jl:337 [inlined]
    │      [8] #pdb_make_ipm#69
    │        @ ~/Projects/ipm/IntegralProjectionModels.jl/ext/IntegralProjectionModelsPadrinoExt.jl:80 [inlined]
    │      [9] top-level scope
    │        @ ~/Projects/ipm/IntegralProjectionModels.jl/vignettes/10_padrino/10_padrino.qmd:95
    │     [10] eval(m::Module, e::Any)
    │        @ Core ./boot.jl:489
    │     [11] (::QuartoNotebookWorker.var"#21#22"{Module, Expr})()
    │        @ QuartoNotebookWorker ~/.julia/packages/QuartoNotebookRunner/evCNi/src/QuartoNotebookWorker/src/render.jl:222
    │     [12] (::QuartoNotebookWorker.Packages.IOCapture.var"#12#13"{Type{InterruptException}, QuartoNotebookWorker.var"#21#22"{Module, Expr}, IOContext{Base.PipeEndpoint}, IOContext{Base.PipeEndpoint}, IOContext{Base.PipeEndpoint}, IOContext{Base.PipeEndpoint}})()
    │        @ QuartoNotebookWorker.Packages.IOCapture ~/.julia/packages/QuartoNotebookRunner/evCNi/src/QuartoNotebookWorker/src/vendor/IOCapture/src/IOCapture.jl:170
    │     [13] with_logstate(f::QuartoNotebookWorker.Packages.IOCapture.var"#12#13"{Type{InterruptException}, QuartoNotebookWorker.var"#21#22"{Module, Expr}, IOContext{Base.PipeEndpoint}, IOContext{Base.PipeEndpoint}, IOContext{Base.PipeEndpoint}, IOContext{Base.PipeEndpoint}}, logstate::Base.CoreLogging.LogState)
    │        @ Base.CoreLogging ./logging/logging.jl:542
    │     [14] with_logger(f::Function, logger::Base.CoreLogging.ConsoleLogger)
    │        @ Base.CoreLogging ./logging/logging.jl:653
    │     [15] capture(f::QuartoNotebookWorker.var"#21#22"{Module, Expr}; rethrow::Type, color::Bool, passthrough::Bool, capture_buffer::IOBuffer, io_context::Vector{Pair{Symbol, Any}})
    │        @ QuartoNotebookWorker.Packages.IOCapture ~/.julia/packages/QuartoNotebookRunner/evCNi/src/QuartoNotebookWorker/src/vendor/IOCapture/src/IOCapture.jl:167
    │     [16] capture
    │        @ ~/.julia/packages/QuartoNotebookRunner/evCNi/src/QuartoNotebookWorker/src/render.jl:248 [inlined]
    │     [17] io_capture(f::Function; cell_options::Dict{String, Any}, kws::@Kwargs{rethrow::DataType, color::Bool, io_context::Vector{Pair{Symbol, Any}}})
    │        @ QuartoNotebookWorker ~/.julia/packages/QuartoNotebookRunner/evCNi/src/QuartoNotebookWorker/src/render.jl:250
    │     [18] io_capture
    │        @ ~/.julia/packages/QuartoNotebookRunner/evCNi/src/QuartoNotebookWorker/src/render.jl:246 [inlined]
    │     [19] include_str(mod::Module, code::String; file::String, line::Int64, cell_options::Dict{String, Any})
    │        @ QuartoNotebookWorker ~/.julia/packages/QuartoNotebookRunner/evCNi/src/QuartoNotebookWorker/src/render.jl:201
    │     [20] #invokelatest_gr#232
    │        @ ./reflection.jl:1297 [inlined]
    │     [21] invokelatest_gr
    │        @ ./reflection.jl:1289 [inlined]
    │     [22] #7
    │        @ ~/.julia/packages/QuartoNotebookRunner/evCNi/src/QuartoNotebookWorker/src/render.jl:18 [inlined]
    │     [23] with_inline_display(f::QuartoNotebookWorker.var"#7#8"{String, String, Int64, Dict{String, Any}}, cell_options::Dict{String, Any})
    │        @ QuartoNotebookWorker ~/.julia/packages/QuartoNotebookRunner/evCNi/src/QuartoNotebookWorker/src/InlineDisplay.jl:31
    │     [24] _render_thunk(thunk::Function, code::String, cell_options::Dict{String, Any}, is_expansion_ref::Base.RefValue{Bool}; inline::Bool)
    │        @ QuartoNotebookWorker ~/.julia/packages/QuartoNotebookRunner/evCNi/src/QuartoNotebookWorker/src/render.jl:43
    │     [25] _render_thunk
    │        @ ~/.julia/packages/QuartoNotebookRunner/evCNi/src/QuartoNotebookWorker/src/render.jl:35 [inlined]
    │     [26] #render#4
    │        @ ~/.julia/packages/QuartoNotebookRunner/evCNi/src/QuartoNotebookWorker/src/render.jl:15 [inlined]
    │     [27] render(code::String, file::String, line::Int64, cell_options::Dict{String, Any})
    │        @ QuartoNotebookWorker ~/.julia/packages/QuartoNotebookRunner/evCNi/src/QuartoNotebookWorker/src/render.jl:1
    │     [28] render(::String, ::Vararg{Any}; kwargs::@Kwargs{})
    │        @ Main ~/.julia/packages/QuartoNotebookRunner/evCNi/src/startup.jl:145
    │     [29] top-level scope
    │        @ none:1
    │     [30] eval(m::Module, e::Any)
    │        @ Core ./boot.jl:489
    │     [31] (::Main.var"#22#23")(chan::Channel{Any})
    │        @ Main ~/.julia/packages/QuartoNotebookRunner/evCNi/src/startup.jl:158
    │     [32] (::Base.var"#562#563"{Main.var"#22#23", Channel{Any}})()
    │        @ Base ./channels.jl:141
    └ @ IntegralProjectionModelsPadrinoExt ~/Projects/ipm/IntegralProjectionModels.jl/ext/IntegralProjectionModelsPadrinoExt.jl:345
    ┌ Warning: Failed to build IPM for aaaa30
    │   exception =
    │    UndefVarError: `gi_le` not defined in `IntegralProjectionModelsPadrinoExt.PadrinoEvalModule`
    │    Suggestion: check for spelling errors or missing imports.
    │    Stacktrace:
    │      [1] (::IntegralProjectionModelsPadrinoExt.PadrinoEvalModule.var"#32#33")(z_new::Float64, z::Float64, h::Float64)
    │        @ IntegralProjectionModelsPadrinoExt.PadrinoEvalModule ~/Projects/ipm/IntegralProjectionModels.jl/ext/PadrinoExt/compiler.jl:484
    │      [2] build_kernel_matrix(kernel_fn::Function, domain::ContinuousDomain{Float64}, family::String3)
    │        @ IntegralProjectionModelsPadrinoExt ~/Projects/ipm/IntegralProjectionModels.jl/ext/PadrinoExt/compiler.jl:549
    │      [3] build_kernel_matrix(kernel_fn::Expr, vectorized_fn::Nothing, domain::ContinuousDomain{Float64}, family::String3)
    │        @ IntegralProjectionModelsPadrinoExt ~/Projects/ipm/IntegralProjectionModels.jl/ext/PadrinoExt/compiler.jl:609
    │      [4] _build_general_det(pm::PadrinoModel, tspan::Tuple{Int64, Int64}, mtype::IntegralProjectionModelsPadrinoExt.ModelTypeInfo)
    │        @ IntegralProjectionModelsPadrinoExt ~/Projects/ipm/IntegralProjectionModels.jl/ext/PadrinoExt/builder.jl:482
    │      [5] build_ipm(pm::PadrinoModel; tspan::Tuple{Int64, Int64})
    │        @ IntegralProjectionModelsPadrinoExt ~/Projects/ipm/IntegralProjectionModels.jl/ext/PadrinoExt/builder.jl:27
    │      [6] _pdb_make_ipm_from_protos(protos::Dict{String, PadrinoModel}; tspan::Tuple{Int64, Int64})
    │        @ IntegralProjectionModelsPadrinoExt ~/Projects/ipm/IntegralProjectionModels.jl/ext/IntegralProjectionModelsPadrinoExt.jl:337
    │      [7] _pdb_make_ipm_from_protos
    │        @ ~/Projects/ipm/IntegralProjectionModels.jl/ext/IntegralProjectionModelsPadrinoExt.jl:337 [inlined]
    │      [8] #pdb_make_ipm#69
    │        @ ~/Projects/ipm/IntegralProjectionModels.jl/ext/IntegralProjectionModelsPadrinoExt.jl:80 [inlined]
    │      [9] top-level scope
    │        @ ~/Projects/ipm/IntegralProjectionModels.jl/vignettes/10_padrino/10_padrino.qmd:95
    │     [10] eval(m::Module, e::Any)
    │        @ Core ./boot.jl:489
    │     [11] (::QuartoNotebookWorker.var"#21#22"{Module, Expr})()
    │        @ QuartoNotebookWorker ~/.julia/packages/QuartoNotebookRunner/evCNi/src/QuartoNotebookWorker/src/render.jl:222
    │     [12] (::QuartoNotebookWorker.Packages.IOCapture.var"#12#13"{Type{InterruptException}, QuartoNotebookWorker.var"#21#22"{Module, Expr}, IOContext{Base.PipeEndpoint}, IOContext{Base.PipeEndpoint}, IOContext{Base.PipeEndpoint}, IOContext{Base.PipeEndpoint}})()
    │        @ QuartoNotebookWorker.Packages.IOCapture ~/.julia/packages/QuartoNotebookRunner/evCNi/src/QuartoNotebookWorker/src/vendor/IOCapture/src/IOCapture.jl:170
    │     [13] with_logstate(f::QuartoNotebookWorker.Packages.IOCapture.var"#12#13"{Type{InterruptException}, QuartoNotebookWorker.var"#21#22"{Module, Expr}, IOContext{Base.PipeEndpoint}, IOContext{Base.PipeEndpoint}, IOContext{Base.PipeEndpoint}, IOContext{Base.PipeEndpoint}}, logstate::Base.CoreLogging.LogState)
    │        @ Base.CoreLogging ./logging/logging.jl:542
    │     [14] with_logger(f::Function, logger::Base.CoreLogging.ConsoleLogger)
    │        @ Base.CoreLogging ./logging/logging.jl:653
    │     [15] capture(f::QuartoNotebookWorker.var"#21#22"{Module, Expr}; rethrow::Type, color::Bool, passthrough::Bool, capture_buffer::IOBuffer, io_context::Vector{Pair{Symbol, Any}})
    │        @ QuartoNotebookWorker.Packages.IOCapture ~/.julia/packages/QuartoNotebookRunner/evCNi/src/QuartoNotebookWorker/src/vendor/IOCapture/src/IOCapture.jl:167
    │     [16] capture
    │        @ ~/.julia/packages/QuartoNotebookRunner/evCNi/src/QuartoNotebookWorker/src/render.jl:248 [inlined]
    │     [17] io_capture(f::Function; cell_options::Dict{String, Any}, kws::@Kwargs{rethrow::DataType, color::Bool, io_context::Vector{Pair{Symbol, Any}}})
    │        @ QuartoNotebookWorker ~/.julia/packages/QuartoNotebookRunner/evCNi/src/QuartoNotebookWorker/src/render.jl:250
    │     [18] io_capture
    │        @ ~/.julia/packages/QuartoNotebookRunner/evCNi/src/QuartoNotebookWorker/src/render.jl:246 [inlined]
    │     [19] include_str(mod::Module, code::String; file::String, line::Int64, cell_options::Dict{String, Any})
    │        @ QuartoNotebookWorker ~/.julia/packages/QuartoNotebookRunner/evCNi/src/QuartoNotebookWorker/src/render.jl:201
    │     [20] #invokelatest_gr#232
    │        @ ./reflection.jl:1297 [inlined]
    │     [21] invokelatest_gr
    │        @ ./reflection.jl:1289 [inlined]
    │     [22] #7
    │        @ ~/.julia/packages/QuartoNotebookRunner/evCNi/src/QuartoNotebookWorker/src/render.jl:18 [inlined]
    │     [23] with_inline_display(f::QuartoNotebookWorker.var"#7#8"{String, String, Int64, Dict{String, Any}}, cell_options::Dict{String, Any})
    │        @ QuartoNotebookWorker ~/.julia/packages/QuartoNotebookRunner/evCNi/src/QuartoNotebookWorker/src/InlineDisplay.jl:31
    │     [24] _render_thunk(thunk::Function, code::String, cell_options::Dict{String, Any}, is_expansion_ref::Base.RefValue{Bool}; inline::Bool)
    │        @ QuartoNotebookWorker ~/.julia/packages/QuartoNotebookRunner/evCNi/src/QuartoNotebookWorker/src/render.jl:43
    │     [25] _render_thunk
    │        @ ~/.julia/packages/QuartoNotebookRunner/evCNi/src/QuartoNotebookWorker/src/render.jl:35 [inlined]
    │     [26] #render#4
    │        @ ~/.julia/packages/QuartoNotebookRunner/evCNi/src/QuartoNotebookWorker/src/render.jl:15 [inlined]
    │     [27] render(code::String, file::String, line::Int64, cell_options::Dict{String, Any})
    │        @ QuartoNotebookWorker ~/.julia/packages/QuartoNotebookRunner/evCNi/src/QuartoNotebookWorker/src/render.jl:1
    │     [28] render(::String, ::Vararg{Any}; kwargs::@Kwargs{})
    │        @ Main ~/.julia/packages/QuartoNotebookRunner/evCNi/src/startup.jl:145
    │     [29] top-level scope
    │        @ none:1
    │     [30] eval(m::Module, e::Any)
    │        @ Core ./boot.jl:489
    │     [31] (::Main.var"#22#23")(chan::Channel{Any})
    │        @ Main ~/.julia/packages/QuartoNotebookRunner/evCNi/src/startup.jl:158
    │     [32] (::Base.var"#562#563"{Main.var"#22#23", Channel{Any}})()
    │        @ Base ./channels.jl:141
    └ @ IntegralProjectionModelsPadrinoExt ~/Projects/ipm/IntegralProjectionModels.jl/ext/IntegralProjectionModelsPadrinoExt.jl:345
    ┌ Warning: Failed to build IPM for aaaa17
    │   exception =
    │    UndefVarError: `crp_yr` not defined in `IntegralProjectionModelsPadrinoExt.PadrinoEvalModule`
    │    Suggestion: check for spelling errors or missing imports.
    │    Stacktrace:
    │      [1] (::IntegralProjectionModelsPadrinoExt.PadrinoEvalModule.var"#35#36")(z_new::Float64, z::Float64, h::Float64)
    │        @ IntegralProjectionModelsPadrinoExt.PadrinoEvalModule ~/Projects/ipm/IntegralProjectionModels.jl/ext/PadrinoExt/compiler.jl:484
    │      [2] build_kernel_matrix(kernel_fn::Function, domain::ContinuousDomain{Float64}, family::String3)
    │        @ IntegralProjectionModelsPadrinoExt ~/Projects/ipm/IntegralProjectionModels.jl/ext/PadrinoExt/compiler.jl:549
    │      [3] build_kernel_matrix(kernel_fn::Expr, vectorized_fn::Nothing, domain::ContinuousDomain{Float64}, family::String3)
    │        @ IntegralProjectionModelsPadrinoExt ~/Projects/ipm/IntegralProjectionModels.jl/ext/PadrinoExt/compiler.jl:609
    │      [4] _build_general_det(pm::PadrinoModel, tspan::Tuple{Int64, Int64}, mtype::IntegralProjectionModelsPadrinoExt.ModelTypeInfo)
    │        @ IntegralProjectionModelsPadrinoExt ~/Projects/ipm/IntegralProjectionModels.jl/ext/PadrinoExt/builder.jl:482
    │      [5] build_ipm(pm::PadrinoModel; tspan::Tuple{Int64, Int64})
    │        @ IntegralProjectionModelsPadrinoExt ~/Projects/ipm/IntegralProjectionModels.jl/ext/PadrinoExt/builder.jl:27
    │      [6] _pdb_make_ipm_from_protos(protos::Dict{String, PadrinoModel}; tspan::Tuple{Int64, Int64})
    │        @ IntegralProjectionModelsPadrinoExt ~/Projects/ipm/IntegralProjectionModels.jl/ext/IntegralProjectionModelsPadrinoExt.jl:337
    │      [7] _pdb_make_ipm_from_protos
    │        @ ~/Projects/ipm/IntegralProjectionModels.jl/ext/IntegralProjectionModelsPadrinoExt.jl:337 [inlined]
    │      [8] #pdb_make_ipm#69
    │        @ ~/Projects/ipm/IntegralProjectionModels.jl/ext/IntegralProjectionModelsPadrinoExt.jl:80 [inlined]
    │      [9] top-level scope
    │        @ ~/Projects/ipm/IntegralProjectionModels.jl/vignettes/10_padrino/10_padrino.qmd:95
    │     [10] eval(m::Module, e::Any)
    │        @ Core ./boot.jl:489
    │     [11] (::QuartoNotebookWorker.var"#21#22"{Module, Expr})()
    │        @ QuartoNotebookWorker ~/.julia/packages/QuartoNotebookRunner/evCNi/src/QuartoNotebookWorker/src/render.jl:222
    │     [12] (::QuartoNotebookWorker.Packages.IOCapture.var"#12#13"{Type{InterruptException}, QuartoNotebookWorker.var"#21#22"{Module, Expr}, IOContext{Base.PipeEndpoint}, IOContext{Base.PipeEndpoint}, IOContext{Base.PipeEndpoint}, IOContext{Base.PipeEndpoint}})()
    │        @ QuartoNotebookWorker.Packages.IOCapture ~/.julia/packages/QuartoNotebookRunner/evCNi/src/QuartoNotebookWorker/src/vendor/IOCapture/src/IOCapture.jl:170
    │     [13] with_logstate(f::QuartoNotebookWorker.Packages.IOCapture.var"#12#13"{Type{InterruptException}, QuartoNotebookWorker.var"#21#22"{Module, Expr}, IOContext{Base.PipeEndpoint}, IOContext{Base.PipeEndpoint}, IOContext{Base.PipeEndpoint}, IOContext{Base.PipeEndpoint}}, logstate::Base.CoreLogging.LogState)
    │        @ Base.CoreLogging ./logging/logging.jl:542
    │     [14] with_logger(f::Function, logger::Base.CoreLogging.ConsoleLogger)
    │        @ Base.CoreLogging ./logging/logging.jl:653
    │     [15] capture(f::QuartoNotebookWorker.var"#21#22"{Module, Expr}; rethrow::Type, color::Bool, passthrough::Bool, capture_buffer::IOBuffer, io_context::Vector{Pair{Symbol, Any}})
    │        @ QuartoNotebookWorker.Packages.IOCapture ~/.julia/packages/QuartoNotebookRunner/evCNi/src/QuartoNotebookWorker/src/vendor/IOCapture/src/IOCapture.jl:167
    │     [16] capture
    │        @ ~/.julia/packages/QuartoNotebookRunner/evCNi/src/QuartoNotebookWorker/src/render.jl:248 [inlined]
    │     [17] io_capture(f::Function; cell_options::Dict{String, Any}, kws::@Kwargs{rethrow::DataType, color::Bool, io_context::Vector{Pair{Symbol, Any}}})
    │        @ QuartoNotebookWorker ~/.julia/packages/QuartoNotebookRunner/evCNi/src/QuartoNotebookWorker/src/render.jl:250
    │     [18] io_capture
    │        @ ~/.julia/packages/QuartoNotebookRunner/evCNi/src/QuartoNotebookWorker/src/render.jl:246 [inlined]
    │     [19] include_str(mod::Module, code::String; file::String, line::Int64, cell_options::Dict{String, Any})
    │        @ QuartoNotebookWorker ~/.julia/packages/QuartoNotebookRunner/evCNi/src/QuartoNotebookWorker/src/render.jl:201
    │     [20] #invokelatest_gr#232
    │        @ ./reflection.jl:1297 [inlined]
    │     [21] invokelatest_gr
    │        @ ./reflection.jl:1289 [inlined]
    │     [22] #7
    │        @ ~/.julia/packages/QuartoNotebookRunner/evCNi/src/QuartoNotebookWorker/src/render.jl:18 [inlined]
    │     [23] with_inline_display(f::QuartoNotebookWorker.var"#7#8"{String, String, Int64, Dict{String, Any}}, cell_options::Dict{String, Any})
    │        @ QuartoNotebookWorker ~/.julia/packages/QuartoNotebookRunner/evCNi/src/QuartoNotebookWorker/src/InlineDisplay.jl:31
    │     [24] _render_thunk(thunk::Function, code::String, cell_options::Dict{String, Any}, is_expansion_ref::Base.RefValue{Bool}; inline::Bool)
    │        @ QuartoNotebookWorker ~/.julia/packages/QuartoNotebookRunner/evCNi/src/QuartoNotebookWorker/src/render.jl:43
    │     [25] _render_thunk
    │        @ ~/.julia/packages/QuartoNotebookRunner/evCNi/src/QuartoNotebookWorker/src/render.jl:35 [inlined]
    │     [26] #render#4
    │        @ ~/.julia/packages/QuartoNotebookRunner/evCNi/src/QuartoNotebookWorker/src/render.jl:15 [inlined]
    │     [27] render(code::String, file::String, line::Int64, cell_options::Dict{String, Any})
    │        @ QuartoNotebookWorker ~/.julia/packages/QuartoNotebookRunner/evCNi/src/QuartoNotebookWorker/src/render.jl:1
    │     [28] render(::String, ::Vararg{Any}; kwargs::@Kwargs{})
    │        @ Main ~/.julia/packages/QuartoNotebookRunner/evCNi/src/startup.jl:145
    │     [29] top-level scope
    │        @ none:1
    │     [30] eval(m::Module, e::Any)
    │        @ Core ./boot.jl:489
    │     [31] (::Main.var"#22#23")(chan::Channel{Any})
    │        @ Main ~/.julia/packages/QuartoNotebookRunner/evCNi/src/startup.jl:158
    │     [32] (::Base.var"#562#563"{Main.var"#22#23", Channel{Any}})()
    │        @ Base ./channels.jl:141
    └ @ IntegralProjectionModelsPadrinoExt ~/Projects/ipm/IntegralProjectionModels.jl/ext/IntegralProjectionModelsPadrinoExt.jl:345
    ┌ Warning: Failed to build env sampler for rain
    │   exception =
    │    MethodError: no method matching translate_env_function(::String31, ::Vector{String})
    │    The function `translate_env_function` exists, but no method is defined for this combination of argument types.
    │    
    │    Closest candidates are:
    │      translate_env_function(::String, ::Vector{String})
    │       @ IntegralProjectionModelsPadrinoExt ~/Projects/ipm/IntegralProjectionModels.jl/ext/PadrinoExt/translator.jl:370
    │    
    └ @ IntegralProjectionModelsPadrinoExt ~/Projects/ipm/IntegralProjectionModels.jl/ext/PadrinoExt/env_state.jl:44
    ┌ Warning: Failed to build IPM for aaaa28
    │   exception =
    │    UndefVarError: `gi_le` not defined in `IntegralProjectionModelsPadrinoExt.PadrinoEvalModule`
    │    Suggestion: check for spelling errors or missing imports.
    │    Stacktrace:
    │      [1] (::IntegralProjectionModelsPadrinoExt.PadrinoEvalModule.var"#38#39")(z_new::Float64, z::Float64, h::Float64)
    │        @ IntegralProjectionModelsPadrinoExt.PadrinoEvalModule ~/Projects/ipm/IntegralProjectionModels.jl/ext/PadrinoExt/compiler.jl:484
    │      [2] build_kernel_matrix(kernel_fn::Function, domain::ContinuousDomain{Float64}, family::String3)
    │        @ IntegralProjectionModelsPadrinoExt ~/Projects/ipm/IntegralProjectionModels.jl/ext/PadrinoExt/compiler.jl:549
    │      [3] build_kernel_matrix(kernel_fn::Expr, vectorized_fn::Nothing, domain::ContinuousDomain{Float64}, family::String3)
    │        @ IntegralProjectionModelsPadrinoExt ~/Projects/ipm/IntegralProjectionModels.jl/ext/PadrinoExt/compiler.jl:609
    │      [4] _build_general_det(pm::PadrinoModel, tspan::Tuple{Int64, Int64}, mtype::IntegralProjectionModelsPadrinoExt.ModelTypeInfo)
    │        @ IntegralProjectionModelsPadrinoExt ~/Projects/ipm/IntegralProjectionModels.jl/ext/PadrinoExt/builder.jl:482
    │      [5] build_ipm(pm::PadrinoModel; tspan::Tuple{Int64, Int64})
    │        @ IntegralProjectionModelsPadrinoExt ~/Projects/ipm/IntegralProjectionModels.jl/ext/PadrinoExt/builder.jl:27
    │      [6] _pdb_make_ipm_from_protos(protos::Dict{String, PadrinoModel}; tspan::Tuple{Int64, Int64})
    │        @ IntegralProjectionModelsPadrinoExt ~/Projects/ipm/IntegralProjectionModels.jl/ext/IntegralProjectionModelsPadrinoExt.jl:337
    │      [7] _pdb_make_ipm_from_protos
    │        @ ~/Projects/ipm/IntegralProjectionModels.jl/ext/IntegralProjectionModelsPadrinoExt.jl:337 [inlined]
    │      [8] #pdb_make_ipm#69
    │        @ ~/Projects/ipm/IntegralProjectionModels.jl/ext/IntegralProjectionModelsPadrinoExt.jl:80 [inlined]
    │      [9] top-level scope
    │        @ ~/Projects/ipm/IntegralProjectionModels.jl/vignettes/10_padrino/10_padrino.qmd:95
    │     [10] eval(m::Module, e::Any)
    │        @ Core ./boot.jl:489
    │     [11] (::QuartoNotebookWorker.var"#21#22"{Module, Expr})()
    │        @ QuartoNotebookWorker ~/.julia/packages/QuartoNotebookRunner/evCNi/src/QuartoNotebookWorker/src/render.jl:222
    │     [12] (::QuartoNotebookWorker.Packages.IOCapture.var"#12#13"{Type{InterruptException}, QuartoNotebookWorker.var"#21#22"{Module, Expr}, IOContext{Base.PipeEndpoint}, IOContext{Base.PipeEndpoint}, IOContext{Base.PipeEndpoint}, IOContext{Base.PipeEndpoint}})()
    │        @ QuartoNotebookWorker.Packages.IOCapture ~/.julia/packages/QuartoNotebookRunner/evCNi/src/QuartoNotebookWorker/src/vendor/IOCapture/src/IOCapture.jl:170
    │     [13] with_logstate(f::QuartoNotebookWorker.Packages.IOCapture.var"#12#13"{Type{InterruptException}, QuartoNotebookWorker.var"#21#22"{Module, Expr}, IOContext{Base.PipeEndpoint}, IOContext{Base.PipeEndpoint}, IOContext{Base.PipeEndpoint}, IOContext{Base.PipeEndpoint}}, logstate::Base.CoreLogging.LogState)
    │        @ Base.CoreLogging ./logging/logging.jl:542
    │     [14] with_logger(f::Function, logger::Base.CoreLogging.ConsoleLogger)
    │        @ Base.CoreLogging ./logging/logging.jl:653
    │     [15] capture(f::QuartoNotebookWorker.var"#21#22"{Module, Expr}; rethrow::Type, color::Bool, passthrough::Bool, capture_buffer::IOBuffer, io_context::Vector{Pair{Symbol, Any}})
    │        @ QuartoNotebookWorker.Packages.IOCapture ~/.julia/packages/QuartoNotebookRunner/evCNi/src/QuartoNotebookWorker/src/vendor/IOCapture/src/IOCapture.jl:167
    │     [16] capture
    │        @ ~/.julia/packages/QuartoNotebookRunner/evCNi/src/QuartoNotebookWorker/src/render.jl:248 [inlined]
    │     [17] io_capture(f::Function; cell_options::Dict{String, Any}, kws::@Kwargs{rethrow::DataType, color::Bool, io_context::Vector{Pair{Symbol, Any}}})
    │        @ QuartoNotebookWorker ~/.julia/packages/QuartoNotebookRunner/evCNi/src/QuartoNotebookWorker/src/render.jl:250
    │     [18] io_capture
    │        @ ~/.julia/packages/QuartoNotebookRunner/evCNi/src/QuartoNotebookWorker/src/render.jl:246 [inlined]
    │     [19] include_str(mod::Module, code::String; file::String, line::Int64, cell_options::Dict{String, Any})
    │        @ QuartoNotebookWorker ~/.julia/packages/QuartoNotebookRunner/evCNi/src/QuartoNotebookWorker/src/render.jl:201
    │     [20] #invokelatest_gr#232
    │        @ ./reflection.jl:1297 [inlined]
    │     [21] invokelatest_gr
    │        @ ./reflection.jl:1289 [inlined]
    │     [22] #7
    │        @ ~/.julia/packages/QuartoNotebookRunner/evCNi/src/QuartoNotebookWorker/src/render.jl:18 [inlined]
    │     [23] with_inline_display(f::QuartoNotebookWorker.var"#7#8"{String, String, Int64, Dict{String, Any}}, cell_options::Dict{String, Any})
    │        @ QuartoNotebookWorker ~/.julia/packages/QuartoNotebookRunner/evCNi/src/QuartoNotebookWorker/src/InlineDisplay.jl:31
    │     [24] _render_thunk(thunk::Function, code::String, cell_options::Dict{String, Any}, is_expansion_ref::Base.RefValue{Bool}; inline::Bool)
    │        @ QuartoNotebookWorker ~/.julia/packages/QuartoNotebookRunner/evCNi/src/QuartoNotebookWorker/src/render.jl:43
    │     [25] _render_thunk
    │        @ ~/.julia/packages/QuartoNotebookRunner/evCNi/src/QuartoNotebookWorker/src/render.jl:35 [inlined]
    │     [26] #render#4
    │        @ ~/.julia/packages/QuartoNotebookRunner/evCNi/src/QuartoNotebookWorker/src/render.jl:15 [inlined]
    │     [27] render(code::String, file::String, line::Int64, cell_options::Dict{String, Any})
    │        @ QuartoNotebookWorker ~/.julia/packages/QuartoNotebookRunner/evCNi/src/QuartoNotebookWorker/src/render.jl:1
    │     [28] render(::String, ::Vararg{Any}; kwargs::@Kwargs{})
    │        @ Main ~/.julia/packages/QuartoNotebookRunner/evCNi/src/startup.jl:145
    │     [29] top-level scope
    │        @ none:1
    │     [30] eval(m::Module, e::Any)
    │        @ Core ./boot.jl:489
    │     [31] (::Main.var"#22#23")(chan::Channel{Any})
    │        @ Main ~/.julia/packages/QuartoNotebookRunner/evCNi/src/startup.jl:158
    │     [32] (::Base.var"#562#563"{Main.var"#22#23", Channel{Any}})()
    │        @ Base ./channels.jl:141
    └ @ IntegralProjectionModelsPadrinoExt ~/Projects/ipm/IntegralProjectionModels.jl/ext/IntegralProjectionModelsPadrinoExt.jl:345
    ┌ Warning: Failed to build IPM for aaaa31
    │   exception =
    │    UndefVarError: `gi_le` not defined in `IntegralProjectionModelsPadrinoExt.PadrinoEvalModule`
    │    Suggestion: check for spelling errors or missing imports.
    │    Stacktrace:
    │      [1] (::IntegralProjectionModelsPadrinoExt.PadrinoEvalModule.var"#41#42")(z_new::Float64, z::Float64, h::Float64)
    │        @ IntegralProjectionModelsPadrinoExt.PadrinoEvalModule ~/Projects/ipm/IntegralProjectionModels.jl/ext/PadrinoExt/compiler.jl:484
    │      [2] build_kernel_matrix(kernel_fn::Function, domain::ContinuousDomain{Float64}, family::String3)
    │        @ IntegralProjectionModelsPadrinoExt ~/Projects/ipm/IntegralProjectionModels.jl/ext/PadrinoExt/compiler.jl:549
    │      [3] build_kernel_matrix(kernel_fn::Expr, vectorized_fn::Nothing, domain::ContinuousDomain{Float64}, family::String3)
    │        @ IntegralProjectionModelsPadrinoExt ~/Projects/ipm/IntegralProjectionModels.jl/ext/PadrinoExt/compiler.jl:609
    │      [4] _build_general_det(pm::PadrinoModel, tspan::Tuple{Int64, Int64}, mtype::IntegralProjectionModelsPadrinoExt.ModelTypeInfo)
    │        @ IntegralProjectionModelsPadrinoExt ~/Projects/ipm/IntegralProjectionModels.jl/ext/PadrinoExt/builder.jl:482
    │      [5] build_ipm(pm::PadrinoModel; tspan::Tuple{Int64, Int64})
    │        @ IntegralProjectionModelsPadrinoExt ~/Projects/ipm/IntegralProjectionModels.jl/ext/PadrinoExt/builder.jl:27
    │      [6] _pdb_make_ipm_from_protos(protos::Dict{String, PadrinoModel}; tspan::Tuple{Int64, Int64})
    │        @ IntegralProjectionModelsPadrinoExt ~/Projects/ipm/IntegralProjectionModels.jl/ext/IntegralProjectionModelsPadrinoExt.jl:337
    │      [7] _pdb_make_ipm_from_protos
    │        @ ~/Projects/ipm/IntegralProjectionModels.jl/ext/IntegralProjectionModelsPadrinoExt.jl:337 [inlined]
    │      [8] #pdb_make_ipm#69
    │        @ ~/Projects/ipm/IntegralProjectionModels.jl/ext/IntegralProjectionModelsPadrinoExt.jl:80 [inlined]
    │      [9] top-level scope
    │        @ ~/Projects/ipm/IntegralProjectionModels.jl/vignettes/10_padrino/10_padrino.qmd:95
    │     [10] eval(m::Module, e::Any)
    │        @ Core ./boot.jl:489
    │     [11] (::QuartoNotebookWorker.var"#21#22"{Module, Expr})()
    │        @ QuartoNotebookWorker ~/.julia/packages/QuartoNotebookRunner/evCNi/src/QuartoNotebookWorker/src/render.jl:222
    │     [12] (::QuartoNotebookWorker.Packages.IOCapture.var"#12#13"{Type{InterruptException}, QuartoNotebookWorker.var"#21#22"{Module, Expr}, IOContext{Base.PipeEndpoint}, IOContext{Base.PipeEndpoint}, IOContext{Base.PipeEndpoint}, IOContext{Base.PipeEndpoint}})()
    │        @ QuartoNotebookWorker.Packages.IOCapture ~/.julia/packages/QuartoNotebookRunner/evCNi/src/QuartoNotebookWorker/src/vendor/IOCapture/src/IOCapture.jl:170
    │     [13] with_logstate(f::QuartoNotebookWorker.Packages.IOCapture.var"#12#13"{Type{InterruptException}, QuartoNotebookWorker.var"#21#22"{Module, Expr}, IOContext{Base.PipeEndpoint}, IOContext{Base.PipeEndpoint}, IOContext{Base.PipeEndpoint}, IOContext{Base.PipeEndpoint}}, logstate::Base.CoreLogging.LogState)
    │        @ Base.CoreLogging ./logging/logging.jl:542
    │     [14] with_logger(f::Function, logger::Base.CoreLogging.ConsoleLogger)
    │        @ Base.CoreLogging ./logging/logging.jl:653
    │     [15] capture(f::QuartoNotebookWorker.var"#21#22"{Module, Expr}; rethrow::Type, color::Bool, passthrough::Bool, capture_buffer::IOBuffer, io_context::Vector{Pair{Symbol, Any}})
    │        @ QuartoNotebookWorker.Packages.IOCapture ~/.julia/packages/QuartoNotebookRunner/evCNi/src/QuartoNotebookWorker/src/vendor/IOCapture/src/IOCapture.jl:167
    │     [16] capture
    │        @ ~/.julia/packages/QuartoNotebookRunner/evCNi/src/QuartoNotebookWorker/src/render.jl:248 [inlined]
    │     [17] io_capture(f::Function; cell_options::Dict{String, Any}, kws::@Kwargs{rethrow::DataType, color::Bool, io_context::Vector{Pair{Symbol, Any}}})
    │        @ QuartoNotebookWorker ~/.julia/packages/QuartoNotebookRunner/evCNi/src/QuartoNotebookWorker/src/render.jl:250
    │     [18] io_capture
    │        @ ~/.julia/packages/QuartoNotebookRunner/evCNi/src/QuartoNotebookWorker/src/render.jl:246 [inlined]
    │     [19] include_str(mod::Module, code::String; file::String, line::Int64, cell_options::Dict{String, Any})
    │        @ QuartoNotebookWorker ~/.julia/packages/QuartoNotebookRunner/evCNi/src/QuartoNotebookWorker/src/render.jl:201
    │     [20] #invokelatest_gr#232
    │        @ ./reflection.jl:1297 [inlined]
    │     [21] invokelatest_gr
    │        @ ./reflection.jl:1289 [inlined]
    │     [22] #7
    │        @ ~/.julia/packages/QuartoNotebookRunner/evCNi/src/QuartoNotebookWorker/src/render.jl:18 [inlined]
    │     [23] with_inline_display(f::QuartoNotebookWorker.var"#7#8"{String, String, Int64, Dict{String, Any}}, cell_options::Dict{String, Any})
    │        @ QuartoNotebookWorker ~/.julia/packages/QuartoNotebookRunner/evCNi/src/QuartoNotebookWorker/src/InlineDisplay.jl:31
    │     [24] _render_thunk(thunk::Function, code::String, cell_options::Dict{String, Any}, is_expansion_ref::Base.RefValue{Bool}; inline::Bool)
    │        @ QuartoNotebookWorker ~/.julia/packages/QuartoNotebookRunner/evCNi/src/QuartoNotebookWorker/src/render.jl:43
    │     [25] _render_thunk
    │        @ ~/.julia/packages/QuartoNotebookRunner/evCNi/src/QuartoNotebookWorker/src/render.jl:35 [inlined]
    │     [26] #render#4
    │        @ ~/.julia/packages/QuartoNotebookRunner/evCNi/src/QuartoNotebookWorker/src/render.jl:15 [inlined]
    │     [27] render(code::String, file::String, line::Int64, cell_options::Dict{String, Any})
    │        @ QuartoNotebookWorker ~/.julia/packages/QuartoNotebookRunner/evCNi/src/QuartoNotebookWorker/src/render.jl:1
    │     [28] render(::String, ::Vararg{Any}; kwargs::@Kwargs{})
    │        @ Main ~/.julia/packages/QuartoNotebookRunner/evCNi/src/startup.jl:145
    │     [29] top-level scope
    │        @ none:1
    │     [30] eval(m::Module, e::Any)
    │        @ Core ./boot.jl:489
    │     [31] (::Main.var"#22#23")(chan::Channel{Any})
    │        @ Main ~/.julia/packages/QuartoNotebookRunner/evCNi/src/startup.jl:158
    │     [32] (::Base.var"#562#563"{Main.var"#22#23", Channel{Any}})()
    │        @ Base ./channels.jl:141
    └ @ IntegralProjectionModelsPadrinoExt ~/Projects/ipm/IntegralProjectionModels.jl/ext/IntegralProjectionModelsPadrinoExt.jl:345
    ┌ Warning: Failed to build IPM for aaaa20
    │   exception =
    │    UndefVarError: `crp_yr` not defined in `IntegralProjectionModelsPadrinoExt.PadrinoEvalModule`
    │    Suggestion: check for spelling errors or missing imports.
    │    Stacktrace:
    │      [1] (::IntegralProjectionModelsPadrinoExt.PadrinoEvalModule.var"#47#48")(z_new::Float64, z::Float64, h::Float64)
    │        @ IntegralProjectionModelsPadrinoExt.PadrinoEvalModule ~/Projects/ipm/IntegralProjectionModels.jl/ext/PadrinoExt/compiler.jl:484
    │      [2] build_kernel_matrix(kernel_fn::Function, domain::ContinuousDomain{Float64}, family::String3)
    │        @ IntegralProjectionModelsPadrinoExt ~/Projects/ipm/IntegralProjectionModels.jl/ext/PadrinoExt/compiler.jl:549
    │      [3] build_kernel_matrix(kernel_fn::Expr, vectorized_fn::Nothing, domain::ContinuousDomain{Float64}, family::String3)
    │        @ IntegralProjectionModelsPadrinoExt ~/Projects/ipm/IntegralProjectionModels.jl/ext/PadrinoExt/compiler.jl:609
    │      [4] _build_general_det(pm::PadrinoModel, tspan::Tuple{Int64, Int64}, mtype::IntegralProjectionModelsPadrinoExt.ModelTypeInfo)
    │        @ IntegralProjectionModelsPadrinoExt ~/Projects/ipm/IntegralProjectionModels.jl/ext/PadrinoExt/builder.jl:482
    │      [5] build_ipm(pm::PadrinoModel; tspan::Tuple{Int64, Int64})
    │        @ IntegralProjectionModelsPadrinoExt ~/Projects/ipm/IntegralProjectionModels.jl/ext/PadrinoExt/builder.jl:27
    │      [6] _pdb_make_ipm_from_protos(protos::Dict{String, PadrinoModel}; tspan::Tuple{Int64, Int64})
    │        @ IntegralProjectionModelsPadrinoExt ~/Projects/ipm/IntegralProjectionModels.jl/ext/IntegralProjectionModelsPadrinoExt.jl:337
    │      [7] _pdb_make_ipm_from_protos
    │        @ ~/Projects/ipm/IntegralProjectionModels.jl/ext/IntegralProjectionModelsPadrinoExt.jl:337 [inlined]
    │      [8] #pdb_make_ipm#69
    │        @ ~/Projects/ipm/IntegralProjectionModels.jl/ext/IntegralProjectionModelsPadrinoExt.jl:80 [inlined]
    │      [9] top-level scope
    │        @ ~/Projects/ipm/IntegralProjectionModels.jl/vignettes/10_padrino/10_padrino.qmd:95
    │     [10] eval(m::Module, e::Any)
    │        @ Core ./boot.jl:489
    │     [11] (::QuartoNotebookWorker.var"#21#22"{Module, Expr})()
    │        @ QuartoNotebookWorker ~/.julia/packages/QuartoNotebookRunner/evCNi/src/QuartoNotebookWorker/src/render.jl:222
    │     [12] (::QuartoNotebookWorker.Packages.IOCapture.var"#12#13"{Type{InterruptException}, QuartoNotebookWorker.var"#21#22"{Module, Expr}, IOContext{Base.PipeEndpoint}, IOContext{Base.PipeEndpoint}, IOContext{Base.PipeEndpoint}, IOContext{Base.PipeEndpoint}})()
    │        @ QuartoNotebookWorker.Packages.IOCapture ~/.julia/packages/QuartoNotebookRunner/evCNi/src/QuartoNotebookWorker/src/vendor/IOCapture/src/IOCapture.jl:170
    │     [13] with_logstate(f::QuartoNotebookWorker.Packages.IOCapture.var"#12#13"{Type{InterruptException}, QuartoNotebookWorker.var"#21#22"{Module, Expr}, IOContext{Base.PipeEndpoint}, IOContext{Base.PipeEndpoint}, IOContext{Base.PipeEndpoint}, IOContext{Base.PipeEndpoint}}, logstate::Base.CoreLogging.LogState)
    │        @ Base.CoreLogging ./logging/logging.jl:542
    │     [14] with_logger(f::Function, logger::Base.CoreLogging.ConsoleLogger)
    │        @ Base.CoreLogging ./logging/logging.jl:653
    │     [15] capture(f::QuartoNotebookWorker.var"#21#22"{Module, Expr}; rethrow::Type, color::Bool, passthrough::Bool, capture_buffer::IOBuffer, io_context::Vector{Pair{Symbol, Any}})
    │        @ QuartoNotebookWorker.Packages.IOCapture ~/.julia/packages/QuartoNotebookRunner/evCNi/src/QuartoNotebookWorker/src/vendor/IOCapture/src/IOCapture.jl:167
    │     [16] capture
    │        @ ~/.julia/packages/QuartoNotebookRunner/evCNi/src/QuartoNotebookWorker/src/render.jl:248 [inlined]
    │     [17] io_capture(f::Function; cell_options::Dict{String, Any}, kws::@Kwargs{rethrow::DataType, color::Bool, io_context::Vector{Pair{Symbol, Any}}})
    │        @ QuartoNotebookWorker ~/.julia/packages/QuartoNotebookRunner/evCNi/src/QuartoNotebookWorker/src/render.jl:250
    │     [18] io_capture
    │        @ ~/.julia/packages/QuartoNotebookRunner/evCNi/src/QuartoNotebookWorker/src/render.jl:246 [inlined]
    │     [19] include_str(mod::Module, code::String; file::String, line::Int64, cell_options::Dict{String, Any})
    │        @ QuartoNotebookWorker ~/.julia/packages/QuartoNotebookRunner/evCNi/src/QuartoNotebookWorker/src/render.jl:201
    │     [20] #invokelatest_gr#232
    │        @ ./reflection.jl:1297 [inlined]
    │     [21] invokelatest_gr
    │        @ ./reflection.jl:1289 [inlined]
    │     [22] #7
    │        @ ~/.julia/packages/QuartoNotebookRunner/evCNi/src/QuartoNotebookWorker/src/render.jl:18 [inlined]
    │     [23] with_inline_display(f::QuartoNotebookWorker.var"#7#8"{String, String, Int64, Dict{String, Any}}, cell_options::Dict{String, Any})
    │        @ QuartoNotebookWorker ~/.julia/packages/QuartoNotebookRunner/evCNi/src/QuartoNotebookWorker/src/InlineDisplay.jl:31
    │     [24] _render_thunk(thunk::Function, code::String, cell_options::Dict{String, Any}, is_expansion_ref::Base.RefValue{Bool}; inline::Bool)
    │        @ QuartoNotebookWorker ~/.julia/packages/QuartoNotebookRunner/evCNi/src/QuartoNotebookWorker/src/render.jl:43
    │     [25] _render_thunk
    │        @ ~/.julia/packages/QuartoNotebookRunner/evCNi/src/QuartoNotebookWorker/src/render.jl:35 [inlined]
    │     [26] #render#4
    │        @ ~/.julia/packages/QuartoNotebookRunner/evCNi/src/QuartoNotebookWorker/src/render.jl:15 [inlined]
    │     [27] render(code::String, file::String, line::Int64, cell_options::Dict{String, Any})
    │        @ QuartoNotebookWorker ~/.julia/packages/QuartoNotebookRunner/evCNi/src/QuartoNotebookWorker/src/render.jl:1
    │     [28] render(::String, ::Vararg{Any}; kwargs::@Kwargs{})
    │        @ Main ~/.julia/packages/QuartoNotebookRunner/evCNi/src/startup.jl:145
    │     [29] top-level scope
    │        @ none:1
    │     [30] eval(m::Module, e::Any)
    │        @ Core ./boot.jl:489
    │     [31] (::Main.var"#22#23")(chan::Channel{Any})
    │        @ Main ~/.julia/packages/QuartoNotebookRunner/evCNi/src/startup.jl:158
    │     [32] (::Base.var"#562#563"{Main.var"#22#23", Channel{Any}})()
    │        @ Base ./channels.jl:141
    └ @ IntegralProjectionModelsPadrinoExt ~/Projects/ipm/IntegralProjectionModels.jl/ext/IntegralProjectionModelsPadrinoExt.jl:345
    ┌ Warning: Failed to build IPM for aaaa23
    │   exception =
    │    TypeError: in typeassert, expected Float64, got a value of type Vector{Vector{Float64}}
    │    Stacktrace:
    │      [1] build_kernel_matrix(kernel_fn::Function, domain::ContinuousDomain{Float64}, family::String3)
    │        @ IntegralProjectionModelsPadrinoExt ~/Projects/ipm/IntegralProjectionModels.jl/ext/PadrinoExt/compiler.jl:549
    │      [2] build_kernel_matrix(kernel_fn::Expr, vectorized_fn::Nothing, domain::ContinuousDomain{Float64}, family::String3)
    │        @ IntegralProjectionModelsPadrinoExt ~/Projects/ipm/IntegralProjectionModels.jl/ext/PadrinoExt/compiler.jl:609
    │      [3] _build_det_par_sets(pm::PadrinoModel, tspan::Tuple{Int64, Int64}, mtype::IntegralProjectionModelsPadrinoExt.ModelTypeInfo)
    │        @ IntegralProjectionModelsPadrinoExt ~/Projects/ipm/IntegralProjectionModels.jl/ext/PadrinoExt/builder.jl:272
    │      [4] build_ipm(pm::PadrinoModel; tspan::Tuple{Int64, Int64})
    │        @ IntegralProjectionModelsPadrinoExt ~/Projects/ipm/IntegralProjectionModels.jl/ext/PadrinoExt/builder.jl:17
    │      [5] _pdb_make_ipm_from_protos(protos::Dict{String, PadrinoModel}; tspan::Tuple{Int64, Int64})
    │        @ IntegralProjectionModelsPadrinoExt ~/Projects/ipm/IntegralProjectionModels.jl/ext/IntegralProjectionModelsPadrinoExt.jl:337
    │      [6] _pdb_make_ipm_from_protos
    │        @ ~/Projects/ipm/IntegralProjectionModels.jl/ext/IntegralProjectionModelsPadrinoExt.jl:337 [inlined]
    │      [7] #pdb_make_ipm#69
    │        @ ~/Projects/ipm/IntegralProjectionModels.jl/ext/IntegralProjectionModelsPadrinoExt.jl:80 [inlined]
    │      [8] top-level scope
    │        @ ~/Projects/ipm/IntegralProjectionModels.jl/vignettes/10_padrino/10_padrino.qmd:95
    │      [9] eval(m::Module, e::Any)
    │        @ Core ./boot.jl:489
    │     [10] (::QuartoNotebookWorker.var"#21#22"{Module, Expr})()
    │        @ QuartoNotebookWorker ~/.julia/packages/QuartoNotebookRunner/evCNi/src/QuartoNotebookWorker/src/render.jl:222
    │     [11] (::QuartoNotebookWorker.Packages.IOCapture.var"#12#13"{Type{InterruptException}, QuartoNotebookWorker.var"#21#22"{Module, Expr}, IOContext{Base.PipeEndpoint}, IOContext{Base.PipeEndpoint}, IOContext{Base.PipeEndpoint}, IOContext{Base.PipeEndpoint}})()
    │        @ QuartoNotebookWorker.Packages.IOCapture ~/.julia/packages/QuartoNotebookRunner/evCNi/src/QuartoNotebookWorker/src/vendor/IOCapture/src/IOCapture.jl:170
    │     [12] with_logstate(f::QuartoNotebookWorker.Packages.IOCapture.var"#12#13"{Type{InterruptException}, QuartoNotebookWorker.var"#21#22"{Module, Expr}, IOContext{Base.PipeEndpoint}, IOContext{Base.PipeEndpoint}, IOContext{Base.PipeEndpoint}, IOContext{Base.PipeEndpoint}}, logstate::Base.CoreLogging.LogState)
    │        @ Base.CoreLogging ./logging/logging.jl:542
    │     [13] with_logger(f::Function, logger::Base.CoreLogging.ConsoleLogger)
    │        @ Base.CoreLogging ./logging/logging.jl:653
    │     [14] capture(f::QuartoNotebookWorker.var"#21#22"{Module, Expr}; rethrow::Type, color::Bool, passthrough::Bool, capture_buffer::IOBuffer, io_context::Vector{Pair{Symbol, Any}})
    │        @ QuartoNotebookWorker.Packages.IOCapture ~/.julia/packages/QuartoNotebookRunner/evCNi/src/QuartoNotebookWorker/src/vendor/IOCapture/src/IOCapture.jl:167
    │     [15] capture
    │        @ ~/.julia/packages/QuartoNotebookRunner/evCNi/src/QuartoNotebookWorker/src/render.jl:248 [inlined]
    │     [16] io_capture(f::Function; cell_options::Dict{String, Any}, kws::@Kwargs{rethrow::DataType, color::Bool, io_context::Vector{Pair{Symbol, Any}}})
    │        @ QuartoNotebookWorker ~/.julia/packages/QuartoNotebookRunner/evCNi/src/QuartoNotebookWorker/src/render.jl:250
    │     [17] io_capture
    │        @ ~/.julia/packages/QuartoNotebookRunner/evCNi/src/QuartoNotebookWorker/src/render.jl:246 [inlined]
    │     [18] include_str(mod::Module, code::String; file::String, line::Int64, cell_options::Dict{String, Any})
    │        @ QuartoNotebookWorker ~/.julia/packages/QuartoNotebookRunner/evCNi/src/QuartoNotebookWorker/src/render.jl:201
    │     [19] #invokelatest_gr#232
    │        @ ./reflection.jl:1297 [inlined]
    │     [20] invokelatest_gr
    │        @ ./reflection.jl:1289 [inlined]
    │     [21] #7
    │        @ ~/.julia/packages/QuartoNotebookRunner/evCNi/src/QuartoNotebookWorker/src/render.jl:18 [inlined]
    │     [22] with_inline_display(f::QuartoNotebookWorker.var"#7#8"{String, String, Int64, Dict{String, Any}}, cell_options::Dict{String, Any})
    │        @ QuartoNotebookWorker ~/.julia/packages/QuartoNotebookRunner/evCNi/src/QuartoNotebookWorker/src/InlineDisplay.jl:31
    │     [23] _render_thunk(thunk::Function, code::String, cell_options::Dict{String, Any}, is_expansion_ref::Base.RefValue{Bool}; inline::Bool)
    │        @ QuartoNotebookWorker ~/.julia/packages/QuartoNotebookRunner/evCNi/src/QuartoNotebookWorker/src/render.jl:43
    │     [24] _render_thunk
    │        @ ~/.julia/packages/QuartoNotebookRunner/evCNi/src/QuartoNotebookWorker/src/render.jl:35 [inlined]
    │     [25] #render#4
    │        @ ~/.julia/packages/QuartoNotebookRunner/evCNi/src/QuartoNotebookWorker/src/render.jl:15 [inlined]
    │     [26] render(code::String, file::String, line::Int64, cell_options::Dict{String, Any})
    │        @ QuartoNotebookWorker ~/.julia/packages/QuartoNotebookRunner/evCNi/src/QuartoNotebookWorker/src/render.jl:1
    │     [27] render(::String, ::Vararg{Any}; kwargs::@Kwargs{})
    │        @ Main ~/.julia/packages/QuartoNotebookRunner/evCNi/src/startup.jl:145
    │     [28] top-level scope
    │        @ none:1
    │     [29] eval(m::Module, e::Any)
    │        @ Core ./boot.jl:489
    │     [30] (::Main.var"#22#23")(chan::Channel{Any})
    │        @ Main ~/.julia/packages/QuartoNotebookRunner/evCNi/src/startup.jl:158
    │     [31] (::Base.var"#562#563"{Main.var"#22#23", Channel{Any}})()
    │        @ Base ./channels.jl:141
    └ @ IntegralProjectionModelsPadrinoExt ~/Projects/ipm/IntegralProjectionModels.jl/ext/IntegralProjectionModelsPadrinoExt.jl:345
    ┌ Warning: Failed to build IPM for aaaa25
    │   exception =
    │    UndefVarError: `gi_le` not defined in `IntegralProjectionModelsPadrinoExt.PadrinoEvalModule`
    │    Suggestion: check for spelling errors or missing imports.
    │    Stacktrace:
    │      [1] (::IntegralProjectionModelsPadrinoExt.PadrinoEvalModule.var"#53#54")(z_new::Float64, z::Float64, h::Float64)
    │        @ IntegralProjectionModelsPadrinoExt.PadrinoEvalModule ~/Projects/ipm/IntegralProjectionModels.jl/ext/PadrinoExt/compiler.jl:484
    │      [2] build_kernel_matrix(kernel_fn::Function, domain::ContinuousDomain{Float64}, family::String3)
    │        @ IntegralProjectionModelsPadrinoExt ~/Projects/ipm/IntegralProjectionModels.jl/ext/PadrinoExt/compiler.jl:549
    │      [3] build_kernel_matrix(kernel_fn::Expr, vectorized_fn::Nothing, domain::ContinuousDomain{Float64}, family::String3)
    │        @ IntegralProjectionModelsPadrinoExt ~/Projects/ipm/IntegralProjectionModels.jl/ext/PadrinoExt/compiler.jl:609
    │      [4] _build_general_det(pm::PadrinoModel, tspan::Tuple{Int64, Int64}, mtype::IntegralProjectionModelsPadrinoExt.ModelTypeInfo)
    │        @ IntegralProjectionModelsPadrinoExt ~/Projects/ipm/IntegralProjectionModels.jl/ext/PadrinoExt/builder.jl:482
    │      [5] build_ipm(pm::PadrinoModel; tspan::Tuple{Int64, Int64})
    │        @ IntegralProjectionModelsPadrinoExt ~/Projects/ipm/IntegralProjectionModels.jl/ext/PadrinoExt/builder.jl:27
    │      [6] _pdb_make_ipm_from_protos(protos::Dict{String, PadrinoModel}; tspan::Tuple{Int64, Int64})
    │        @ IntegralProjectionModelsPadrinoExt ~/Projects/ipm/IntegralProjectionModels.jl/ext/IntegralProjectionModelsPadrinoExt.jl:337
    │      [7] _pdb_make_ipm_from_protos
    │        @ ~/Projects/ipm/IntegralProjectionModels.jl/ext/IntegralProjectionModelsPadrinoExt.jl:337 [inlined]
    │      [8] #pdb_make_ipm#69
    │        @ ~/Projects/ipm/IntegralProjectionModels.jl/ext/IntegralProjectionModelsPadrinoExt.jl:80 [inlined]
    │      [9] top-level scope
    │        @ ~/Projects/ipm/IntegralProjectionModels.jl/vignettes/10_padrino/10_padrino.qmd:95
    │     [10] eval(m::Module, e::Any)
    │        @ Core ./boot.jl:489
    │     [11] (::QuartoNotebookWorker.var"#21#22"{Module, Expr})()
    │        @ QuartoNotebookWorker ~/.julia/packages/QuartoNotebookRunner/evCNi/src/QuartoNotebookWorker/src/render.jl:222
    │     [12] (::QuartoNotebookWorker.Packages.IOCapture.var"#12#13"{Type{InterruptException}, QuartoNotebookWorker.var"#21#22"{Module, Expr}, IOContext{Base.PipeEndpoint}, IOContext{Base.PipeEndpoint}, IOContext{Base.PipeEndpoint}, IOContext{Base.PipeEndpoint}})()
    │        @ QuartoNotebookWorker.Packages.IOCapture ~/.julia/packages/QuartoNotebookRunner/evCNi/src/QuartoNotebookWorker/src/vendor/IOCapture/src/IOCapture.jl:170
    │     [13] with_logstate(f::QuartoNotebookWorker.Packages.IOCapture.var"#12#13"{Type{InterruptException}, QuartoNotebookWorker.var"#21#22"{Module, Expr}, IOContext{Base.PipeEndpoint}, IOContext{Base.PipeEndpoint}, IOContext{Base.PipeEndpoint}, IOContext{Base.PipeEndpoint}}, logstate::Base.CoreLogging.LogState)
    │        @ Base.CoreLogging ./logging/logging.jl:542
    │     [14] with_logger(f::Function, logger::Base.CoreLogging.ConsoleLogger)
    │        @ Base.CoreLogging ./logging/logging.jl:653
    │     [15] capture(f::QuartoNotebookWorker.var"#21#22"{Module, Expr}; rethrow::Type, color::Bool, passthrough::Bool, capture_buffer::IOBuffer, io_context::Vector{Pair{Symbol, Any}})
    │        @ QuartoNotebookWorker.Packages.IOCapture ~/.julia/packages/QuartoNotebookRunner/evCNi/src/QuartoNotebookWorker/src/vendor/IOCapture/src/IOCapture.jl:167
    │     [16] capture
    │        @ ~/.julia/packages/QuartoNotebookRunner/evCNi/src/QuartoNotebookWorker/src/render.jl:248 [inlined]
    │     [17] io_capture(f::Function; cell_options::Dict{String, Any}, kws::@Kwargs{rethrow::DataType, color::Bool, io_context::Vector{Pair{Symbol, Any}}})
    │        @ QuartoNotebookWorker ~/.julia/packages/QuartoNotebookRunner/evCNi/src/QuartoNotebookWorker/src/render.jl:250
    │     [18] io_capture
    │        @ ~/.julia/packages/QuartoNotebookRunner/evCNi/src/QuartoNotebookWorker/src/render.jl:246 [inlined]
    │     [19] include_str(mod::Module, code::String; file::String, line::Int64, cell_options::Dict{String, Any})
    │        @ QuartoNotebookWorker ~/.julia/packages/QuartoNotebookRunner/evCNi/src/QuartoNotebookWorker/src/render.jl:201
    │     [20] #invokelatest_gr#232
    │        @ ./reflection.jl:1297 [inlined]
    │     [21] invokelatest_gr
    │        @ ./reflection.jl:1289 [inlined]
    │     [22] #7
    │        @ ~/.julia/packages/QuartoNotebookRunner/evCNi/src/QuartoNotebookWorker/src/render.jl:18 [inlined]
    │     [23] with_inline_display(f::QuartoNotebookWorker.var"#7#8"{String, String, Int64, Dict{String, Any}}, cell_options::Dict{String, Any})
    │        @ QuartoNotebookWorker ~/.julia/packages/QuartoNotebookRunner/evCNi/src/QuartoNotebookWorker/src/InlineDisplay.jl:31
    │     [24] _render_thunk(thunk::Function, code::String, cell_options::Dict{String, Any}, is_expansion_ref::Base.RefValue{Bool}; inline::Bool)
    │        @ QuartoNotebookWorker ~/.julia/packages/QuartoNotebookRunner/evCNi/src/QuartoNotebookWorker/src/render.jl:43
    │     [25] _render_thunk
    │        @ ~/.julia/packages/QuartoNotebookRunner/evCNi/src/QuartoNotebookWorker/src/render.jl:35 [inlined]
    │     [26] #render#4
    │        @ ~/.julia/packages/QuartoNotebookRunner/evCNi/src/QuartoNotebookWorker/src/render.jl:15 [inlined]
    │     [27] render(code::String, file::String, line::Int64, cell_options::Dict{String, Any})
    │        @ QuartoNotebookWorker ~/.julia/packages/QuartoNotebookRunner/evCNi/src/QuartoNotebookWorker/src/render.jl:1
    │     [28] render(::String, ::Vararg{Any}; kwargs::@Kwargs{})
    │        @ Main ~/.julia/packages/QuartoNotebookRunner/evCNi/src/startup.jl:145
    │     [29] top-level scope
    │        @ none:1
    │     [30] eval(m::Module, e::Any)
    │        @ Core ./boot.jl:489
    │     [31] (::Main.var"#22#23")(chan::Channel{Any})
    │        @ Main ~/.julia/packages/QuartoNotebookRunner/evCNi/src/startup.jl:158
    │     [32] (::Base.var"#562#563"{Main.var"#22#23", Channel{Any}})()
    │        @ Base ./channels.jl:141
    └ @ IntegralProjectionModelsPadrinoExt ~/Projects/ipm/IntegralProjectionModels.jl/ext/IntegralProjectionModelsPadrinoExt.jl:345
    ┌ Warning: Failed to build IPM for aaaa29
    │   exception =
    │    UndefVarError: `gi_le` not defined in `IntegralProjectionModelsPadrinoExt.PadrinoEvalModule`
    │    Suggestion: check for spelling errors or missing imports.
    │    Stacktrace:
    │      [1] (::IntegralProjectionModelsPadrinoExt.PadrinoEvalModule.var"#56#57")(z_new::Float64, z::Float64, h::Float64)
    │        @ IntegralProjectionModelsPadrinoExt.PadrinoEvalModule ~/Projects/ipm/IntegralProjectionModels.jl/ext/PadrinoExt/compiler.jl:484
    │      [2] build_kernel_matrix(kernel_fn::Function, domain::ContinuousDomain{Float64}, family::String3)
    │        @ IntegralProjectionModelsPadrinoExt ~/Projects/ipm/IntegralProjectionModels.jl/ext/PadrinoExt/compiler.jl:549
    │      [3] build_kernel_matrix(kernel_fn::Expr, vectorized_fn::Nothing, domain::ContinuousDomain{Float64}, family::String3)
    │        @ IntegralProjectionModelsPadrinoExt ~/Projects/ipm/IntegralProjectionModels.jl/ext/PadrinoExt/compiler.jl:609
    │      [4] _build_general_det(pm::PadrinoModel, tspan::Tuple{Int64, Int64}, mtype::IntegralProjectionModelsPadrinoExt.ModelTypeInfo)
    │        @ IntegralProjectionModelsPadrinoExt ~/Projects/ipm/IntegralProjectionModels.jl/ext/PadrinoExt/builder.jl:482
    │      [5] build_ipm(pm::PadrinoModel; tspan::Tuple{Int64, Int64})
    │        @ IntegralProjectionModelsPadrinoExt ~/Projects/ipm/IntegralProjectionModels.jl/ext/PadrinoExt/builder.jl:27
    │      [6] _pdb_make_ipm_from_protos(protos::Dict{String, PadrinoModel}; tspan::Tuple{Int64, Int64})
    │        @ IntegralProjectionModelsPadrinoExt ~/Projects/ipm/IntegralProjectionModels.jl/ext/IntegralProjectionModelsPadrinoExt.jl:337
    │      [7] _pdb_make_ipm_from_protos
    │        @ ~/Projects/ipm/IntegralProjectionModels.jl/ext/IntegralProjectionModelsPadrinoExt.jl:337 [inlined]
    │      [8] #pdb_make_ipm#69
    │        @ ~/Projects/ipm/IntegralProjectionModels.jl/ext/IntegralProjectionModelsPadrinoExt.jl:80 [inlined]
    │      [9] top-level scope
    │        @ ~/Projects/ipm/IntegralProjectionModels.jl/vignettes/10_padrino/10_padrino.qmd:95
    │     [10] eval(m::Module, e::Any)
    │        @ Core ./boot.jl:489
    │     [11] (::QuartoNotebookWorker.var"#21#22"{Module, Expr})()
    │        @ QuartoNotebookWorker ~/.julia/packages/QuartoNotebookRunner/evCNi/src/QuartoNotebookWorker/src/render.jl:222
    │     [12] (::QuartoNotebookWorker.Packages.IOCapture.var"#12#13"{Type{InterruptException}, QuartoNotebookWorker.var"#21#22"{Module, Expr}, IOContext{Base.PipeEndpoint}, IOContext{Base.PipeEndpoint}, IOContext{Base.PipeEndpoint}, IOContext{Base.PipeEndpoint}})()
    │        @ QuartoNotebookWorker.Packages.IOCapture ~/.julia/packages/QuartoNotebookRunner/evCNi/src/QuartoNotebookWorker/src/vendor/IOCapture/src/IOCapture.jl:170
    │     [13] with_logstate(f::QuartoNotebookWorker.Packages.IOCapture.var"#12#13"{Type{InterruptException}, QuartoNotebookWorker.var"#21#22"{Module, Expr}, IOContext{Base.PipeEndpoint}, IOContext{Base.PipeEndpoint}, IOContext{Base.PipeEndpoint}, IOContext{Base.PipeEndpoint}}, logstate::Base.CoreLogging.LogState)
    │        @ Base.CoreLogging ./logging/logging.jl:542
    │     [14] with_logger(f::Function, logger::Base.CoreLogging.ConsoleLogger)
    │        @ Base.CoreLogging ./logging/logging.jl:653
    │     [15] capture(f::QuartoNotebookWorker.var"#21#22"{Module, Expr}; rethrow::Type, color::Bool, passthrough::Bool, capture_buffer::IOBuffer, io_context::Vector{Pair{Symbol, Any}})
    │        @ QuartoNotebookWorker.Packages.IOCapture ~/.julia/packages/QuartoNotebookRunner/evCNi/src/QuartoNotebookWorker/src/vendor/IOCapture/src/IOCapture.jl:167
    │     [16] capture
    │        @ ~/.julia/packages/QuartoNotebookRunner/evCNi/src/QuartoNotebookWorker/src/render.jl:248 [inlined]
    │     [17] io_capture(f::Function; cell_options::Dict{String, Any}, kws::@Kwargs{rethrow::DataType, color::Bool, io_context::Vector{Pair{Symbol, Any}}})
    │        @ QuartoNotebookWorker ~/.julia/packages/QuartoNotebookRunner/evCNi/src/QuartoNotebookWorker/src/render.jl:250
    │     [18] io_capture
    │        @ ~/.julia/packages/QuartoNotebookRunner/evCNi/src/QuartoNotebookWorker/src/render.jl:246 [inlined]
    │     [19] include_str(mod::Module, code::String; file::String, line::Int64, cell_options::Dict{String, Any})
    │        @ QuartoNotebookWorker ~/.julia/packages/QuartoNotebookRunner/evCNi/src/QuartoNotebookWorker/src/render.jl:201
    │     [20] #invokelatest_gr#232
    │        @ ./reflection.jl:1297 [inlined]
    │     [21] invokelatest_gr
    │        @ ./reflection.jl:1289 [inlined]
    │     [22] #7
    │        @ ~/.julia/packages/QuartoNotebookRunner/evCNi/src/QuartoNotebookWorker/src/render.jl:18 [inlined]
    │     [23] with_inline_display(f::QuartoNotebookWorker.var"#7#8"{String, String, Int64, Dict{String, Any}}, cell_options::Dict{String, Any})
    │        @ QuartoNotebookWorker ~/.julia/packages/QuartoNotebookRunner/evCNi/src/QuartoNotebookWorker/src/InlineDisplay.jl:31
    │     [24] _render_thunk(thunk::Function, code::String, cell_options::Dict{String, Any}, is_expansion_ref::Base.RefValue{Bool}; inline::Bool)
    │        @ QuartoNotebookWorker ~/.julia/packages/QuartoNotebookRunner/evCNi/src/QuartoNotebookWorker/src/render.jl:43
    │     [25] _render_thunk
    │        @ ~/.julia/packages/QuartoNotebookRunner/evCNi/src/QuartoNotebookWorker/src/render.jl:35 [inlined]
    │     [26] #render#4
    │        @ ~/.julia/packages/QuartoNotebookRunner/evCNi/src/QuartoNotebookWorker/src/render.jl:15 [inlined]
    │     [27] render(code::String, file::String, line::Int64, cell_options::Dict{String, Any})
    │        @ QuartoNotebookWorker ~/.julia/packages/QuartoNotebookRunner/evCNi/src/QuartoNotebookWorker/src/render.jl:1
    │     [28] render(::String, ::Vararg{Any}; kwargs::@Kwargs{})
    │        @ Main ~/.julia/packages/QuartoNotebookRunner/evCNi/src/startup.jl:145
    │     [29] top-level scope
    │        @ none:1
    │     [30] eval(m::Module, e::Any)
    │        @ Core ./boot.jl:489
    │     [31] (::Main.var"#22#23")(chan::Channel{Any})
    │        @ Main ~/.julia/packages/QuartoNotebookRunner/evCNi/src/startup.jl:158
    │     [32] (::Base.var"#562#563"{Main.var"#22#23", Channel{Any}})()
    │        @ Base ./channels.jl:141
    └ @ IntegralProjectionModelsPadrinoExt ~/Projects/ipm/IntegralProjectionModels.jl/ext/IntegralProjectionModelsPadrinoExt.jl:345
    ┌ Warning: Failed to build IPM for aaaa18
    │   exception =
    │    UndefVarError: `crp_yr` not defined in `IntegralProjectionModelsPadrinoExt.PadrinoEvalModule`
    │    Suggestion: check for spelling errors or missing imports.
    │    Stacktrace:
    │      [1] (::IntegralProjectionModelsPadrinoExt.PadrinoEvalModule.var"#62#63")(z_new::Float64, z::Float64, h::Float64)
    │        @ IntegralProjectionModelsPadrinoExt.PadrinoEvalModule ~/Projects/ipm/IntegralProjectionModels.jl/ext/PadrinoExt/compiler.jl:484
    │      [2] build_kernel_matrix(kernel_fn::Function, domain::ContinuousDomain{Float64}, family::String3)
    │        @ IntegralProjectionModelsPadrinoExt ~/Projects/ipm/IntegralProjectionModels.jl/ext/PadrinoExt/compiler.jl:549
    │      [3] build_kernel_matrix(kernel_fn::Expr, vectorized_fn::Nothing, domain::ContinuousDomain{Float64}, family::String3)
    │        @ IntegralProjectionModelsPadrinoExt ~/Projects/ipm/IntegralProjectionModels.jl/ext/PadrinoExt/compiler.jl:609
    │      [4] _build_general_det(pm::PadrinoModel, tspan::Tuple{Int64, Int64}, mtype::IntegralProjectionModelsPadrinoExt.ModelTypeInfo)
    │        @ IntegralProjectionModelsPadrinoExt ~/Projects/ipm/IntegralProjectionModels.jl/ext/PadrinoExt/builder.jl:482
    │      [5] build_ipm(pm::PadrinoModel; tspan::Tuple{Int64, Int64})
    │        @ IntegralProjectionModelsPadrinoExt ~/Projects/ipm/IntegralProjectionModels.jl/ext/PadrinoExt/builder.jl:27
    │      [6] _pdb_make_ipm_from_protos(protos::Dict{String, PadrinoModel}; tspan::Tuple{Int64, Int64})
    │        @ IntegralProjectionModelsPadrinoExt ~/Projects/ipm/IntegralProjectionModels.jl/ext/IntegralProjectionModelsPadrinoExt.jl:337
    │      [7] _pdb_make_ipm_from_protos
    │        @ ~/Projects/ipm/IntegralProjectionModels.jl/ext/IntegralProjectionModelsPadrinoExt.jl:337 [inlined]
    │      [8] #pdb_make_ipm#69
    │        @ ~/Projects/ipm/IntegralProjectionModels.jl/ext/IntegralProjectionModelsPadrinoExt.jl:80 [inlined]
    │      [9] top-level scope
    │        @ ~/Projects/ipm/IntegralProjectionModels.jl/vignettes/10_padrino/10_padrino.qmd:95
    │     [10] eval(m::Module, e::Any)
    │        @ Core ./boot.jl:489
    │     [11] (::QuartoNotebookWorker.var"#21#22"{Module, Expr})()
    │        @ QuartoNotebookWorker ~/.julia/packages/QuartoNotebookRunner/evCNi/src/QuartoNotebookWorker/src/render.jl:222
    │     [12] (::QuartoNotebookWorker.Packages.IOCapture.var"#12#13"{Type{InterruptException}, QuartoNotebookWorker.var"#21#22"{Module, Expr}, IOContext{Base.PipeEndpoint}, IOContext{Base.PipeEndpoint}, IOContext{Base.PipeEndpoint}, IOContext{Base.PipeEndpoint}})()
    │        @ QuartoNotebookWorker.Packages.IOCapture ~/.julia/packages/QuartoNotebookRunner/evCNi/src/QuartoNotebookWorker/src/vendor/IOCapture/src/IOCapture.jl:170
    │     [13] with_logstate(f::QuartoNotebookWorker.Packages.IOCapture.var"#12#13"{Type{InterruptException}, QuartoNotebookWorker.var"#21#22"{Module, Expr}, IOContext{Base.PipeEndpoint}, IOContext{Base.PipeEndpoint}, IOContext{Base.PipeEndpoint}, IOContext{Base.PipeEndpoint}}, logstate::Base.CoreLogging.LogState)
    │        @ Base.CoreLogging ./logging/logging.jl:542
    │     [14] with_logger(f::Function, logger::Base.CoreLogging.ConsoleLogger)
    │        @ Base.CoreLogging ./logging/logging.jl:653
    │     [15] capture(f::QuartoNotebookWorker.var"#21#22"{Module, Expr}; rethrow::Type, color::Bool, passthrough::Bool, capture_buffer::IOBuffer, io_context::Vector{Pair{Symbol, Any}})
    │        @ QuartoNotebookWorker.Packages.IOCapture ~/.julia/packages/QuartoNotebookRunner/evCNi/src/QuartoNotebookWorker/src/vendor/IOCapture/src/IOCapture.jl:167
    │     [16] capture
    │        @ ~/.julia/packages/QuartoNotebookRunner/evCNi/src/QuartoNotebookWorker/src/render.jl:248 [inlined]
    │     [17] io_capture(f::Function; cell_options::Dict{String, Any}, kws::@Kwargs{rethrow::DataType, color::Bool, io_context::Vector{Pair{Symbol, Any}}})
    │        @ QuartoNotebookWorker ~/.julia/packages/QuartoNotebookRunner/evCNi/src/QuartoNotebookWorker/src/render.jl:250
    │     [18] io_capture
    │        @ ~/.julia/packages/QuartoNotebookRunner/evCNi/src/QuartoNotebookWorker/src/render.jl:246 [inlined]
    │     [19] include_str(mod::Module, code::String; file::String, line::Int64, cell_options::Dict{String, Any})
    │        @ QuartoNotebookWorker ~/.julia/packages/QuartoNotebookRunner/evCNi/src/QuartoNotebookWorker/src/render.jl:201
    │     [20] #invokelatest_gr#232
    │        @ ./reflection.jl:1297 [inlined]
    │     [21] invokelatest_gr
    │        @ ./reflection.jl:1289 [inlined]
    │     [22] #7
    │        @ ~/.julia/packages/QuartoNotebookRunner/evCNi/src/QuartoNotebookWorker/src/render.jl:18 [inlined]
    │     [23] with_inline_display(f::QuartoNotebookWorker.var"#7#8"{String, String, Int64, Dict{String, Any}}, cell_options::Dict{String, Any})
    │        @ QuartoNotebookWorker ~/.julia/packages/QuartoNotebookRunner/evCNi/src/QuartoNotebookWorker/src/InlineDisplay.jl:31
    │     [24] _render_thunk(thunk::Function, code::String, cell_options::Dict{String, Any}, is_expansion_ref::Base.RefValue{Bool}; inline::Bool)
    │        @ QuartoNotebookWorker ~/.julia/packages/QuartoNotebookRunner/evCNi/src/QuartoNotebookWorker/src/render.jl:43
    │     [25] _render_thunk
    │        @ ~/.julia/packages/QuartoNotebookRunner/evCNi/src/QuartoNotebookWorker/src/render.jl:35 [inlined]
    │     [26] #render#4
    │        @ ~/.julia/packages/QuartoNotebookRunner/evCNi/src/QuartoNotebookWorker/src/render.jl:15 [inlined]
    │     [27] render(code::String, file::String, line::Int64, cell_options::Dict{String, Any})
    │        @ QuartoNotebookWorker ~/.julia/packages/QuartoNotebookRunner/evCNi/src/QuartoNotebookWorker/src/render.jl:1
    │     [28] render(::String, ::Vararg{Any}; kwargs::@Kwargs{})
    │        @ Main ~/.julia/packages/QuartoNotebookRunner/evCNi/src/startup.jl:145
    │     [29] top-level scope
    │        @ none:1
    │     [30] eval(m::Module, e::Any)
    │        @ Core ./boot.jl:489
    │     [31] (::Main.var"#22#23")(chan::Channel{Any})
    │        @ Main ~/.julia/packages/QuartoNotebookRunner/evCNi/src/startup.jl:158
    │     [32] (::Base.var"#562#563"{Main.var"#22#23", Channel{Any}})()
    │        @ Base ./channels.jl:141
    └ @ IntegralProjectionModelsPadrinoExt ~/Projects/ipm/IntegralProjectionModels.jl/ext/IntegralProjectionModelsPadrinoExt.jl:345
    ┌ Warning: Failed to build IPM for aaaa19
    │   exception =
    │    UndefVarError: `crp_yr` not defined in `IntegralProjectionModelsPadrinoExt.PadrinoEvalModule`
    │    Suggestion: check for spelling errors or missing imports.
    │    Stacktrace:
    │      [1] (::IntegralProjectionModelsPadrinoExt.PadrinoEvalModule.var"#68#69")(z_new::Float64, z::Float64, h::Float64)
    │        @ IntegralProjectionModelsPadrinoExt.PadrinoEvalModule ~/Projects/ipm/IntegralProjectionModels.jl/ext/PadrinoExt/compiler.jl:484
    │      [2] build_kernel_matrix(kernel_fn::Function, domain::ContinuousDomain{Float64}, family::String3)
    │        @ IntegralProjectionModelsPadrinoExt ~/Projects/ipm/IntegralProjectionModels.jl/ext/PadrinoExt/compiler.jl:549
    │      [3] build_kernel_matrix(kernel_fn::Expr, vectorized_fn::Nothing, domain::ContinuousDomain{Float64}, family::String3)
    │        @ IntegralProjectionModelsPadrinoExt ~/Projects/ipm/IntegralProjectionModels.jl/ext/PadrinoExt/compiler.jl:609
    │      [4] _build_general_det(pm::PadrinoModel, tspan::Tuple{Int64, Int64}, mtype::IntegralProjectionModelsPadrinoExt.ModelTypeInfo)
    │        @ IntegralProjectionModelsPadrinoExt ~/Projects/ipm/IntegralProjectionModels.jl/ext/PadrinoExt/builder.jl:482
    │      [5] build_ipm(pm::PadrinoModel; tspan::Tuple{Int64, Int64})
    │        @ IntegralProjectionModelsPadrinoExt ~/Projects/ipm/IntegralProjectionModels.jl/ext/PadrinoExt/builder.jl:27
    │      [6] _pdb_make_ipm_from_protos(protos::Dict{String, PadrinoModel}; tspan::Tuple{Int64, Int64})
    │        @ IntegralProjectionModelsPadrinoExt ~/Projects/ipm/IntegralProjectionModels.jl/ext/IntegralProjectionModelsPadrinoExt.jl:337
    │      [7] _pdb_make_ipm_from_protos
    │        @ ~/Projects/ipm/IntegralProjectionModels.jl/ext/IntegralProjectionModelsPadrinoExt.jl:337 [inlined]
    │      [8] #pdb_make_ipm#69
    │        @ ~/Projects/ipm/IntegralProjectionModels.jl/ext/IntegralProjectionModelsPadrinoExt.jl:80 [inlined]
    │      [9] top-level scope
    │        @ ~/Projects/ipm/IntegralProjectionModels.jl/vignettes/10_padrino/10_padrino.qmd:95
    │     [10] eval(m::Module, e::Any)
    │        @ Core ./boot.jl:489
    │     [11] (::QuartoNotebookWorker.var"#21#22"{Module, Expr})()
    │        @ QuartoNotebookWorker ~/.julia/packages/QuartoNotebookRunner/evCNi/src/QuartoNotebookWorker/src/render.jl:222
    │     [12] (::QuartoNotebookWorker.Packages.IOCapture.var"#12#13"{Type{InterruptException}, QuartoNotebookWorker.var"#21#22"{Module, Expr}, IOContext{Base.PipeEndpoint}, IOContext{Base.PipeEndpoint}, IOContext{Base.PipeEndpoint}, IOContext{Base.PipeEndpoint}})()
    │        @ QuartoNotebookWorker.Packages.IOCapture ~/.julia/packages/QuartoNotebookRunner/evCNi/src/QuartoNotebookWorker/src/vendor/IOCapture/src/IOCapture.jl:170
    │     [13] with_logstate(f::QuartoNotebookWorker.Packages.IOCapture.var"#12#13"{Type{InterruptException}, QuartoNotebookWorker.var"#21#22"{Module, Expr}, IOContext{Base.PipeEndpoint}, IOContext{Base.PipeEndpoint}, IOContext{Base.PipeEndpoint}, IOContext{Base.PipeEndpoint}}, logstate::Base.CoreLogging.LogState)
    │        @ Base.CoreLogging ./logging/logging.jl:542
    │     [14] with_logger(f::Function, logger::Base.CoreLogging.ConsoleLogger)
    │        @ Base.CoreLogging ./logging/logging.jl:653
    │     [15] capture(f::QuartoNotebookWorker.var"#21#22"{Module, Expr}; rethrow::Type, color::Bool, passthrough::Bool, capture_buffer::IOBuffer, io_context::Vector{Pair{Symbol, Any}})
    │        @ QuartoNotebookWorker.Packages.IOCapture ~/.julia/packages/QuartoNotebookRunner/evCNi/src/QuartoNotebookWorker/src/vendor/IOCapture/src/IOCapture.jl:167
    │     [16] capture
    │        @ ~/.julia/packages/QuartoNotebookRunner/evCNi/src/QuartoNotebookWorker/src/render.jl:248 [inlined]
    │     [17] io_capture(f::Function; cell_options::Dict{String, Any}, kws::@Kwargs{rethrow::DataType, color::Bool, io_context::Vector{Pair{Symbol, Any}}})
    │        @ QuartoNotebookWorker ~/.julia/packages/QuartoNotebookRunner/evCNi/src/QuartoNotebookWorker/src/render.jl:250
    │     [18] io_capture
    │        @ ~/.julia/packages/QuartoNotebookRunner/evCNi/src/QuartoNotebookWorker/src/render.jl:246 [inlined]
    │     [19] include_str(mod::Module, code::String; file::String, line::Int64, cell_options::Dict{String, Any})
    │        @ QuartoNotebookWorker ~/.julia/packages/QuartoNotebookRunner/evCNi/src/QuartoNotebookWorker/src/render.jl:201
    │     [20] #invokelatest_gr#232
    │        @ ./reflection.jl:1297 [inlined]
    │     [21] invokelatest_gr
    │        @ ./reflection.jl:1289 [inlined]
    │     [22] #7
    │        @ ~/.julia/packages/QuartoNotebookRunner/evCNi/src/QuartoNotebookWorker/src/render.jl:18 [inlined]
    │     [23] with_inline_display(f::QuartoNotebookWorker.var"#7#8"{String, String, Int64, Dict{String, Any}}, cell_options::Dict{String, Any})
    │        @ QuartoNotebookWorker ~/.julia/packages/QuartoNotebookRunner/evCNi/src/QuartoNotebookWorker/src/InlineDisplay.jl:31
    │     [24] _render_thunk(thunk::Function, code::String, cell_options::Dict{String, Any}, is_expansion_ref::Base.RefValue{Bool}; inline::Bool)
    │        @ QuartoNotebookWorker ~/.julia/packages/QuartoNotebookRunner/evCNi/src/QuartoNotebookWorker/src/render.jl:43
    │     [25] _render_thunk
    │        @ ~/.julia/packages/QuartoNotebookRunner/evCNi/src/QuartoNotebookWorker/src/render.jl:35 [inlined]
    │     [26] #render#4
    │        @ ~/.julia/packages/QuartoNotebookRunner/evCNi/src/QuartoNotebookWorker/src/render.jl:15 [inlined]
    │     [27] render(code::String, file::String, line::Int64, cell_options::Dict{String, Any})
    │        @ QuartoNotebookWorker ~/.julia/packages/QuartoNotebookRunner/evCNi/src/QuartoNotebookWorker/src/render.jl:1
    │     [28] render(::String, ::Vararg{Any}; kwargs::@Kwargs{})
    │        @ Main ~/.julia/packages/QuartoNotebookRunner/evCNi/src/startup.jl:145
    │     [29] top-level scope
    │        @ none:1
    │     [30] eval(m::Module, e::Any)
    │        @ Core ./boot.jl:489
    │     [31] (::Main.var"#22#23")(chan::Channel{Any})
    │        @ Main ~/.julia/packages/QuartoNotebookRunner/evCNi/src/startup.jl:158
    │     [32] (::Base.var"#562#563"{Main.var"#22#23", Channel{Any}})()
    │        @ Base ./channels.jl:141
    └ @ IntegralProjectionModelsPadrinoExt ~/Projects/ipm/IntegralProjectionModels.jl/ext/IntegralProjectionModelsPadrinoExt.jl:345
    Successfully built: 6 / 20

``` julia
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

      aaaa34: λ = 0.8859
      aaaa35: λ = 1.0641
      aaaa36: λ = 1.0975

## Subsetting the Database

You can subset the database to work with specific models:

``` julia
sub_pdb = pdb_subset(pdb, ["aaa310", "aaa341"])
println("Models in subset: ", nrow(sub_pdb.tables["Metadata"]))
species = pdb_species(sub_pdb)
println("Species: ", species)
```

    Models in subset: 2
    Species: ["Aconitum_noveboracense", "Lonicera_maackii"]

## Saving and Loading

Save the database locally for offline use:

``` julia
pdb_save(pdb, joinpath(tempdir(), "padrino_data"))
pdb2 = pdb_load(joinpath(tempdir(), "padrino_data"))
```

## Notes

- **Expression parsing**: PADRINO stores vital rates as R expressions.
  The extension includes a complete R expression parser that translates
  these to Julia functions.
- **Model types**: The extension handles simple deterministic,
  stochastic kernel-resampled, stochastic parameter-resampled,
  density-dependent, and general (multi-state) models.
- **Performance**: Kernel matrices are computed using the midpoint rule,
  matching the ipmr R package approach.
- **Eviction**: Truncated distribution and discrete extrema eviction
  corrections are supported.
