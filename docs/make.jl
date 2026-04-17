using Documenter
using IntegralProjectionModels
using StructuredPopulationCore

makedocs(;
    modules = [IntegralProjectionModels, StructuredPopulationCore],
    warnonly = true,
    authors = "Simon Frost",
    sitename = "IntegralProjectionModels.jl",
    format = Documenter.HTML(;
        prettyurls = get(ENV, "CI", "false") == "true",
        canonical = "https://ecorecipes.github.io/IntegralProjectionModels.jl"),
    pages = [
        "Home" => "index.md",
        "Tutorials" => [
            "Getting Started" => "tutorials/getting_started.md",
            "Monocarp Example" => "tutorials/monocarp_example.md",
            "Stochastic Models" => "tutorials/stochastic_models.md",
            "Bayesian Fitting" => "tutorials/bayesian_fitting.md",
        ],
        "API Reference" => [
            "Domains" => "api/domains.md",
            "Vital Rates" => "api/vital_rates.md",
            "Kernels" => "api/kernels.md",
            "Types & Traits" => "api/types.md",
            "Problem & Solution" => "api/problems.md",
            "Analysis" => "api/analysis.md",
            "Time-Lag Models" => "api/time_lag.md",
            "Age×Size Models" => "api/age_size.md",
            "Categorical Composition" => "api/categorical.md",
            "SciML Interface" => "api/sciml.md",
            "Utilities" => "api/utilities.md",
            "PADRINO Integration" => "api/padrino.md",
        ],
    ])

deploydocs(;
    repo = "github.com/ecorecipes/IntegralProjectionModels.jl.git",
)
