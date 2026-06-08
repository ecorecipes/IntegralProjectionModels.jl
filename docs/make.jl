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
            "Introduction to Integral Projection Models" => "tutorials/01_introduction.md",
            "Ungulate IPM: Soay Sheep" => "tutorials/02_ungulate.md",
            "Eviction Correction" => "tutorials/03_eviction.md",
            "Sensitivity and Elasticity Analysis" => "tutorials/04_sensitivity_elasticity.md",
            "Density-Dependent IPMs" => "tutorials/05_density_dependent.md",
            "Stochastic IPMs" => "tutorials/06_stochastic.md",
            "General IPM with Seed Bank" => "tutorials/07_general_ipm.md",
            "Age x Size IPM" => "tutorials/08_age_size.md",
            "Categorical Composition of IPMs" => "tutorials/09_categorical.md",
            "PADRINO Database Integration" => "tutorials/10_padrino.md",
            "Evolving IPMs: Evolutionary Demography" => "tutorials/11_evolving_ipm.md",
            "Time-Lagged Integral Projection Models" => "tutorials/12_time_lag.md",
            "Spectral Diagnostics and Type Hierarchy" => "tutorials/13_diagnostics.md",
            "Kernel Composition and Helpers" => "tutorials/14_kernels_and_helpers.md",
            "PADRINO Validation: Test Targets and λ Reproducibility" => "tutorials/15_padrino_validation.md",
            "Time-Lag Metrics: Augmentation, R0, and Generation Time" => "tutorials/16_lag_metrics.md",
            "IPM ↔ MPM Discretization Bridge" => "tutorials/17_ipm_mpm_bridge.md",
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
