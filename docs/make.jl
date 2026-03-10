using Documenter
using IntegralProjectionModels

makedocs(;
    modules = [IntegralProjectionModels],
    authors = "Simon Frost <sdwfrost@gmail.com>",
    sitename = "IntegralProjectionModels.jl",
    format = Documenter.HTML(;
        prettyurls = get(ENV, "CI", "false") == "true",
        canonical = "https://sdwfrost.github.io/IntegralProjectionModels.jl"),
    pages = [
        "Home" => "index.md",
        "Tutorials" => [
            "Getting Started" => "tutorials/getting_started.md",
            "Monocarp Example" => "tutorials/monocarp_example.md",
            "Stochastic Models" => "tutorials/stochastic_models.md",
            "Bayesian Fitting" => "tutorials/bayesian_fitting.md",
        ],
    ])
