using IntegralProjectionModels
using Test
using LinearAlgebra
using Distributions
using Statistics

@testset "IntegralProjectionModels.jl" begin
    include("test_domains.jl")
    include("test_vital_rates.jl")
    include("test_kernels.jl")
    include("test_materialize.jl")
    include("test_eviction.jl")
    include("test_solve_simple_di_det.jl")
    include("test_solve_simple_di_stoch.jl")
    include("test_solve_simple_dd.jl")
    include("test_solve_general.jl")
    include("test_age_size.jl")
    include("test_analysis.jl")
    include("test_ad_compatibility.jl")
    include("test_monocarp.jl")
    include("test_padrino.jl")
    include("test_time_lag.jl")
    include("test_demographic.jl")
end
