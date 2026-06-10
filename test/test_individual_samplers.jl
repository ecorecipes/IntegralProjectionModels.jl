using Random
using Statistics
using StructuredPopulationCore: sample_survives, sample_growth, sample_recruit,
    expected_offspring, offspring_count, sample_from_density

@testset "Individual samplers (IPM vital-rate types)" begin
    rng = Random.Xoshiro(7)
    dom = ContinuousDomain(0.0, 20.0, 200)

    @testset "NormalGrowth analytic sampler" begin
        g = NormalGrowth(0.2, 1.02, 0.7)
        z = 5.0
        gs = [sample_growth(rng, g, z, dom) for _ in 1:20000]
        @test isapprox(mean(gs), mean_size(g, z); atol=0.03)   # μ = 0.2 + 1.02·z
        @test isapprox(std(gs), 0.7; rtol=0.05)
    end

    @testset "FecundityRate count + recruit" begin
        f = FecundityRate(0.003, 0.015, 2.0, 0.3, 1.0)
        z = 5.0
        @test expected_offspring(f, z, dom) ≈ 1.0 * exp(0.003 + 0.015 * z)
        ks = [offspring_count(rng, f, z, dom) for _ in 1:20000]
        @test isapprox(mean(ks), expected_offspring(f, z, dom); rtol=0.05)
        rs = [sample_recruit(rng, f, z, dom) for _ in 1:20000]
        @test isapprox(mean(rs), 2.0; atol=0.02)
        @test isapprox(std(rs), 0.3; rtol=0.05)
    end

    @testset "survival uses the Core Bernoulli default" begin
        s = LinearSurvival(2.2, 0.25)
        xs = [sample_survives(rng, s, 1.0) for _ in 1:20000]
        @test isapprox(mean(xs), s(1.0); atol=0.02)
    end

    @testset "analytic agrees with generic density fallback" begin
        g = NormalGrowth(0.2, 1.02, 0.7); z = 5.0
        fb = [sample_from_density(rng, x -> g(x, z), dom) for _ in 1:20000]
        @test isapprox(mean(fb), mean_size(g, z); atol=0.05)
    end
end
