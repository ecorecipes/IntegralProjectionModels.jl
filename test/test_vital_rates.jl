@testset "Vital Rates" begin
    @testset "LinearSurvival" begin
        s = LinearSurvival(2.2, 0.25)
        @test s(0.0) ≈ 1 / (1 + exp(-2.2))
        @test s(1.0) ≈ 1 / (1 + exp(-(2.2 + 0.25)))
        # Survival should be in (0, 1)
        @test 0 < s(0.0) < 1
        @test 0 < s(10.0) < 1
    end

    @testset "QuadraticSurvival" begin
        s = QuadraticSurvival(1.0, 0.5, -0.01)
        @test s(0.0) ≈ 1 / (1 + exp(-1.0))
        @test 0 < s(5.0) < 1
    end

    @testset "ConstantSurvival" begin
        s = ConstantSurvival(0.9)
        @test s(0.0) == 0.9
        @test s(100.0) == 0.9
    end

    @testset "NormalGrowth" begin
        g = NormalGrowth(0.2, 1.02, 0.7)
        # Growth density should be positive
        @test g(5.0, 5.0) > 0
        # Mean size
        @test mean_size(g, 5.0) ≈ 0.2 + 1.02 * 5.0
        # Symmetric around mean
        mu = mean_size(g, 5.0)
        @test g(mu + 0.5, 5.0) ≈ g(mu - 0.5, 5.0)
    end

    @testset "LogNormalGrowth" begin
        g = LogNormalGrowth(1.0, 0.5, 0.3)
        @test g(3.0, 2.0) > 0
        @test g(-1.0, 2.0) == 0.0  # log-normal is 0 for negative values
    end

    @testset "FecundityRate" begin
        f = FecundityRate(0.01, 0.0005, 3.0, 0.7, 1.0)
        # Fecundity should be non-negative
        @test f(3.0, 5.0) ≥ 0
        @test f(3.0, 5.0) > 0  # should be positive at recruit mean
    end

    @testset "LogisticFecundityRate" begin
        f = LogisticFecundityRate(0.09, 0.05, 0.01, 0.0005, 3.0, 0.7)
        @test f(3.0, 5.0) > 0
    end

    @testset "RecruitmentDistribution" begin
        r = RecruitmentDistribution(3.0, 0.7)
        @test r(3.0) ≈ pdf(Normal(3.0, 0.7), 3.0)
        @test r(3.0, 10.0) ≈ r(3.0)  # independent of parent size
    end

    @testset "CustomVitalRate" begin
        vr = CustomVitalRate(z -> 0.5 * z)
        @test vr(2.0) ≈ 1.0

        vr2 = CustomVitalRate((z, p) -> p.a + p.b * z, (a = 1.0, b = 0.5))
        @test vr2(2.0) ≈ 2.0
    end
end
