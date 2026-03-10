@testset "Kernels" begin
    domain = ContinuousDomain(0.0, 50.0, 100)

    @testset "PKernel construction" begin
        s = LinearSurvival(2.2, 0.25)
        g = NormalGrowth(0.2, 1.02, 0.7)
        P = PKernel(s, g, domain)
        @test P.survival === s
        @test P.growth === g
        @test P.domain === domain
        @test P.eviction == NoCorrection

        P2 = PKernel(s, g, domain; eviction = TruncatedDistributions)
        @test P2.eviction == TruncatedDistributions
    end

    @testset "FKernel construction" begin
        f = FecundityRate(0.01, 0.0005, 3.0, 0.7, 1.0)
        F = FKernel(f, domain)
        @test F.fecundity === f
        @test F.eviction == NoCorrection
    end

    @testset "CustomKernel construction" begin
        ck = CustomKernel((z_prime, z) -> exp(-(z_prime - z)^2), domain)
        @test ck.family == CC
        @test ck.eviction == NoCorrection
    end

    @testset "ComposedKernel" begin
        s = LinearSurvival(2.2, 0.25)
        g = NormalGrowth(0.2, 1.02, 0.7)
        f = FecundityRate(0.01, 0.0005, 3.0, 0.7, 1.0)
        P = PKernel(s, g, domain)
        F = FKernel(f, domain)

        K = P + F
        @test K isa ComposedKernel
        @test length(K.subkernels) == 2

        # Adding more kernels
        ck = CustomKernel((z_prime, z) -> 0.0, domain)
        K2 = K + ck
        @test length(K2.subkernels) == 3
    end

    @testset "MegaKernel construction" begin
        s = LinearSurvival(2.2, 0.25)
        g = NormalGrowth(0.2, 1.02, 0.7)
        f = FecundityRate(0.01, 0.0005, 3.0, 0.7, 1.0)

        size_domain = ContinuousDomain(0.0, 50.0, 50)
        states = (size = size_domain,)

        P = PKernel(s, g, size_domain)
        F = FKernel(f, size_domain)

        mk = MegaKernel(; states = states, size_to_size = P + F)
        @test mk.states === states
        @test haskey(mk.kernels, (:size, :size))
    end
end
