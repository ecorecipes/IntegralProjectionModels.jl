@testset "Eviction Correction" begin
    @testset "Truncated distributions" begin
        domain = ContinuousDomain(0.0, 10.0, 100)
        g = NormalGrowth(5.0, 0.0, 2.0)  # constant mean=5, sd=2

        # Without correction, mass is lost outside [0,10]
        z = meshpoints(domain)
        h = step_size(domain)

        # Raw growth column sum (should be < 1 due to truncation)
        raw_sum = sum(g(zi, 5.0) * h for zi in z)
        @test raw_sum < 1.0

        # With truncated correction, should integrate to ≈ 1
        corrected_sum = sum(truncated_growth(g, zi, 5.0, domain) * h for zi in z)
        @test corrected_sum ≈ 1.0 atol = 0.02

        # PKernel with truncation
        s = ConstantSurvival(1.0)
        P_trunc = PKernel(s, g, domain; eviction = TruncatedDistributions)
        K = materialize(P_trunc)
        # Column sums should be close to 1 (survival = 1)
        for j in 1:100
            @test sum(K[:, j]) ≈ 1.0 atol = 0.02
        end
    end

    @testset "Discrete extrema" begin
        K = rand(10, 10) .* 0.08  # small values, columns sum < 1
        K_copy = copy(K)
        apply_discrete_extrema!(K)

        # Columns in first half should have mass added to row 1
        # Columns in second half should have mass added to last row
        for j in 1:5
            deficit = 1.0 - sum(K_copy[:, j])
            @test K[1, j] ≈ K_copy[1, j] + deficit
        end
        for j in 6:10
            deficit = 1.0 - sum(K_copy[:, j])
            @test K[10, j] ≈ K_copy[10, j] + deficit
        end
    end

    @testset "Discrete extrema does not renormalize fecundity/custom kernels" begin
        domain = ContinuousDomain(0.0, 1.0, 5)
        fec = FKernel(CustomVitalRate((z_prime, z) -> 3.0), domain; eviction = DiscreteExtrema)
        custom = CustomKernel((z_prime, z) -> 3.0, domain; eviction = DiscreteExtrema)

        @test materialize(fec) ≈ materialize(FKernel(CustomVitalRate((z_prime, z) -> 3.0), domain))
        @test materialize(custom) ≈ materialize(CustomKernel((z_prime, z) -> 3.0, domain))
    end
end
