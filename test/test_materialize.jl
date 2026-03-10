@testset "Materialize" begin
    domain = ContinuousDomain(0.0, 50.0, 100)
    h = step_size(domain)
    z = meshpoints(domain)
    m = length(z)

    @testset "PKernel materialization" begin
        s = LinearSurvival(2.2, 0.25)
        g = NormalGrowth(0.2, 1.02, 0.7)
        P = PKernel(s, g, domain)
        K = materialize(P)

        @test size(K) == (m, m)
        @test all(K .≥ 0)  # all entries non-negative

        # Check a specific entry: K[i,j] = h * s(z[j]) * g(z[i], z[j])
        i, j = 50, 30
        expected = h * s(z[j]) * g(z[i], z[j])
        @test K[i, j] ≈ expected

        # Column sums should be ≤ 1 (survival probability bounds)
        # In practice columns in the bulk of the domain should be close to s(z_j)
    end

    @testset "FKernel materialization" begin
        f = FecundityRate(0.01, 0.0005, 3.0, 0.7, 1.0)
        F = FKernel(f, domain)
        K = materialize(F)

        @test size(K) == (m, m)
        @test all(K .≥ 0)

        # Check a specific entry
        i, j = 10, 50
        expected = h * f(z[i], z[j])
        @test K[i, j] ≈ expected
    end

    @testset "CustomKernel materialization" begin
        func = (z_prime, z) -> exp(-(z_prime - z)^2 / 2)
        ck = CustomKernel(func, domain)
        K = materialize(ck)

        @test size(K) == (m, m)
        i, j = 50, 50
        @test K[i, j] ≈ h * func(z[i], z[j])
    end

    @testset "ComposedKernel materialization" begin
        s = LinearSurvival(2.2, 0.25)
        g = NormalGrowth(0.2, 1.02, 0.7)
        f = FecundityRate(0.01, 0.0005, 3.0, 0.7, 1.0)
        P = PKernel(s, g, domain)
        F = FKernel(f, domain)

        K_composed = materialize(P + F)
        K_P = materialize(P)
        K_F = materialize(F)

        @test K_composed ≈ K_P + K_F
    end

    @testset "MegaKernel materialization" begin
        s = LinearSurvival(2.2, 0.25)
        g = NormalGrowth(0.2, 1.02, 0.7)

        small_domain = ContinuousDomain(0.0, 10.0, 20)
        states = (size = small_domain,)

        P = PKernel(s, g, small_domain)
        mk = MegaKernel(; states = states, size_to_size = P)
        K = materialize(mk)

        @test size(K) == (20, 20)
        @test K ≈ materialize(P)
    end

    @testset "Midpoint rule integration accuracy" begin
        # Test that midpoint rule integrates a constant function correctly
        # If g(z',z) = 1/L (uniform) and s(z) = 1, then column sum * h should ≈ 1
        L = 10.0
        d = ContinuousDomain(0.0, L, 200)
        s_const = ConstantSurvival(1.0)
        # Growth that's uniform over domain
        g_uniform = CustomVitalRate((z_prime, z) -> 1.0 / L)
        P = PKernel(s_const, g_uniform, d)
        K = materialize(P)
        # Each column sum should ≈ 1 (survival=1, growth integrates to 1)
        for j in 1:200
            @test sum(K[:, j]) ≈ 1.0 atol = 0.01
        end
    end
end
