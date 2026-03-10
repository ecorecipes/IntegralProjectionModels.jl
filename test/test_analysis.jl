@testset "Analysis" begin
    domain = ContinuousDomain(0.0, 50.0, 100)
    s = LinearSurvival(2.2, 0.25)
    g = NormalGrowth(0.2, 1.02, 0.7)
    f = FecundityRate(0.003, 0.015, 2.0, 0.3, 1.0)

    P = PKernel(s, g, domain)
    F = FKernel(f, domain)
    kernel = P + F
    n0 = uniform_population(domain)

    sol_iter = solve(IPMProblem(kernel, domain, n0, (0, 200)), DirectIteration())
    sol_eigen = solve(IPMProblem(kernel, domain, n0, (0, 200)), EigenAnalysis())

    @testset "lambda" begin
        λ_iter = lambda(sol_iter)
        λ_eigen = lambda(sol_eigen)
        @test λ_eigen > 0
        @test λ_iter ≈ λ_eigen atol = 1e-4
    end

    @testset "stable_distribution" begin
        w = stable_distribution(sol_eigen)
        @test length(w) == 100
        @test sum(w) ≈ 1.0 atol = 1e-10
        @test all(w .≥ 0)  # should be non-negative (Perron-Frobenius)
    end

    @testset "reproductive_value" begin
        v = reproductive_value(sol_eigen)
        w = stable_distribution(sol_eigen)
        @test length(v) == 100
        @test dot(v, w) ≈ 1.0 atol = 1e-6  # normalized so dot(v,w) = 1
    end

    @testset "sensitivity" begin
        S = sensitivity(sol_eigen)
        @test size(S) == (100, 100)
        # All entries should be non-negative for a non-negative matrix
        @test all(S .≥ -1e-10)
    end

    @testset "elasticity" begin
        E = elasticity(sol_eigen)
        @test size(E) == (100, 100)
        # Elasticities should sum to ≈ 1
        @test sum(E) ≈ 1.0 atol = 0.01
    end

    @testset "stochastic_growth_rate" begin
        # Build a stochastic model
        small_domain = ContinuousDomain(0.0, 50.0, 50)
        n0_s = uniform_population(small_domain)

        function build_kernel(params)
            s = LinearSurvival(params.s_int, 0.25)
            g = NormalGrowth(0.2, 1.02, 0.7)
            f = FecundityRate(0.003, 0.015, 2.0, 0.3, 1.0)
            PKernel(s, g, small_domain) + FKernel(f, small_domain)
        end

        env_sampler(t) = (s_int = 2.0 + 0.3 * randn(),)

        prob = IPMProblem(StochasticParameterResampled(),
            build_kernel, small_domain, n0_s, (0, 200);
            env_state = env_sampler)
        sol = solve(prob)

        λ_s = stochastic_growth_rate(sol; burn_in = 50)
        @test λ_s > 0
        @test isfinite(λ_s)
    end
end
