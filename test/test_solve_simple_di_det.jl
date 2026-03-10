@testset "Solve: simple_di_det" begin
    # Parameters from ipmr test suite
    domain = ContinuousDomain(0.0, 50.0, 100)
    s = LinearSurvival(2.2, 0.25)
    g = NormalGrowth(0.2, 1.02, 0.7)
    f = FecundityRate(0.003, 0.015, 2.0, 0.3, 1.0)

    P = PKernel(s, g, domain)
    F = FKernel(f, domain)
    kernel = P + F

    n0 = uniform_population(domain)

    @testset "DirectIteration" begin
        prob = IPMProblem(kernel, domain, n0, (0, 100))
        sol = solve(prob, DirectIteration())

        @test sol.retcode == :Success
        @test length(sol.t) == 101
        @test length(sol.u) == 101
        @test all(u -> all(u .≥ 0), sol.u)
        @test length(sol.lambdas) == 100

        # Lambda should converge
        λ_final = sol.lambdas[end]
        λ_penult = sol.lambdas[end - 1]
        @test abs(λ_final - λ_penult) < 1e-6

        # Eigenanalysis should be performed
        @test sol.eigenanalysis !== nothing
        @test sol.eigenanalysis.lambda > 0
    end

    @testset "EigenAnalysis" begin
        prob = IPMProblem(kernel, domain, n0, (0, 100))
        sol = solve(prob, EigenAnalysis())

        @test sol.retcode == :Success
        @test sol.eigenanalysis !== nothing
        @test sol.eigenanalysis.lambda > 0

        # Stable distribution should sum to 1
        @test sum(sol.eigenanalysis.stable_dist) ≈ 1.0 atol = 1e-10

        # Reproductive value normalized so dot(v, w) = 1
        @test dot(sol.eigenanalysis.repro_value, sol.eigenanalysis.stable_dist) ≈ 1.0 atol = 1e-6
    end

    @testset "DirectIteration matches EigenAnalysis" begin
        prob = IPMProblem(kernel, domain, n0, (0, 200))
        sol_iter = solve(prob, DirectIteration())
        sol_eigen = solve(prob, EigenAnalysis())

        # Lambda from iteration should converge to eigenvalue
        @test sol_iter.lambdas[end] ≈ sol_eigen.eigenanalysis.lambda atol = 1e-4
    end

    @testset "Normalized iteration" begin
        prob = IPMProblem(kernel, domain, n0, (0, 50); normalize = true)
        sol = solve(prob)

        # Population should remain normalized
        for u in sol.u
            @test sum(u) ≈ 1.0 atol = 1e-10
        end
    end

    @testset "Default algorithm" begin
        prob = IPMProblem(kernel, domain, n0, (0, 10))
        sol = solve(prob)  # should use DirectIteration by default
        @test sol.retcode == :Success
        @test length(sol.u) == 11
    end
end
