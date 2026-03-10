@testset "Time Lag" begin

    # Set up a simple IPM with known P and F kernels
    domain = ContinuousDomain(0.0, 5.0, 20)
    mesh = meshpoints(domain)
    h = step_size(domain)

    # Simple survival: logistic with size
    survival = LinearSurvival(-0.5, 0.3)
    # Simple growth: normal with linear mean
    growth = NormalGrowth(0.5, 0.8, 0.5)
    # Simple fecundity
    fec = FecundityRate(0.003, 0.015, 2.0, 0.3, 1.0)

    P = PKernel(survival, growth, domain)
    F = FKernel(fec, domain)

    @testset "LaggedKernel construction" begin
        lk = LaggedKernel(P, F)
        @test lk isa LaggedKernel
        @test lk.lag_structure.max_lag == 1
        @test lk.immediate === P
        @test haskey(lk.lagged, 1)
        @test lk.lagged[1] === F
    end

    @testset "LaggedKernel materialize" begin
        lk = LaggedKernel(P, F)
        K_aug = materialize(lk)
        m = length(mesh)
        @test size(K_aug) == (2m, 2m)

        # Top-left block should be materialized P
        P_mat = materialize(P)
        @test K_aug[1:m, 1:m] ≈ P_mat

        # Top-right block should be materialized F
        F_mat = materialize(F)
        @test K_aug[1:m, m+1:2m] ≈ F_mat

        # Bottom-left block should be identity
        @test K_aug[m+1:2m, 1:m] ≈ Matrix{Float64}(I, m, m)

        # Bottom-right block should be zero
        @test K_aug[m+1:2m, m+1:2m] ≈ zeros(m, m)
    end

    @testset "LaggedKernel show" begin
        lk = LaggedKernel(P, F)
        s = sprint(show, lk)
        @test occursin("LaggedKernel", s)
        @test occursin("max_lag=1", s)
    end

    @testset "expand_lag_kernels convenience" begin
        K_aug = expand_lag_kernels(P, F, domain)
        m = length(mesh)
        @test size(K_aug) == (2m, 2m)

        # Should match LaggedKernel materialize
        lk = LaggedKernel(P, F)
        @test K_aug ≈ materialize(lk)
    end

    @testset "DirectIteration with LaggedKernel" begin
        lk = LaggedKernel(P, F)
        m = length(mesh)
        n0 = ones(m) ./ m
        prob = IPMProblem(SimpleIPM(), DensityIndependent(), Deterministic(),
                          lk, domain, n0, (0, 200))
        sol = solve(prob, DirectIteration())

        @test sol.retcode == :Success
        @test length(sol.t) == 201
        @test length(sol.u) == 201
        @test length(sol.u[end]) == m  # physical state, not augmented
        @test length(sol.lambdas) == 200

        # Converged lambda should match eigenanalysis
        λ_iter = sol.lambdas[end]
        K_aug = materialize(lk)
        λ_eigen = lambda(K_aug)
        @test λ_iter ≈ λ_eigen atol=0.01
    end

    @testset "EigenAnalysis with LaggedKernel" begin
        lk = LaggedKernel(P, F)
        m = length(mesh)
        n0 = ones(m) ./ m
        prob = IPMProblem(SimpleIPM(), DensityIndependent(), Deterministic(),
                          lk, domain, n0, (0, 50))
        sol = solve(prob, EigenAnalysis())

        @test sol.retcode == :Success
        @test sol.eigenanalysis !== nothing
        @test sol.eigenanalysis.lambda > 0
    end

    @testset "Lagged vs standard eigenvalue" begin
        # Lagged model should have different lambda from standard P+F model
        lk = LaggedKernel(P, F)
        K_aug = materialize(lk)
        λ_lagged = lambda(K_aug)

        K_standard = materialize(P) .+ materialize(F)
        λ_standard = lambda(K_standard)

        @test λ_lagged > 0
        @test isfinite(λ_lagged)
        @test !isapprox(λ_lagged, λ_standard; atol=1e-4)
    end

    @testset "Multi-lag L=2" begin
        # F at lag 2
        lagged_dict = Dict{Int, typeof(F)}(2 => F)
        lk = LaggedKernel(P, lagged_dict, TimeLagStructure(2))
        K_aug = materialize(lk)
        m = length(mesh)
        @test size(K_aug) == (3m, 3m)

        λ = lambda(K_aug)
        @test λ > 0
        @test isfinite(λ)

        # DirectIteration should work
        n0 = ones(m) ./ m
        prob = IPMProblem(SimpleIPM(), DensityIndependent(), Deterministic(),
                          lk, domain, n0, (0, 50))
        sol = solve(prob, DirectIteration())
        @test sol.retcode == :Success
    end

    @testset "Normalized iteration" begin
        lk = LaggedKernel(P, F)
        m = length(mesh)
        n0 = ones(m)
        prob = IPMProblem(SimpleIPM(), DensityIndependent(), Deterministic(),
                          lk, domain, n0, (0, 20); normalize=true)
        sol = solve(prob, DirectIteration())
        @test sol.retcode == :Success
        # With normalization, population should stay bounded
        @test all(sum.(sol.u) .< 1e10)
    end

end
