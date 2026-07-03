@testset "Solve: simple_dd" begin
    domain = ContinuousDomain(0.0, 50.0, 50)
    n0 = ones(50) .* 10.0  # start with a population

    @testset "Density-dependent deterministic" begin
        # Kernel function that depends on population state
        # As population grows, survival decreases (negative DD)
        function dd_kernel(n_t, t, p)
            total_N = sum(n_t)
            # Survival decreases with population density
            s_int = 2.2 - 0.001 * total_N
            s = LinearSurvival(s_int, 0.25)
            g = NormalGrowth(0.2, 1.02, 0.7)
            f = FecundityRate(0.003, 0.015, 2.0, 0.3, 1.0)
            PKernel(s, g, domain) + FKernel(f, domain)
        end

        prob = IPMProblem(DensityDependent(),
            dd_kernel, domain, n0, (0, 50))
        sol = solve(prob)

        @test sol.retcode == :Success
        @test length(sol.u) == 51
        @test length(sol.lambdas) == 50

        # With negative DD, population should stabilize (lambda → 1)
        # or decline depending on parameters
        @test sol.lambdas[end] < sol.lambdas[1]  # growth rate should decrease

        # Kernel matrices should be stored per step
        @test sol.kernel_matrices isa Vector
        @test length(sol.kernel_matrices) == 50
    end

    @testset "Density-dependent parameter-resampled stochasticity" begin
        function dd_stoch_kernel(n_t, t, params)
            total_N = sum(n_t)
            s = ConstantSurvival(max(0.05, params.base_survival - params.dd_strength * total_N))
            g = CustomVitalRate((z_prime, z) -> z_prime == z ? 1 / step_size(domain) : 0.0)
            PKernel(s, g, domain)
        end

        env_sampler(t) = (base_survival = 1.3, dd_strength = 0.002)

        prob = IPMProblem(DensityDependent(), StochasticParameterResampled(),
            dd_stoch_kernel, domain, n0, (0, 6); env_state = env_sampler)
        sol = solve(prob)

        @test sol.retcode == :Success
        @test length(sol.kernel_matrices) == 6
        @test maximum(abs.(sol.kernel_matrices[1] .- sol.kernel_matrices[2])) > 0
        @test sol.lambdas[2] != sol.lambdas[1]
    end
end
