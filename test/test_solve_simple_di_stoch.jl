@testset "Solve: simple_di_stoch" begin
    domain = ContinuousDomain(0.0, 50.0, 50)
    n0 = uniform_population(domain)

    @testset "Kernel-resampled stochasticity" begin
        # Build a set of kernels with different parameters
        s1 = LinearSurvival(2.0, 0.2)
        s2 = LinearSurvival(2.5, 0.3)
        g = NormalGrowth(0.2, 1.02, 0.7)
        f = FecundityRate(0.003, 0.015, 2.0, 0.3, 1.0)

        kern1 = PKernel(s1, g, domain) + FKernel(f, domain)
        kern2 = PKernel(s2, g, domain) + FKernel(f, domain)

        prob = IPMProblem(StochasticKernelResampled(),
            [kern1, kern2], domain, n0, (0, 50))
        sol = solve(prob; kernel_seq = rand(1:2, 50))

        @test sol.retcode == :Success
        @test length(sol.u) == 51
        @test length(sol.lambdas) == 50
        @test all(sol.lambdas .> 0)
    end

    @testset "Parameter-resampled stochasticity" begin
        # Kernel builder function: takes sampled params, returns kernel
        function build_kernel(params)
            s = LinearSurvival(params.s_int, params.s_slope)
            g = NormalGrowth(0.2, 1.02, 0.7)
            f = FecundityRate(0.003, 0.015, 2.0, 0.3, 1.0)
            PKernel(s, g, domain) + FKernel(f, domain)
        end

        # Environment sampler: returns parameter NamedTuple each step
        function env_sampler(t)
            (s_int = 2.0 + 0.5 * randn(), s_slope = 0.25 + 0.05 * randn())
        end

        prob = IPMProblem(StochasticParameterResampled(),
            build_kernel, domain, n0, (0, 50);
            env_state = env_sampler)
        sol = solve(prob)

        @test sol.retcode == :Success
        @test length(sol.u) == 51
        @test length(sol.lambdas) == 50

        # Stochastic growth rate
        λ_s = stochastic_growth_rate(sol; burn_in = 10)
        @test λ_s > 0
    end
end
