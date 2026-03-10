@testset "Solve: general models" begin
    @testset "General DI deterministic with MegaKernel" begin
        # Two continuous states: size and reproductive_effort
        size_domain = ContinuousDomain(0.0, 20.0, 30)
        states = (size = size_domain,)

        s = LinearSurvival(2.2, 0.25)
        g = NormalGrowth(0.2, 1.02, 0.7)
        f = FecundityRate(0.003, 0.015, 2.0, 0.3, 1.0)

        mk = MegaKernel(; states = states,
            size_to_size = PKernel(s, g, size_domain) + FKernel(f, size_domain))

        n0 = uniform_population(size_domain)
        prob = IPMProblem(GeneralIPM(), mk, states, n0, (0, 50))
        sol = solve(prob)

        @test sol.retcode == :Success
        @test length(sol.u) == 51
    end

    @testset "General model: two-state (size + seedbank)" begin
        size_domain = ContinuousDomain(0.0, 10.0, 20)
        seedbank_domain = DiscreteDomain([:seedbank])
        states = (size = size_domain, seedbank = seedbank_domain)

        m_size = n_states(size_domain)  # 20
        m_seed = n_states(seedbank_domain)  # 1

        # Kernels for transitions between states
        s = LinearSurvival(1.0, 0.1)
        g = NormalGrowth(0.5, 0.9, 0.5)
        f = FecundityRate(0.01, 0.005, 2.0, 0.5, 0.5)

        # size → size: survival-growth + fecundity (20×20)
        P_ss = PKernel(s, g, size_domain)
        F_ss = FKernel(f, size_domain)

        # seedbank → size: germination (20×1 matrix, DC family)
        z = meshpoints(size_domain)
        h = step_size(size_domain)
        germ_mat = reshape([0.3 * pdf(Normal(2.0, 0.5), zi) * h for zi in z], m_size, m_seed)
        germ_kernel = MatrixKernel(germ_mat)

        # size → seedbank: seed production (1×20 matrix, CD family)
        seed_mat = reshape([0.1 * exp(0.01 * zi) * h for zi in z], m_seed, m_size)
        seed_kernel = MatrixKernel(seed_mat)

        # seedbank → seedbank: survival in seedbank (1×1)
        sb_surv = MatrixKernel(fill(0.5, 1, 1))

        mk = MegaKernel(;
            states = states,
            size_to_size = P_ss + F_ss,
            seedbank_to_size = germ_kernel,
            size_to_seedbank = seed_kernel,
            seedbank_to_seedbank = sb_surv)

        total_n = m_size + m_seed
        n0 = ones(total_n) ./ total_n
        prob = IPMProblem(GeneralIPM(), mk, states, n0, (0, 30))
        sol = solve(prob)

        @test sol.retcode == :Success
        @test length(sol.u) == 31
        @test length(sol.u[1]) == total_n
    end
end
