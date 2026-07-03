@testset "Age×Size Models" begin
    @testset "AgeStructure" begin
        age = AgeStructure(5)
        @test n_ages(age) == 5
        @test ages(age) == 1:5
        @test_throws ArgumentError AgeStructure(0)
    end

    @testset "Age kernel expansion" begin
        domain = ContinuousDomain(0.0, 10.0, 20)
        age_struct = AgeStructure(3)
        m = n_states(domain)

        # Age-specific P and F kernels
        function p_func(a)
            # Survival increases with age (up to max)
            s = LinearSurvival(1.0 + 0.2 * a, 0.1)
            g = NormalGrowth(0.5, 0.9, 0.5)
            PKernel(s, g, domain)
        end

        function f_func(a)
            # Only adults (age > 1) reproduce, fecundity increases with age
            if a <= 1
                return CustomKernel((z_prime, z) -> 0.0, domain)
            end
            f = FecundityRate(0.01 * a, 0.005, 2.0, 0.5, 1.0)
            FKernel(f, domain)
        end

        K = expand_age_kernels(p_func, f_func, age_struct, domain)

        # Total dimension should be n_ages * n_meshpoints
        total = n_ages(age_struct) * m
        @test size(K) == (total, total)

        # Sub-diagonal blocks (aging): P_a goes from block a to block a+1
        # Verify P_1 is in block (2,1)
        P1 = materialize(p_func(1))
        @test K[(m + 1):(2m), 1:m] ≈ P1

        # Max age stays: P_3 is in block (3,3) diagonal
        P3 = materialize(p_func(3))
        @test K[(2m + 1):(3m), (2m + 1):(3m)] ≈ P3

        # Fecundity: F_2 goes from block 2 to block 1
        F2 = materialize(f_func(2))
        @test K[1:m, (m + 1):(2m)] ≈ F2

        # Age 1 doesn't reproduce: F_1 block should be zero
        @test all(K[1:m, 1:m] .== 0.0) || maximum(abs.(K[1:m, 1:m])) < 1e-15
    end

    @testset "Age kernel expansion preserves age-1 fecundity at max age" begin
        domain = ContinuousDomain(0.0, 10.0, 10)
        age_struct = AgeStructure(1)
        m = n_states(domain)

        p_func(a) = MatrixKernel(fill(0.2, m, m))
        f_func(a) = MatrixKernel(fill(0.1, m, m))

        K = expand_age_kernels(p_func, f_func, age_struct, domain)
        @test K ≈ fill(0.3, m, m)
    end

    @testset "Age×size IPM solve" begin
        domain = ContinuousDomain(0.0, 10.0, 20)
        age_struct = AgeStructure(3)
        m = n_states(domain)
        total = n_ages(age_struct) * m

        p_func(a) = PKernel(LinearSurvival(1.0 + 0.2 * a, 0.1),
            NormalGrowth(0.5, 0.9, 0.5), domain)
        f_func(a) = a > 1 ?
                     FKernel(FecundityRate(0.01 * a, 0.005, 2.0, 0.5, 1.0), domain) :
                     CustomKernel((z_prime, z) -> 0.0, domain)

        K = expand_age_kernels(p_func, f_func, age_struct, domain)

        # Wrap the matrix in a CustomKernel for the solver
        # (since expand_age_kernels returns a matrix directly)
        n0 = ones(total) ./ total
        prob = IPMProblem(SimpleIPM(), DensityIndependent(), Deterministic(),
            CustomKernel((i, j) -> K[i, j], domain), domain, n0, (0, 50);
            uses_age = true)

        # Direct matrix iteration using the expanded kernel
        u = Vector{Vector{Float64}}(undef, 51)
        u[1] = copy(n0)
        for t in 1:50
            u[t + 1] = K * u[t]
        end

        # Lambda should converge
        λ_t = [sum(u[t + 1]) / sum(u[t]) for t in 1:50]
        @test abs(λ_t[end] - λ_t[end - 1]) < 1e-6
    end
end
