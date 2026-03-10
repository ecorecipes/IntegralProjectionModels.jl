@testset "AD Compatibility" begin
    using ForwardDiff

    domain = ContinuousDomain(0.0, 10.0, 20)

    @testset "ForwardDiff through materialize" begin
        # Test that we can differentiate lambda w.r.t. vital rate parameters
        # using power iteration (eigen is not AD-compatible)
        function compute_lambda_power(params)
            s = LinearSurvival(params[1], params[2])
            g = NormalGrowth(params[3], params[4], params[5])
            f = FecundityRate(params[6], params[7], 2.0, 0.3, 1.0)
            P = PKernel(s, g, domain)
            F = FKernel(f, domain)
            K = materialize(P + F)
            # Power iteration to find dominant eigenvalue
            m = size(K, 1)
            v = ones(eltype(K), m) / m
            λ = one(eltype(K))
            for _ in 1:200
                w = K * v
                λ = sum(w)
                v = w / λ
            end
            return λ
        end

        params = [2.2, 0.25, 0.2, 1.02, 0.7, 0.003, 0.015]
        λ = compute_lambda_power(params)
        @test λ > 0

        # Gradient should be computable
        grad = ForwardDiff.gradient(compute_lambda_power, params)
        @test length(grad) == 7
        @test all(isfinite.(grad))

        # Sensitivity: lambda should increase with survival intercept
        @test grad[1] > 0
    end

    @testset "ForwardDiff through vital rates" begin
        # Survival
        f_surv(x) = LinearSurvival(x[1], x[2])(3.0)
        g = ForwardDiff.gradient(f_surv, [2.2, 0.25])
        @test length(g) == 2
        @test all(isfinite.(g))

        # Growth
        f_growth(x) = NormalGrowth(x[1], x[2], x[3])(5.0, 3.0)
        g = ForwardDiff.gradient(f_growth, [0.2, 1.02, 0.7])
        @test length(g) == 3
        @test all(isfinite.(g))
    end
end
