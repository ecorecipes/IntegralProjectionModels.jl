@testset "Monocarp Integration Test" begin
    # Reproduce the monocarp (Oenothera glazioviana) example from ipmr
    # Parameters from ipmr's data-raw/create_example_data.R

    # Parameters
    s_int = 1.03
    s_slope = 0.19
    g_int = 8.0
    g_slope = 0.92
    sd_g = 0.9
    f_r_int = 0.09
    f_r_slope = 0.05
    f_s_int = 0.01
    f_s_slope = 0.0005
    mu_fd = 3.0
    sd_fd = 0.7

    # Domain
    domain = ContinuousDomain(0.3, 200.0, 500)
    z = meshpoints(domain)
    h = step_size(domain)

    @testset "Manual kernel construction matches materialize" begin
        # Build P kernel: s(z) * (1 - f_r(z)) * g(z', z)
        # In the monocarp model, survival is modified by flowering probability
        # s_total(z) = logistic(s_int + s_slope * z) * (1 - logistic(f_r_int + f_r_slope * z))
        function total_survival(z_val)
            s = 1 / (1 + exp(-(s_int + s_slope * z_val)))
            f_r = 1 / (1 + exp(-(f_r_int + f_r_slope * z_val)))
            return s * (1 - f_r)
        end

        function growth(z_prime, z_val)
            mu = g_int + g_slope * z_val
            return pdf(Normal(mu, sd_g), z_prime)
        end

        # Build F kernel: f_r(z) * f_s(z) * f_d(z')
        function fecundity(z_prime, z_val)
            f_r = 1 / (1 + exp(-(f_r_int + f_r_slope * z_val)))
            f_s = exp(f_s_int + f_s_slope * z_val)
            f_d = pdf(Normal(mu_fd, sd_fd), z_prime)
            return f_r * f_s * f_d
        end

        # Manual kernel construction
        m = length(z)
        P_manual = zeros(m, m)
        F_manual = zeros(m, m)

        for j in 1:m
            s_j = total_survival(z[j])
            for i in 1:m
                P_manual[i, j] = s_j * growth(z[i], z[j]) * h
                F_manual[i, j] = fecundity(z[i], z[j]) * h
            end
        end

        K_manual = P_manual + F_manual

        # Now build with IntegralProjectionModels API
        surv = CustomVitalRate(total_survival)
        grow = CustomVitalRate(growth)
        fecund = CustomVitalRate(fecundity)

        P = PKernel(surv, grow, domain)
        F = FKernel(fecund, domain)
        K_api = materialize(P + F)

        @test K_api ≈ K_manual atol = 1e-12
    end

    @testset "Lambda from eigenanalysis" begin
        function total_survival(z_val)
            s = 1 / (1 + exp(-(s_int + s_slope * z_val)))
            f_r = 1 / (1 + exp(-(f_r_int + f_r_slope * z_val)))
            return s * (1 - f_r)
        end

        function growth(z_prime, z_val)
            return pdf(Normal(g_int + g_slope * z_val, sd_g), z_prime)
        end

        function fecundity(z_prime, z_val)
            f_r = 1 / (1 + exp(-(f_r_int + f_r_slope * z_val)))
            f_s = exp(f_s_int + f_s_slope * z_val)
            f_d = pdf(Normal(mu_fd, sd_fd), z_prime)
            return f_r * f_s * f_d
        end

        surv = CustomVitalRate(total_survival)
        grow = CustomVitalRate(growth)
        fecund = CustomVitalRate(fecundity)

        P = PKernel(surv, grow, domain)
        F = FKernel(fecund, domain)
        kernel = P + F

        n0 = uniform_population(domain)
        prob = IPMProblem(kernel, domain, n0, (0, 100))

        sol_eigen = solve(prob, EigenAnalysis())
        λ = lambda(sol_eigen)

        # Lambda should be reasonable for this model (positive, close to 1)
        @test λ > 0
        @test 0.5 < λ < 5.0  # reasonable range for a plant IPM

        # Iteration should converge to same lambda
        sol_iter = solve(prob, DirectIteration())
        @test sol_iter.lambdas[end] ≈ λ atol = 1e-4
    end

    @testset "Using built-in vital rate types" begin
        # Alternative: use the LogisticFecundityRate type
        # (doesn't exactly match monocarp because survival includes flowering prob,
        #  but demonstrates the API)
        surv = LinearSurvival(s_int, s_slope)
        grow = NormalGrowth(g_int, g_slope, sd_g)
        fecund = LogisticFecundityRate(f_r_int, f_r_slope, f_s_int, f_s_slope, mu_fd, sd_fd)

        # P kernel with just linear survival (not accounting for flowering)
        P = PKernel(surv, grow, domain)
        F = FKernel(fecund, domain)
        kernel = P + F

        n0 = uniform_population(domain)
        prob = IPMProblem(kernel, domain, n0, (0, 100))
        sol = solve(prob, EigenAnalysis())

        # This won't match the monocarp exactly (different survival formulation)
        # but lambda should still be positive and reasonable
        @test lambda(sol) > 0
    end

    @testset "Sensitivity and elasticity" begin
        function total_survival(z_val)
            s = 1 / (1 + exp(-(s_int + s_slope * z_val)))
            f_r = 1 / (1 + exp(-(f_r_int + f_r_slope * z_val)))
            return s * (1 - f_r)
        end

        function growth(z_prime, z_val)
            return pdf(Normal(g_int + g_slope * z_val, sd_g), z_prime)
        end

        function fecundity(z_prime, z_val)
            f_r = 1 / (1 + exp(-(f_r_int + f_r_slope * z_val)))
            f_s = exp(f_s_int + f_s_slope * z_val)
            return f_r * f_s * pdf(Normal(mu_fd, sd_fd), z_prime)
        end

        P = PKernel(CustomVitalRate(total_survival), CustomVitalRate(growth), domain)
        F = FKernel(CustomVitalRate(fecundity), domain)
        kernel = P + F

        n0 = uniform_population(domain)
        sol = solve(IPMProblem(kernel, domain, n0, (0, 100)), EigenAnalysis())

        S = sensitivity(sol)
        E = elasticity(sol)

        @test size(S) == (500, 500)
        @test size(E) == (500, 500)
        @test sum(E) ≈ 1.0 atol = 0.01
    end
end
