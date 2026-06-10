using Random
using StructuredPopulationCore: quasi_extinction

@testset "Demographic stochasticity (binned IPM)" begin
    rng = Random.Xoshiro(2024)

    domain = ContinuousDomain(0.0, 10.0, 25)
    surv(z) = 0.7                                              # constant survival < 1
    grow(z2, z1) = exp(-((z2 - (z1 + 0.5))^2) / 2) / sqrt(2π) # N(z1+0.5, 1) growth
    fec(z2, z1) = (z1 > 3.0 ? 0.8 : 0.0) *
                  exp(-((z2 - 1.0)^2) / (2 * 0.25)) / (0.5 * sqrt(2π))   # recruits near z'=1

    P = PKernel(surv, grow, domain)
    F = FKernel(fec, domain)
    kernel = P + F
    z = meshpoints(domain)
    n0 = round.(Int, 20 .* exp.(-((z .- 5.0) .^ 2) ./ 2))     # integer initial counts

    @testset "ensemble mean tracks deterministic K-iteration" begin
        detsol = solve(IPMProblem(kernel, domain, Float64.(n0), (0, 3)), DirectIteration())

        dprob = IPMProblem(Demographic(), kernel, domain, Float64.(n0), (0, 3))
        reps = 5000
        acc = [zeros(length(z)) for _ in 1:4]
        for _ in 1:reps
            s = solve(dprob, DirectIteration(); rng=rng)
            for tt in 1:4
                acc[tt] .+= s.u[tt]
            end
        end
        for tt in 1:4
            @test isapprox(acc[tt] ./ reps, detsol.u[tt]; rtol=0.06, atol=0.25)
        end
    end

    @testset "integer counts, variance, and split requirement" begin
        dprob = IPMProblem(Demographic(), kernel, domain, Float64.(n0), (0, 3))
        s1 = solve(dprob, DirectIteration(); rng=Random.Xoshiro(1))
        s2 = solve(dprob, DirectIteration(); rng=Random.Xoshiro(2))
        @test all(x -> x == round(x) && x >= 0, s1.u[end])
        @test s1.u[end] != s2.u[end]
        @test length(s1.t) == 4

        # A CustomKernel cannot be split into survival vs fecundity -> error
        badprob = IPMProblem(Demographic(), CustomKernel((z2, z1) -> 0.5, domain),
            domain, Float64.(n0), (0, 2))
        @test_throws ErrorException solve(badprob, DirectIteration())
    end

    @testset "subcritical model goes extinct (ensemble + quasi_extinction)" begin
        kernel2 = PKernel(z -> 0.2, grow, domain) +
                  FKernel((z2, z1) -> (z1 > 3.0 ? 0.1 : 0.0) *
                          exp(-((z2 - 1.0)^2) / (2 * 0.25)) / (0.5 * sqrt(2π)), domain)
        dprob = IPMProblem(Demographic(), kernel2, domain, Float64.(n0), (0, 40))
        totals, _ = demographic_ensemble(dprob; n_reps=200, rng=rng)
        @test size(totals) == (41, 200)
        qe = quasi_extinction(totals; threshold=1.0)
        @test qe.prob_extinct > 0.5
    end
end
