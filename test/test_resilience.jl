import ResilienceCore
import StructuredPopulationCore
using StructuredPopulationCore: ContinuousDomain, DiscreteDomain, step_size, meshpoints

@testset "linearise: resilience operators" begin
    s = LinearSurvival(2.2, 0.25)
    g = NormalGrowth(0.2, 1.02, 0.7)
    f = FecundityRate(0.003, 0.015, 2.0, 0.3, 1.0)
    function ipm(n)
        domain = ContinuousDomain(0.0, 50.0, n)
        return IPMProblem(PKernel(s, g, domain) + FKernel(f, domain), domain, ones(n), (0, 10))
    end

    op = linearise(ipm(100))
    @test op isa ResilienceCore.ProjectionMatrix
    @test op.weights ≈ fill(0.5, 100)
    @test ResilienceCore.growth_rate(op) ≈ StructuredPopulationCore.lambda(materialize(ipm(100).kernel))

    # the metrics approximate those of the integral operator, so they converge with the mesh
    fine = linearise(ipm(400))
    @test ResilienceCore.growth_rate(op) ≈ ResilienceCore.growth_rate(fine) rtol = 1e-3
    @test ResilienceCore.reactivity(op) ≈ ResilienceCore.reactivity(fine) rtol = 1e-2
    @test ResilienceCore.max_amplification(op)[1] ≈ ResilienceCore.max_amplification(fine)[1] rtol = 1e-2

    # a general IPM with a discrete seedbank: weights h for the mesh, 1 for the seedbank
    size_domain = ContinuousDomain(0.0, 10.0, 20)
    states = (size = size_domain, seedbank = DiscreteDomain([:seedbank]))
    h = step_size(size_domain)
    z = meshpoints(size_domain)
    mk = MegaKernel(; states,
                    size_to_size = PKernel(LinearSurvival(1.0, 0.1), NormalGrowth(0.5, 0.9, 0.5), size_domain),
                    seedbank_to_size = MatrixKernel(reshape(0.3 .* exp.(-(z .- 2) .^ 2), 20, 1)),
                    size_to_seedbank = MatrixKernel(reshape(0.1 .* exp.(0.01 .* z) .* h, 1, 20)),
                    seedbank_to_seedbank = MatrixKernel(fill(0.5, 1, 1)))
    general = linearise(IPMProblem(GeneralIPM(), mk, states, ones(21), (0, 10)))
    @test general.weights ≈ [fill(h, 20); 1.0]
    @test ResilienceCore.growth_rate(general) ≈ StructuredPopulationCore.lambda(materialize(mk))
end
