# Functions shared with the rest of the stack: one function each, extended here for this package's types.
import StructuredPopulationCore
import SciMLBase

@testset "functions shared with the stack" begin
    core = StructuredPopulationCore
    for name in (:left_kan_extension, :right_kan_extension, :stratify, :coarsen, :compose_from_uwd,
                 :demographic_ensemble, :to_discrete_problem, :meshpoints, :step_size, :bounds, :n_states)
        @test getfield(IntegralProjectionModels, name) === getfield(core, name)
    end
    @test remake === SciMLBase.remake
    # this package's methods are on those functions
    @test hasmethod(demographic_ensemble, Tuple{IPMProblem})
    @test hasmethod(remake, Tuple{IPMProblem})
    @test any(m -> m.module === IntegralProjectionModels, methods(to_discrete_problem))

    # the categorical operations on the shared domain type
    domain = ContinuousDomain(0.0, 10.0, 20)
    z, h = meshpoints(domain), step_size(domain)
    k(y, x) = exp(-(y - x)^2)
    A = left_kan_extension(k, domain)
    @test A == [h * k(y, x) for y in z, x in z]
    @test compose_kernels((a = k, b = k), domain) ≈ 2A
    @test right_kan_extension(A, domain)(z[4], z[9]) == A[4, 9] / h
    D = [0.9 0.1; 0.1 0.9]
    @test stratify(A, D) == kron(D, A)
    C = coarsen(A, domain, ContinuousDomain(0.0, 10.0, 5))
    @test size(C) == (5, 5) && C[1, 2] ≈ sum(A[1:4, 5:8]) / 4
end
