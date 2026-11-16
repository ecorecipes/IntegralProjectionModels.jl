using Random
import MatrixProjectionModels
using StructuredPopulationCore: rand_poisson, demographic_step!

@testset "Coupled IPMs" begin
    domain = ContinuousDomain(0.0, 10.0, 30)
    z = meshpoints(domain)
    surv(x) = 0.2 + 0.07x
    grow(y, x) = exp(-((y - (1.0 + 0.8x))^2) / 2) / sqrt(2π)
    fec(y, x) = (x > 4.0 ? 0.3 : 0.0) * exp(-((y - 1.0)^2) / 0.5) / sqrt(0.5π)
    kernel = PKernel(surv, grow, domain) + FKernel(fec, domain)
    P, F = materialize(PKernel(surv, grow, domain)), materialize(FKernel(fec, domain))
    n0 = 10.0 .* exp.(-((z .- 3.0) .^ 2))

    @testset "a discretised kernel" begin
        states = (seed = DiscreteDomain([:seed]), plant = domain)
        U = zeros(31, 31); U[2:end, 2:end] = P; U[1, 1] = 0.4
        Fs = zeros(31, 31); Fs[1, 2:end] .= 0.1
        k = DiscretizedKernel(U, Fs, states)
        @test materialize(k) == U .+ Fs
        @test k.states === states
        @test IntegralProjectionModels.component_stage_names(k, 31)[1:3] == [:seed, :plant1, :plant2]
        @test_throws DimensionMismatch DiscretizedKernel(U, Fs, (plant = domain,))
        @test_throws DimensionMismatch DiscretizedKernel(U, zeros(30, 30), states)
        @test DiscretizedKernel(ones(Int, 2, 2), zeros(2, 2), (a = DiscreteDomain([:x, :y]),)).survival isa Matrix{Float64}
        # a single population can be solved with it, deterministically and with whole individuals
        whole = DiscretizedKernel(P, F, (plant = domain,))
        @test solve(IPMProblem(whole, domain, n0, (0, 5))).u == solve(IPMProblem(kernel, domain, n0, (0, 5))).u
        d = solve(IPMProblem(Demographic(), whole, domain, round.(n0), (0, 3)), DirectIteration(); rng = Xoshiro(3))
        @test all(isinteger, d.u[end])
    end

    @testset "one population: the same as an IPMProblem" begin
        single = solve(IPMProblem(kernel, domain, n0, (0, 20)), DirectIteration())
        sys = PopulationSystem(:plant => PopulationComponent(kernel, n0))
        prob = CoupledIPMProblem(sys, (0, 20))
        sol = solve(prob, DirectIteration())
        @test sol isa CoupledIPMSolution && sol.t == 0:20 && sol.retcode === :Success
        @test all(sol.u[k][:plant] == single.u[k] for k in 1:21)
        @test sol[:plant] ≈ sum.(single.u)
        @test sol.component_matrices[:plant][1] == materialize(kernel)
        @test all(isempty, values(solve(prob; save_matrices = false).component_matrices))
        # the kernel can be given by a function of the system, a MegaKernel, or a discretised kernel
        for model in ((s, day, p) -> kernel, MegaKernel(Dict((:plant, :plant) => kernel), (plant = domain,)),
                      DiscretizedKernel(P, F, (plant = domain,)), materialize(kernel))
            other = solve(CoupledIPMProblem(PopulationSystem(:plant => PopulationComponent(model, n0)), (0, 20)))
            @test other.u[end][:plant] ≈ single.u[end] rtol = 1e-12
        end
        @test_throws ArgumentError solve(prob, EigenAnalysis())
        @test remake(prob; tspan = (0, 3)).tspan == (0, 3)
        @test occursin("CoupledIPMProblem(1 components", sprint(show, prob))
        @test_throws ArgumentError CoupledIPMProblem(sys, (0, 2); inflows = Dict(:tree => (s, d, p) -> zeros(30)))
        @test_throws ArgumentError CoupledIPMProblem(sys, (0, 2); stochasticity = StochasticKernelResampled())
        @test_throws DimensionMismatch solve(CoupledIPMProblem(PopulationSystem(:plant => PopulationComponent(
            (s, day, p) -> zeros(3, 3), n0)), (0, 1)))
    end

    # two populations whose kernels read each other's numbers through the system
    crowding(s) = 1 / (1 + 0.01 * (sum(s[:plant].population) + sum(s[:shrub].population)))
    stages = (seedling = DiscreteDomain([:seedling]), adult = DiscreteDomain([:adult]))
    Us = [0.0 0.0; 0.3 0.9]
    Fs = [0.0 2.0; 0.0 0.0]
    function coupled(target)
        plant = (s, day, p) -> target === :ipm ? DiscretizedKernel(P, crowding(s) .* F, (plant = domain,)) :
                                                MatrixProjectionModels.MatrixProjectionModel(P, crowding(s) .* F)
        shrub = (s, day, p) -> target === :ipm ? DiscretizedKernel(Us, p .* crowding(s) .* Fs, stages) :
                                                MatrixProjectionModels.MatrixProjectionModel(Us, p .* crowding(s) .* Fs)
        sys = PopulationSystem(:plant => PopulationComponent(plant, n0), :shrub => PopulationComponent(shrub, [5.0, 1.0]);
                               state = (seen = 0.0,))
        recording = CustomSubstep(:recording, (s, day, p) -> (crowding = crowding(s),))
        browsing = CustomRule(:browsing, (s, day, p) -> [StageFlux(:shrub, 2, -0.1 * s[:shrub].population[2])])
        total = Observable(:total, (s, t, p) -> total_population(s))
        args = (sys, (0, 40))
        kw = (p = 1.5, substeps = [recording], rules = [browsing], observables = [total],
              events = [SingleDayRelease(:shrub, 1.0, 10)])
        return target === :ipm ? CoupledIPMProblem(args...; kw...) : MatrixProjectionModels.CoupledMPMProblem(args...; kw...)
    end

    @testset "two populations: the same as a CoupledMPMProblem" begin
        a = solve(coupled(:ipm), DirectIteration())
        b = MatrixProjectionModels.solve(coupled(:mpm), MatrixProjectionModels.DirectIteration())
        @test a.u == b.u
        @test a.component_totals == b.component_totals
        @test a.substep_log == b.substep_log && a.rule_log == b.rule_log && a.event_log == b.event_log
        @test a.observables == b.observables
        @test a.component_matrices == b.component_matrices
    end

    @testset "inflows" begin
        sys = PopulationSystem(:plant => PopulationComponent(kernel, n0))
        rain = zeros(30); rain[1:3] .= [1.0, 2.0, 0.5]
        sol = solve(CoupledIPMProblem(sys, (0, 5); inflows = Dict(:plant => (s, day, p) -> rain)))
        K = materialize(kernel)
        n = n0
        for k in 1:5
            n = K * n .+ rain
            @test sol.u[k + 1][:plant] == n
        end
    end

    @testset "whole individuals" begin
        start = round.(n0)
        rain = zeros(30); rain[2] = 1.5
        whole = DiscretizedKernel(P, F, (plant = domain,))
        sys = PopulationSystem(:plant => PopulationComponent(whole, start), :shrub => PopulationComponent(
            DiscretizedKernel(Us, Fs, stages), [5.0, 1.0]))
        inflows = Dict(:plant => (s, day, p) -> rain)
        prob = CoupledIPMProblem(sys, (0, 8); inflows)
        sol = solve(prob; rng = Xoshiro(7))
        # by hand: for each population in turn, survival and births, then its inflows
        rng = Xoshiro(7)
        plant, shrub = start, [5.0, 1.0]
        for k in 1:8
            next = Float64.(demographic_step!(rng, zeros(Int, 30), round.(Int, plant), P, F))
            for (i, r) in enumerate(rain)
                r > 0 && (next[i] += rand_poisson(rng, r))
            end
            shrub = Float64.(demographic_step!(rng, zeros(Int, 2), round.(Int, shrub), Us, Fs))
            plant = next
            @test sol.u[k + 1][:plant] == plant && sol.u[k + 1][:shrub] == shrub
        end
        @test all(all(isinteger, u[:plant]) for u in sol.u)
        # Demographic() without an rng uses the default generator, and the draws differ with the seed
        @test solve(remake(prob; stochasticity = Demographic())).u[1] == sol.u[1]
        @test solve(prob; rng = Xoshiro(8)).u != sol.u
        @test sol.component_matrices[:plant][1] == P .+ F
        # the demographic mode needs whole numbers, survival and births apart, and survival probabilities
        @test_throws ArgumentError solve(CoupledIPMProblem(PopulationSystem(:plant => PopulationComponent(whole, n0)), (0, 2));
                                         rng = Xoshiro(1))
        @test_throws ArgumentError solve(CoupledIPMProblem(PopulationSystem(:plant => PopulationComponent(P .+ F, start)), (0, 2));
                                         rng = Xoshiro(1))
        @test_throws ArgumentError solve(CoupledIPMProblem(PopulationSystem(:plant => PopulationComponent(
            DiscretizedKernel(2 .* P, F, (plant = domain,)), start)), (0, 2)); rng = Xoshiro(1))
        # PKernel and FKernel sub-kernels keep them apart too
        @test all(isinteger, solve(CoupledIPMProblem(PopulationSystem(:plant => PopulationComponent(kernel, start)), (0, 3));
                                   rng = Xoshiro(1)).u[end][:plant])
    end
end
