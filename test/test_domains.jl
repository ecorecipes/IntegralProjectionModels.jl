@testset "Domains" begin
    @testset "ContinuousDomain" begin
        d = ContinuousDomain(0.0, 50.0, 100)
        @test d.lower == 0.0
        @test d.upper == 50.0
        @test d.n_meshpoints == 100
        @test step_size(d) ≈ 0.5
        @test n_states(d) == 100

        z = meshpoints(d)
        @test length(z) == 100
        @test z[1] ≈ 0.25   # midpoint of first bin [0, 0.5]
        @test z[end] ≈ 49.75 # midpoint of last bin [49.5, 50]

        b = bounds(d)
        @test length(b) == 101
        @test b[1] ≈ 0.0
        @test b[end] ≈ 50.0

        # Midpoint is average of consecutive bounds
        for i in 1:100
            @test z[i] ≈ (b[i] + b[i + 1]) / 2
        end
    end

    @testset "ContinuousDomain validation" begin
        @test_throws ArgumentError ContinuousDomain(50.0, 0.0, 100)
        @test_throws ArgumentError ContinuousDomain(0.0, 50.0, 0)
        @test_throws ArgumentError ContinuousDomain(0.0, 50.0, -1)
    end

    @testset "ContinuousDomain type promotion" begin
        d = ContinuousDomain(0, 50.0, 100)
        @test d.lower isa Float64
        @test d.upper isa Float64
    end

    @testset "DiscreteDomain" begin
        d = DiscreteDomain([:seedbank, :dormant, :active])
        @test n_states(d) == 3
        @test d.labels == [:seedbank, :dormant, :active]
    end
end
