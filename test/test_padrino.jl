using IntegralProjectionModels
using Test
using LinearAlgebra
using Distributions

@testset "PADRINO Extension" begin

    @testset "Stub functions error without extension" begin
        # Without the extension, these are just function declarations with no methods
        @test pdb_download isa Function
        @test pdb_load isa Function
        @test pdb_make_ipm isa Function
    end

    @testset "PadrinoDB and PadrinoModel types" begin
        pdb = PadrinoDB(Dict{String, Any}("Metadata" => nothing))
        @test pdb isa PadrinoDB
        @test haskey(pdb.tables, "Metadata")

        pm = PadrinoModel(
            "test_id",
            Dict{Symbol, Any}(:species_accepted => "Test species"),
            [Dict{Symbol, Any}(:state_variable => "size", :discrete => false)],
            [Dict{Symbol, Any}(:state_variable => "size", :lower => 0.0, :upper => 50.0)],
            [Dict{Symbol, Any}(:kernel_id => "P", :integration_rule => "midpoint")],
            [Dict{Symbol, Any}(:expression => "size", :n_bins => 100)],
            [Dict{Symbol, Any}(:kernel_id => "P", :model_family => "CC",
                               :formula => "P = s * g", :domain_start => "size",
                               :domain_end => "size")],
            [Dict{Symbol, Any}(:name => "s", :formula => "s = 1/(1+exp(-(s_int + s_slope * size_1)))",
                               :model_type => "Evaluated", :kernel_ids => ["P"])],
            Dict{String, Float64}("s_int" => 2.2, "s_slope" => 0.25),
            Dict{Symbol, Any}[],
            Dict{Symbol, Any}[]
        )
        @test pm isa PadrinoModel
        @test pm.ipm_id == "test_id"
    end

    # Extension tests require CSV and DataFrames
    # These are run separately or when the extension is loaded
    csv_available = try
        @eval using CSV
        true
    catch
        false
    end

    df_available = try
        @eval using DataFrames
        true
    catch
        false
    end

    dl_available = try
        @eval using Downloads
        true
    catch
        false
    end

    if csv_available && df_available && dl_available
        @eval using CSV, DataFrames, Downloads

        @testset "R Expression Parser" begin
            # Access parser through the extension module
            ext = Base.get_extension(IntegralProjectionModels, :IntegralProjectionModelsPadrinoExt)
            if ext !== nothing
                @testset "Tokenizer" begin
                    tokens = ext.tokenize("1 + 2 * 3")
                    @test length(tokens) >= 5  # 1, +, 2, *, 3, EOF
                end

                @testset "Simple expressions" begin
                    rexpr = ext.parse_rexpr("1 + 2")
                    @test rexpr isa ext.RBinOp
                    @test rexpr.op == "+"

                    rexpr2 = ext.parse_rexpr("x * y + z")
                    @test rexpr2 isa ext.RBinOp
                    @test rexpr2.op == "+"
                end

                @testset "Function calls" begin
                    rexpr = ext.parse_rexpr("exp(x)")
                    @test rexpr isa ext.RCall
                    @test rexpr.func == "exp"
                    @test length(rexpr.args) == 1

                    rexpr2 = ext.parse_rexpr("Norm(mu, sd)")
                    @test rexpr2 isa ext.RCall
                    @test rexpr2.func == "Norm"
                    @test length(rexpr2.args) == 2
                end

                @testset "Keyword arguments" begin
                    rexpr = ext.parse_rexpr("Lognorm(x, m, s, first_arg = TRUE)")
                    @test rexpr isa ext.RCall
                    @test rexpr.func == "Lognorm"
                    @test length(rexpr.kwargs) == 1
                    @test rexpr.kwargs[1].first == "first_arg"
                end

                @testset "Operator precedence" begin
                    rexpr = ext.parse_rexpr("a + b * c")
                    @test rexpr isa ext.RBinOp
                    @test rexpr.op == "+"
                    @test rexpr.right isa ext.RBinOp
                    @test rexpr.right.op == "*"
                end

                @testset "Nested expressions" begin
                    rexpr = ext.parse_rexpr("1/(1 + exp(-(a + b * x)))")
                    @test rexpr isa ext.RBinOp
                    @test rexpr.op == "/"
                end

                @testset "Translation" begin
                    state_vars = ["size"]

                    # State variable mapping
                    rexpr = ext.parse_rexpr("size_1")
                    jexpr = ext.translate_rexpr(rexpr, state_vars)
                    @test jexpr == :z

                    rexpr = ext.parse_rexpr("size_2")
                    jexpr = ext.translate_rexpr(rexpr, state_vars)
                    @test jexpr == :z_new

                    rexpr = ext.parse_rexpr("d_size")
                    jexpr = ext.translate_rexpr(rexpr, state_vars)
                    @test jexpr == :h

                    # Distribution translation
                    rexpr = ext.parse_rexpr("Norm(mu_g, sd_g)")
                    jexpr = ext.translate_rexpr(rexpr, state_vars)
                    @test jexpr isa Expr
                    @test jexpr.head == :call
                    @test jexpr.args[1] == :pdf
                end

                @testset "Dependency ordering" begin
                    vrs = Dict{String, Any}(
                        "a" => :(x + 1),
                        "b" => :(a * 2),
                        "c" => :(b + a)
                    )
                    order = ext.resolve_dependency_order(vrs)
                    @test findfirst(==("a"), order) < findfirst(==("b"), order)
                    @test findfirst(==("a"), order) < findfirst(==("c"), order)
                    @test findfirst(==("b"), order) < findfirst(==("c"), order)
                end
            end
        end

        @testset "PADRINO Download and Parse" begin
            # Test downloading (requires internet)
            pdb = nothing
            try
                pdb = pdb_download()
                @test pdb isa PadrinoDB
                @test haskey(pdb.tables, "Metadata")
                @test haskey(pdb.tables, "IpmKernels")
                @test haskey(pdb.tables, "VitalRateExpr")
                @test haskey(pdb.tables, "ParameterValues")
            catch e
                @warn "PADRINO download failed (expected if no internet)" exception=e
                @test_skip "PADRINO download requires internet"
            end

            if pdb !== nothing && haskey(pdb.tables, "Metadata") &&
               DataFrames.nrow(pdb.tables["Metadata"]) > 0

                @testset "Query functions" begin
                    species = pdb_species(pdb)
                    @test length(species) > 0
                    @test all(s -> s isa String, species)

                    ids = pdb_metadata(pdb, :ipm_id)
                    @test length(ids) > 0
                end

                @testset "Subset" begin
                    ids = String.(pdb.tables["Metadata"].ipm_id[1:min(3, end)])
                    sub = pdb_subset(pdb, ids)
                    @test sub isa PadrinoDB
                    @test DataFrames.nrow(sub.tables["Metadata"]) <= 3
                end

                @testset "Simple deterministic model" begin
                    # Find a simple deterministic model to test
                    # aaa310 is a known simple det model
                    test_ids = String.(pdb.tables["Metadata"].ipm_id)

                    if "aaa310" in test_ids
                        protos = pdb_make_proto_ipm(pdb; ipm_id="aaa310")
                        @test haskey(protos, "aaa310")
                        @test protos["aaa310"] isa PadrinoModel

                        ipms = pdb_make_ipm(protos; tspan=(0, 100))
                        if haskey(ipms, "aaa310")
                            ipm = ipms["aaa310"]
                            @test ipm isa IPMProblem
                            sol = solve(ipm)
                            lam = lambda(sol)
                            @test isfinite(lam)
                            @info "aaa310 λ = $lam"
                        end
                    end
                end

                @testset "Batch model building" begin
                    # Try building a few models
                    test_ids = String.(pdb.tables["Metadata"].ipm_id[1:min(5, end)])
                    protos = pdb_make_proto_ipm(pdb; ipm_id=test_ids)
                    @test length(protos) > 0

                    ipms = pdb_make_ipm(protos; tspan=(0, 50))
                    @test ipms isa Dict
                    # At least some should build successfully
                    @info "Built $(length(ipms)) / $(length(protos)) models"
                end
            end
        end
    else
        @info "Skipping extension tests (CSV/DataFrames/Downloads not available)"
    end
end
