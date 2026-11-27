using PDEBase
using Test

const TTFX_CHILD = joinpath(pkgdir(PDEBase), "benchmark", "ttfx_child.jl")
ttfx_cmd(stage) = `$(Base.julia_cmd()) --startup-file=no --project=$(Base.active_project()) $TTFX_CHILD $stage`

@testset "benchmark/ttfx_child.jl stages" begin
    mktemp() do path, io
        invalid = run(pipeline(ignorestatus(ttfx_cmd("frist_calls")); stderr = io))
        close(io)
        @test !success(invalid)
        diagnostic = read(path, String)
        @test occursin("unknown stage \"frist_calls\"", diagnostic)
        @test occursin("load, first_calls, import", diagnostic)
    end
    @test success(pipeline(ttfx_cmd("load")))
    @test success(pipeline(ttfx_cmd("first_calls")))
    reply = read(pipeline(ttfx_cmd("import"); stdin = IOBuffer("go\n")), String)
    @test split(reply) == ["ready", split(reply)[2]] && parse(Float64, split(reply)[2]) > 0
end
