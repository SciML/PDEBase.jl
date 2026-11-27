using PDEBase
using Test

const TTFX_CHILD = joinpath(pkgdir(PDEBase), "benchmark", "ttfx_child.jl")
# Pkg.test narrows JULIA_LOAD_PATH to the test project, which hides the standard libraries the
# import stage needs, so the child gets the default load path.
ttfx_cmd(stage) = addenv(
    `$(Base.julia_cmd()) --startup-file=no --project=$(Base.active_project()) $TTFX_CHILD $stage`,
    "JULIA_LOAD_PATH" => nothing
)

# Runs the import stage through a duplex pipe: reads `ready`, sends `trigger` (or closes stdin
# when it is `nothing`) and returns the readiness line, the reply and the exit status.
function import_stage(trigger)
    child = open(pipeline(ttfx_cmd("import"); stderr = devnull), "r+")
    ready = readline(child)
    trigger === nothing ? close(child.in) : println(child, trigger)
    reply = readline(child)
    wait(child)
    return ready, reply, success(child)
end

@testset "benchmark/ttfx_child.jl stages" begin
    mktemp() do path, io
        invalid = run(pipeline(ignorestatus(ttfx_cmd("no_such_stage")); stderr = io))
        close(io)
        @test !success(invalid)
        diagnostic = read(path, String)
        @test occursin("unknown stage \"no_such_stage\"", diagnostic)
        @test occursin("load, first_calls, import", diagnostic)
    end

    # The heat model has one field and two spatial dimensions, and only first_calls analyses it.
    @test !occursin("analysis", read(ttfx_cmd("load"), String))
    @test occursin(r"^analysis nfields=1 dim=2$"m, read(ttfx_cmd("first_calls"), String))

    ready, reply, succeeded = import_stage("go")
    @test ready == "ready"
    @test succeeded
    @test isfinite(parse(Float64, reply)) && parse(Float64, reply) > 0
    for trigger in (nothing, "not_go")
        ready, reply, succeeded = import_stage(trigger)
        @test ready == "ready"
        @test isempty(reply)
        @test !succeeded
    end
end
