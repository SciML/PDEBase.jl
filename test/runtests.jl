using SafeTestsets, Test
using SciMLTesting

const is_APPVEYOR = Sys.iswindows() && haskey(ENV, "APPVEYOR")

const is_TRAVIS = haskey(ENV, "TRAVIS")

const is_CI = haskey(ENV, "CI")

run_tests(;
    core = () -> begin
        # Currently verified by Downstream tests
        @test true

        @safetestset "Allocation Tests" begin
            include("alloc_tests.jl")
        end
    end,
    # `qa` runs under both "All" and "QA" (matches the original
    # `if GROUP == "All" || GROUP == "QA"` branch). The qa env has no
    # [sources] table, so develop (default) reproduces the original
    # `Pkg.develop(path = dirname(@__DIR__))` of the repo root.
    qa = (;
        env = joinpath(@__DIR__, "qa"),
        body = () -> begin
            @safetestset "JET Static Analysis" begin
                include("qa/jet_tests.jl")
            end
        end,
    ),
)
