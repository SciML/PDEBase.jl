using SafeTestsets, Test

const GROUP = get(ENV, "GROUP", "All")

const is_APPVEYOR = Sys.iswindows() && haskey(ENV, "APPVEYOR")

const is_TRAVIS = haskey(ENV, "TRAVIS")

const is_CI = haskey(ENV, "CI")

if GROUP == "All" || GROUP == "Core"
    # Currently verified by Downstream tests
    @test true

    @safetestset "Allocation Tests" begin
        include("alloc_tests.jl")
    end
end

if GROUP == "All" || GROUP == "QA"
    using Pkg
    Pkg.activate(joinpath(@__DIR__, "qa"))
    Pkg.develop(path = dirname(@__DIR__))
    Pkg.instantiate()

    @safetestset "JET Static Analysis" begin
        include("qa/jet_tests.jl")
    end
end
