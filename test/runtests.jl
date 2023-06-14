using SafeTestsets, Test

const GROUP = get(ENV, "GROUP", "All")

const is_APPVEYOR = Sys.iswindows() && haskey(ENV, "APPVEYOR")

const is_TRAVIS = haskey(ENV, "TRAVIS")

const is_CI = haskey(ENV, "CI")

@time begin
#    # Currently verified by Downstream tests
#    @test true
    if GROUP == "All" || GROUP == "VariableMap"
        @time @safetestset "VariableMap" begin include("variable_map_tests.jl") end
    end
end
