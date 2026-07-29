using PDEBase
using ModelingToolkit
using Symbolics
using SymbolicUtils
using Test

@testset "Flattening array dependent variables" begin
    N = 3
    @parameters x
    @variables (u(..))[1:N]

    dvs = collect([u(x)[i] for i in 1:N])
    ch = PDEBase.chain_flatten_array_variables(dvs)
    flattened = map(dv -> ch(PDEBase.safe_unwrap(dv)), dvs)

    # Each `u[i]` needs its own replacement; keying only on `u` collapses them all.
    @test length(unique(string.(flattened))) == N

    for f in flattened
        # The replacement stays a call in the independent variables so downstream
        # `operation`/`arguments` keep working. An `Arr` here means the index was
        # never resolved away.
        @test !(f isa Symbolics.Arr)
        @test iscall(f)
        @test isequal(only(arguments(f)), PDEBase.safe_unwrap(x))
    end

    @test PDEBase.chain_flatten_array_variables([]) === identity
end
