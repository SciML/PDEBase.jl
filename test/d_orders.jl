using PDEBase
using ModelingToolkit
using Test

@testset "d_orders" begin
    @parameters x t
    @variables u(..)
    Dx = Differential(x)

    @test d_orders(x, [Dx(Dx(u(x))) ~ 0, Dx(u(x)) ~ 0]) == [2, 1]
    @test d_orders(x, [Dx(Dx(u(x))) ~ 0]) == [2]
    @test d_orders(t, [Dx(Dx(u(x, t))) ~ 0]) == Int[]

    # A system without boundary conditions has no equations to inspect.
    @testset "empty equation collection" for eqs in (Equation[], Pair[], ())
        orders = d_orders(x, eqs)
        @test orders isa Vector{Int}
        @test isempty(orders)
        @test collect(union(d_orders(x, [Dx(Dx(u(x))) ~ 0]), orders)) == [2]
    end
end
