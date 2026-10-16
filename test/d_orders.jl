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

    @testset "Num-wrapped Pair boundary conditions" begin
        @variables w(..)
        Dt = Differential(t)
        bcs = [w(0, x) => sin(x), Dt(w(0, x)) => 0, Dx(w(t, 0)) => 0, w(t, 1) ~ 0]
        @test d_orders(x, bcs) == [1]
        @test d_orders(t, bcs) == [1]
        @test d_orders(x, [Dx(Dx(w(t, x))) => Dx(w(t, x))]) == [2, 1]
    end
end
