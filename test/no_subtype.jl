using PDEBase
using PDEBase.StrongForm: analyse
using ModelingToolkit
using ModelingToolkit: get_bcs
using Symbolics: unwrap
using SymbolicUtils: operation
using Test
import DomainSets
using DomainSets: Interval, Rectangle, (..)

# Not a subtype of any PDEBase or SciMLBase type: `VariableMap` only asks it for `get_time`.
struct PlainTime
    time::Any
end
PDEBase.get_time(p::PlainTime) = p.time

@parameters t x y
@variables u(..)
Dt, Dx, Dy = Differential(t), Differential(x), Differential(y)
rect = Rectangle(0.0 .. 1.0, 0.0 .. 2.0)
Γ_left, Γ_right = boundarypiece(rect, :left), boundarypiece(rect, :right)

@testset "PDEBase functions without a discretization type" begin
    @named heat = PDESystem(
        [Dt(u(t, x, y)) ~ Dx(Dx(u(t, x, y))) + Dy(Dy(u(t, x, y)))],
        [u(0, x, y) ~ 0, u(t, Γ_left) ~ sin(t), NormalDerivative(Γ_right)(u(t, Γ_right)) ~ 0],
        [(x, y) ∈ rect, t ∈ Interval(0.0, 1.0)], [t, x, y], [u(t, x, y)]
    )
    lowered = to_box_bcs(heat; time = t)
    bcs = get_bcs(lowered)
    @test isequal(bcs[2], u(t, 0, y) ~ sin(t))
    @test isequal(bcs[3], Dx(u(t, 1.0, y)) ~ 0)

    v = VariableMap(lowered, PlainTime(t))
    @test isequal(v.time, t)
    orders = Dict(map(iv -> iv => d_orders(iv, bcs, v), all_ivs(v)))
    boundarymap = PDEBase.parse_bcs(bcs, v, orders)
    @test length(boundarymap[operation(unwrap(u(t, x, y)))][x]) == 2
    @test !isempty(PDEBase.flatten_vardict(boundarymap))

    analysis = analyse(lowered, v, boundarymap)
    @test (analysis.nfields, analysis.dim, length(analysis.strong)) == (1, 2, 1)
end
