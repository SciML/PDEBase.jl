using PDEBase
using PDEBase: ArgumentClass, Free, BoxFace, Piece, Interior, ConditionClass, location, boxface, normal
using ModelingToolkit
using Symbolics: unwrap
using Test
import DomainSets
using DomainSets: Interval, Rectangle, (..)

struct ClassifyDiscretization <: PDEBase.AbstractEquationSystemDiscretization
    time::Any
end
PDEBase.get_time(disc::ClassifyDiscretization) = disc.time

struct ClassifyMesh end
PDEBase.meshdim(::ClassifyMesh) = 2
PDEBase.facetset_names(::ClassifyMesh) = (:inlet, :wall)

function errmsg(f)
    try
        f()
    catch e
        e isa ArgumentError || rethrow()
        return e.msg
    end
    return nothing
end

@parameters t x y
@variables u(..) w(..)
Dt, Dx, Dy = Differential(t), Differential(x), Differential(y)
rect = Rectangle(0.0 .. 1.0, 0.0 .. 2.0)
heat = Dt(u(t, x, y)) ~ Dx(Dx(u(t, x, y))) + Dy(Dy(u(t, x, y)))
@named rectsys = PDESystem(
    [heat], [u(0, x, y) ~ 0], [(x, y) ∈ rect, t ∈ Interval(0.0, 1.0)], [t, x, y], [u(t, x, y)]
)
v = VariableMap(rectsys, ClassifyDiscretization(t))
uop = operation(unwrap(u(t, x, y)))
Γ_left, Γ_right, Γ_bottom, Γ_top = (boundarypiece(rect, n) for n in (:left, :right, :bottom, :top))

loc(call) = location(v, classify_call(v, call)[2])

@testset "classify_argument" begin
    @test classify_argument(v, uop, 1, t) == Free(unwrap(t))
    @test classify_argument(v, uop, 2, 0) == BoxFace(unwrap(x), :lo)
    @test classify_argument(v, uop, 3, 2) == BoxFace(unwrap(y), :hi)
    @test classify_argument(v, uop, 2, Γ_left) == Piece(Γ_left)
    # The parent is compared structurally.
    @test classify_argument(v, uop, 2, boundarypiece(Rectangle(0.0 .. 1.0, 0.0 .. 2.0), :left)) == Piece(Γ_left)
    msg = errmsg(() -> classify_argument(v, uop, 2, boundarypiece(Rectangle(0.0 .. 3.0, 0.0 .. 2.0), :left)))
    @test occursin("not of the spatial domain", msg)
    msg = errmsg(() -> classify_argument(v, uop, 2, y))
    @test occursin("declared argument list is u(t, x, y)", msg)

    @testset "endpoint tolerance" begin
        @variables p(..)
        for (domain, val, side) in (
                (Interval(0.0, π), π, :hi), (Interval(0.0, 0.3), 0.1 + 0.2, :hi),
                (Interval(0.0, 0.3), 1.0e-17, :lo), (Interval(0, 1), 1, :hi),
            )
            @named sys = PDESystem(
                [Dt(p(t, x)) ~ Dx(Dx(p(t, x)))], [p(0, x) ~ 0], [t ∈ Interval(0.0, 1.0), x ∈ domain],
                [t, x], [p(t, x)]
            )
            pv = VariableMap(sys, ClassifyDiscretization(t))
            _, classes = classify_call(pv, p(t, val))
            @test classes == [Free(unwrap(t)), BoxFace(unwrap(x), side)]
        end
        msg = errmsg(() -> classify_argument(v, uop, 2, 0.25))
        @test occursin("not on the boundary", msg)
        @test occursin("nearest endpoint is the lower endpoint 0.0", msg)
        @test occursin("at distance 0.25", msg)
        msg = errmsg(() -> classify_argument(v, uop, 3, 1.5))
        @test occursin("upper endpoint 2.0", msg)
    end
end

@testset "classify_call and location" begin
    @test loc(u(t, x, y)) === Interior()
    @test isequal(loc(u(t, 0, y)), BoxFace(unwrap(x), :lo))
    @test isequal(loc(u(0, x, y)), BoxFace(unwrap(t), :lo))
    # M0 exit criterion: both spellings of a box face give one location.
    for (endpoint, piece) in (
            (u(t, 0, y), Γ_left), (u(t, 1, y), Γ_right),
            (u(t, x, 0), Γ_bottom), (u(t, x, 2), Γ_top),
        )
        @test isequal(loc(endpoint), loc(u(t, piece)))
        @test hash(loc(endpoint)) == hash(loc(u(t, piece)))
    end
    @test boxface(v, Γ_left) == BoxFace(unwrap(x), :lo)
    op, classes = classify_call(v, u(t, Γ_left))
    @test isequal(op, uop)
    @test classes == [Free(unwrap(t)), Piece(Γ_left)]
    @test location(classes) == Piece(Γ_left)

    msg = errmsg(() -> classify_call(v, u(Γ_left)))
    @test occursin("trace is missing the time argument; write `u(t, Γ_left)`", msg)
    @test occursin("has 2 arguments", errmsg(() -> classify_call(v, u(t, x))))
    @test occursin("declared argument list", errmsg(() -> classify_call(v, u(x, t, y))))
    @test occursin("declared positions", errmsg(() -> classify_call(v, u(Γ_left, t))))
    @test occursin("edge or a corner", errmsg(() -> loc(u(t, 0, 0))))
    @test occursin("edge or a corner", errmsg(() -> location([BoxFace(unwrap(x), :lo), Piece(Γ_left)])))

    @test normal(BoxFace(unwrap(x), :lo)) == (unwrap(x), -1)
    @test normal(BoxFace(unwrap(y), :hi)) == (unwrap(y), 1)
    @test occursin("backend supplies it", errmsg(() -> normal(Piece(Γ_left), [0.0, 1.0])))
end

@testset "classify_condition" begin
    Dn_left, Dn_right = NormalDerivative(Γ_left), NormalDerivative(Γ_right)
    class(eq) = classify_condition(v, eq)
    @test class(u(0, x, y) ~ sin(x)) isa ConditionClass
    @test class(u(0, x, y) ~ sin(x)).kind === :initial
    @test class(Dt(u(0, x, y)) => 0).kind === :initial
    @test class(u(0, x, y) ~ sin(x)).locations == [BoxFace(unwrap(t), :lo)]
    @test occursin("not an initial condition", errmsg(() -> class(u(1, x, y) ~ 0)))
    @test occursin("not an initial condition", errmsg(() -> class(u(1.0, x, y) ~ u(0, x, y))))
    for bc in (
            u(t, 0, y) ~ 0, Dx(u(t, 1, y)) ~ 0, u(t, Γ_left) ~ sin(t), Dn_left(u(t, Γ_left)) ~ 1,
            2 * Dn_right(u(t, Γ_right)) + 3 * u(t, Γ_right) ~ 0,
        )
        @test class(bc).kind === :boundary
    end
    @test class(u(t, Γ_left) ~ 0).locations == class(u(t, 0, y) ~ 0).locations
    periodic = class(u(t, 1, y) ~ u(t, 0, y))
    @test periodic.kind === :periodic
    @test periodic.locations == [BoxFace(unwrap(x), :lo), BoxFace(unwrap(x), :hi)]
    @test class(u(t, Γ_bottom) ~ u(t, Γ_top)).kind === :periodic
    @test class(u(t, Γ_bottom) ~ u(t, x, 2)).kind === :periodic
    @test occursin("relates the locations", errmsg(() -> class(u(t, 0, y) ~ u(t, x, 0))))
    @test occursin("mixes interior calls with traces", errmsg(() -> class(u(t, 0, y) ~ u(t, x, y))))
    @test occursin("every call is in the interior", errmsg(() -> class(u(t, x, y) ~ 0)))
    @test occursin("edge or a corner", errmsg(() -> class(u(t, 0, 2) ~ 0)))
    # §3.3 (b): the operator's piece must be the trace's piece.
    msg = errmsg(() -> class(Dn_left(u(t, Γ_right)) ~ 0))
    @test occursin("not on Γ_left", msg)
    @test occursin("not on Γ_left", errmsg(() -> class(Dn_left(u(t, 1, y)) ~ 0)))
    @test class(Dn_left(u(t, 0, y)) ~ 0).kind === :boundary
    @test occursin("missing the time argument", errmsg(() -> class(u(Γ_left) ~ 0)))

    @testset "Brusselator's time-last declaration" begin
        @variables b(..)
        eq = Dt(b(x, y, t)) ~ Dx(Dx(b(x, y, t))) + Dy(Dy(b(x, y, t)))
        bcs = [b(x, y, 0) ~ 0, b(0, y, t) ~ b(1, y, t), b(x, 0, t) ~ b(x, 1, t)]
        @named bruss = PDESystem(
            [eq], bcs, [x ∈ Interval(0.0, 1.0), y ∈ Interval(0.0, 1.0), t ∈ Interval(0.0, 11.5)],
            [x, y, t], [b(x, y, t)]
        )
        bv = VariableMap(bruss, ClassifyDiscretization(t))
        @test [classify_condition(bv, bc).kind for bc in bcs] == [:initial, :periodic, :periodic]
        @test occursin("not an initial condition", errmsg(() -> classify_condition(bv, b(x, y, 11.5) ~ 0)))
        # Split interval pairings: axis i of the box is indvars(bv)[i].
        box = spatial_domain(bv)
        Γ = boundarypiece(box, :left)
        @test classify_condition(bv, b(Γ, t) ~ 0).locations == [BoxFace(unwrap(x), :lo)]
        @test occursin("declared positions", errmsg(() -> classify_call(bv, b(t, Γ))))
    end

    @testset "steady trace u(Γ)" begin
        @variables s(..)
        @named steady = PDESystem(
            [Dx(Dx(s(x, y))) + Dy(Dy(s(x, y))) ~ 0], [s(0, y) ~ 0], [(x, y) ∈ rect], [x, y], [s(x, y)]
        )
        sv = VariableMap(steady, ClassifyDiscretization(nothing))
        @test classify_condition(sv, s(Γ_top) ~ 1).locations == [BoxFace(unwrap(y), :hi)]
    end

    @testset "pieces of a mesh" begin
        @variables m(..)
        Ω = MeshDomain(ClassifyMesh())
        @named meshsys = PDESystem(
            [Dt(m(t, x, y)) ~ Dx(Dx(m(t, x, y)))], [m(0, x, y) ~ 0], [(x, y) ∈ Ω, t ∈ Interval(0.0, 1.0)],
            [t, x, y], [m(t, x, y)]
        )
        mv = VariableMap(meshsys, ClassifyDiscretization(t))
        inlet, wall = boundarypiece(Ω, :inlet), boundarypiece(Ω, :wall)
        @test boxface(mv, inlet) === nothing
        bc = classify_condition(mv, NormalDerivative(inlet)(m(t, inlet)) ~ 0)
        @test bc.kind === :boundary
        @test bc.locations == [Piece(inlet)]
        @test d_orders(inlet, [NormalDerivative(inlet)(m(t, inlet)) ~ 0]) == [1]
        @test d_orders(x, [NormalDerivative(inlet)(m(t, inlet)) ~ 0], mv) == Int[]
        @test occursin("not on Γ_wall", errmsg(() -> classify_condition(mv, NormalDerivative(wall)(m(t, inlet)) ~ 0)))
        @test occursin("number 0 in the place of x", errmsg(() -> classify_call(mv, m(t, 0, y))))
        @test occursin("relates the locations", errmsg(() -> classify_condition(mv, m(t, inlet) ~ m(t, wall))))
    end
end

@testset "substitute_traces" begin
    k, α = 2, 3
    Dn = NormalDerivative(Γ_left)
    rules = Dict(u(t, Γ_left) => u(t, 0, y))
    for (ex, expected) in (
            (Dx(u(t, Γ_left)), Dx(u(t, 0, y))), (Dn(u(t, Γ_left)), Dn(u(t, 0, y))),
            (k * Dn(u(t, Γ_left)) + α * u(t, Γ_left), k * Dn(u(t, 0, y)) + α * u(t, 0, y)),
        )
        @test isequal(substitute_traces(ex, rules), expected)
        # Plain substitute leaves the trace inside the operator.
        @test !isequal(substitute(ex, rules), expected)
    end
    eq = substitute_traces(Dn(u(t, Γ_left)) ~ sin(t), rules)
    @test eq isa Equation
    @test isequal(eq.lhs, Dn(u(t, 0, y)))

    slot = Dict(u(t, Γ_left) => u(t, x, y), Dn(u(t, Γ_left)) => -Dx(u(t, x, y)))
    robin = k * Dn(u(t, Γ_left)) + α * u(t, Γ_left) ~ 0
    replaced = substitute_traces(robin, slot)
    @test assert_no_traces(replaced) === replaced
    @test isequal(replaced.lhs, -k * Dx(u(t, x, y)) + α * u(t, x, y))
    # A deliberately skipped operator: replacing the trace alone leaves `Dn` behind.
    msg = errmsg(() -> assert_no_traces(substitute_traces(robin, Dict(u(t, Γ_left) => u(t, 0, y)))))
    @test occursin("Dn_left", msg)
    @test occursin("must be replaced", msg)
    @test occursin("u(t, Γ_left)", errmsg(() -> assert_no_traces(Dx(u(t, Γ_left)))))
    @test occursin("u(t, Γ_left)", errmsg(() -> assert_no_traces(substitute(Dx(u(t, Γ_left)), rules))))
end

@testset "d_orders with NormalDerivative" begin
    Dn = NormalDerivative(Γ_left)
    bcs = [Dn(u(t, Γ_left)) ~ 0, u(t, Γ_right) ~ 0, Dx(u(t, Γ_right)) ~ 0]
    @test d_orders(x, bcs) == [1]
    @test d_orders(x, bcs[1:2]) == Int[]
    @test d_orders(x, bcs, v) == [1]
    @test d_orders(x, bcs[1:2], v) == [1]
    @test d_orders(y, bcs[1:2], v) == Int[]
    @test d_orders(t, bcs, v) == Int[]
    @test d_orders(Γ_left, bcs) == [1]
    @test d_orders(Γ_right, bcs) == Int[]
    @test d_orders(x, [Dn(u(t, Γ_left)) + Dx(Dx(u(t, 0, y))) ~ 0], v) == [2, 1]
end
