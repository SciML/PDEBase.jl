using PDEBase
using PDEBase: DomainSets
using PDEBase.DomainSets: (..)
using ModelingToolkit
using Symbolics
using SymbolicUtils: unwrap_const
using Test
using Symbolics: unwrap

# The fallbacks lack the parts of the upstream definitions that need changes inside
# Symbolics; those checks are known broken until Symbolics defines `NormalDerivative`.
const UPSTREAM_DN = isdefined(Symbolics, :NormalDerivative)

struct TestPiece <: DomainSets.Domain{Float64}
    name::Symbol
end
Base.show(io::IO, p::TestPiece) = print(io, "Γ_", p.name)

@testset "upstream gating" begin
    @test PDEBase.NormalDerivative ===
        (UPSTREAM_DN ? Symbolics.NormalDerivative : PDEBase.NormalDerivative)
    @test UPSTREAM_DN || parentmodule(NormalDerivative) === PDEBase
    @test PDEBase.BoundaryPiece === (
        isdefined(DomainSets, :BoundaryPiece) ? DomainSets.BoundaryPiece :
            PDEBase.BoundaryPiece
    )
    @test isdefined(DomainSets, :BoundaryPiece) || parentmodule(BoundaryPiece) === PDEBase
end

# Ported from the "NormalDerivative" testset of JuliaSymbolics/Symbolics.jl#2115.
@testset "NormalDerivative" begin
    @variables t x u(..) v(..)
    @variables w(..)[1:2]
    Γ = TestPiece(:wall)
    Dn = NormalDerivative(Γ)

    @testset "construction, printing and recovery" begin
        tr = u(t, Γ)
        ex = Dn(tr)
        @test ex isa Num
        @test operation(unwrap(ex)) isa NormalDerivative
        @test isequal(only(arguments(unwrap(ex))), unwrap(tr))
        @test repr(ex) == "Dn_wall(u(t, Γ_wall))"
        @test repr(Dn) == "Dn_wall"
        @test repr(NormalDerivative(DomainSets.UnitInterval())) ==
            "Dn_" * repr(DomainSets.UnitInterval())

        @test operation(unwrap(ex)).piece === Γ
        @test unwrap_const(last(arguments(only(arguments(unwrap(ex)))))) === Γ

        arr = Dn(w(t, Γ))
        @test arr isa Symbolics.Arr
        @test operation(unwrap(arr)) isa NormalDerivative
        @test size(arr) == (2,)
    end

    @testset "equality and hashing" begin
        Dn2 = NormalDerivative(TestPiece(:wall))
        @test isequal(Dn, Dn2)
        @test hash(Dn) == hash(Dn2)
        @test !isequal(Dn, NormalDerivative(TestPiece(:other)))
        a = Dn(u(t, Γ))
        b = Dn2(u(t, TestPiece(:wall)))
        @test isequal(a, b)
        @test hash(a) == hash(b)
        @test !isequal(a, NormalDerivative(TestPiece(:other))(u(t, Γ)))
        @test isequal(2 * a + a, 3 * b)
    end

    @testset "expand_derivatives leaves it untouched" begin
        ex = Dn(u(t, Γ))
        @test isequal(expand_derivatives(ex), ex)
        ex2 = 2 * Dn(u(t, Γ)) + u(t, Γ)
        @test isequal(expand_derivatives(ex2), ex2)
        Dt = Differential(t)
        @test isequal(expand_derivatives(Dt(ex)), Dt(ex)) broken = !UPSTREAM_DN
        @test isequal(expand_derivatives(Dt(2 * ex)), 2 * Dt(ex)) broken = !UPSTREAM_DN
    end

    @testset "mismatched pieces are not an error here" begin
        ex = NormalDerivative(TestPiece(:a))(u(t, TestPiece(:b)))
        @test operation(unwrap(ex)).piece == TestPiece(:a)
        @test unwrap_const(last(arguments(only(arguments(unwrap(ex)))))) == TestPiece(:b)
    end

    @testset "substitution through operators" begin
        Dt = Differential(t)
        rule = u(t, Γ) => u(t, x)
        dn = Dn(u(t, Γ))
        @test isequal(substitute(dn, rule), dn)
        # Descent into `Dn(...)` needs Symbolics' substitution filter to accept any operator.
        @test isequal(Symbolics.substitute_in_deriv(dn, rule), Dn(u(t, x))) broken = !UPSTREAM_DN
        @test isequal(Symbolics.substitute_in_deriv(2 * dn + 1, rule), 2 * Dn(u(t, x)) + 1) broken = !UPSTREAM_DN
        dd = Dt(u(t, x))
        @test isequal(substitute(dd, u(t, x) => v(t, x)), dd)
        @test isequal(Symbolics.substitute_in_deriv(dd, u(t, x) => v(t, x)), Dt(v(t, x)))
        Dn_other = NormalDerivative(TestPiece(:other))
        @test isequal(substitute(Dn, Dict(Dn => Dn_other)), Dn_other)
        @test isequal(Symbolics.substitute_in_deriv(Dn(u(t, Γ)), Dict(Dn => Dn_other)), Dn_other(u(t, Γ))) broken = !UPSTREAM_DN
        nested = Dt(Dn(u(t, Γ)))
        @test isequal(substitute(nested, rule), nested)
        @test isequal(Symbolics.substitute_in_deriv(nested, rule), Dt(Dn(u(t, x)))) broken = !UPSTREAM_DN
    end
end

@testset "BoundaryPiece" begin
    rect = DomainSets.Rectangle(0.0 .. 1.0, 2.0 .. 3.0)

    @testset "construction, show and parent" begin
        p = BoundaryPiece(rect, :inlet)
        @test p isa DomainSets.Domain{eltype(rect)}
        @test parent(p) === rect
        @test p.name === :inlet
        @test p.geometry === nothing
        @test repr(p) == "Γ_inlet"
        @test isequal(p, BoundaryPiece(DomainSets.Rectangle(0.0 .. 1.0, 2.0 .. 3.0), :inlet))
        @test hash(p) == hash(BoundaryPiece(DomainSets.Rectangle(0.0 .. 1.0, 2.0 .. 3.0), :inlet))
        @test !isequal(p, BoundaryPiece(rect, :outlet))
        @test !isequal(p, BoundaryPiece(DomainSets.Rectangle(0.0 .. 2.0, 2.0 .. 3.0), :inlet))
    end

    @testset "rectangle faces" begin
        pieces = boundarypieces(rect)
        faces = DomainSets.components(DomainSets.boundary(rect))
        @test [p.name for p in pieces] == [:bottom, :right, :top, :left]
        @test all(p -> parent(p) === rect, pieces)
        @test all(isequal(p.geometry, f) for (p, f) in zip(pieces, faces))
        @test isequal(boundarypiece(rect, :top), pieces[3])
        @test [PDEBase._box_face(p) for p in pieces] == [(2, :lo), (1, :hi), (2, :hi), (1, :lo)]
        # each face lies on the side of the axis its name says
        for p in pieces
            axis, side = PDEBase._box_face(p)
            bound = side === :lo ? DomainSets.leftendpoint(rect)[axis] :
                DomainSets.rightendpoint(rect)[axis]
            @test DomainSets.forward_map(p.geometry, 0.5)[axis] ≈ bound
        end
        err = try
            boundarypiece(rect, :front)
            nothing
        catch e
            e
        end
        @test err isa ArgumentError
        @test occursin(":front", err.msg)
        @test all(n -> occursin(repr(n), err.msg), (:bottom, :right, :top, :left))
    end

    @testset "3D box faces" begin
        box = DomainSets.Rectangle(0.0 .. 1.0, 0.0 .. 2.0, 0.0 .. 3.0)
        pieces = boundarypieces(box)
        faces = DomainSets.components(DomainSets.boundary(box))
        @test [p.name for p in pieces] == [:x1_lo, :x1_hi, :x2_lo, :x2_hi, :x3_lo, :x3_hi]
        @test all(isequal(p.geometry, f) for (p, f) in zip(pieces, faces))
        @test [PDEBase._box_face(p) for p in pieces] ==
            [(1, :lo), (1, :hi), (2, :lo), (2, :hi), (3, :lo), (3, :hi)]
        for p in pieces
            axis, side = PDEBase._box_face(p)
            bound = side === :lo ? DomainSets.leftendpoint(box)[axis] :
                DomainSets.rightendpoint(box)[axis]
            @test DomainSets.forward_map(p.geometry, [0.5, 0.5])[axis] ≈ bound
        end
        @test_throws ArgumentError boundarypiece(box, :bottom)
        @test_throws ArgumentError PDEBase._box_face(BoundaryPiece(box, :inlet))
    end

    @testset "trace syntax" begin
        @parameters t x y
        @variables u(..)
        Γ = boundarypiece(rect, :left)
        Dn = NormalDerivative(Γ)
        tr = u(t, Γ)
        @test repr(tr) == "u(t, Γ_left)"
        @test isequal(tr, u(t, boundarypiece(DomainSets.Rectangle(0.0 .. 1.0, 2.0 .. 3.0), :left)))
        @test repr(Dn(tr)) == "Dn_left(u(t, Γ_left))"

        Dt = Differential(t)
        Dxx = Differential(x)^2
        eqs = [Dt(u(t, x, y)) ~ Dxx(u(t, x, y))]
        bcs = [u(t, Γ) ~ 0, Dn(u(t, boundarypiece(rect, :right))) ~ 1]
        domains = [(x, y) ∈ rect, t ∈ DomainSets.Interval(0.0, 1.0)]
        @named pdesys = PDESystem(eqs, bcs, domains, [t, x, y], [u(t, x, y)])
        u_op = [operation(unwrap(u(t, x, y)))]
        @test isequal(PDEBase.get_all_depvars(pdesys, u_op), [unwrap(u(t, x, y))])
        # Over equations and conditions together a trace still counts as a dependent
        # variable; the trace-aware VariableMap must exclude it.
        @test isequal(PDEBase.get_all_depvars([eqs; bcs], u_op), [unwrap(u(t, x, y))]) broken = true
    end
end
