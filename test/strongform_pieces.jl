using PDEBase, ModelingToolkit, Symbolics, SciMLBase, Test
using DomainSets: Interval, Rectangle, (..)
using PDEBase: BoxFace, PieceBoundary
using PDEBase.StrongForm: analyse, StateSlots, SlotContext, StrongFormAnalysis, EssentialBC, NaturalBC, PeriodicBC,
    boundary_conditions, slot_condition, restrict_to_face

struct PieceBoxDisc <: PDEBase.AbstractBoxDiscretization
    t::Any
end
PDEBase.get_time(d::PieceBoxDisc) = d.t

# Stops `symbolic_discretize` at the boundary map and keeps it for the analysis.
struct PieceMeshDisc <: PDEBase.AbstractMeshDiscretization
    t::Any
    captured::Base.RefValue{Any}
end
PieceMeshDisc(t) = PieceMeshDisc(t, Ref{Any}(nothing))
PDEBase.get_time(d::PieceMeshDisc) = d.t
struct StopAfterAnalysis <: Exception end
function PDEBase.check_boundarymap(bmap, v::VariableMap, d::PieceMeshDisc)
    d.captured[] = (bmap, v)
    throw(StopAfterAnalysis())
end

struct PieceMesh end
PDEBase.meshdim(::PieceMesh) = 2
PDEBase.facetset_names(::PieceMesh) = (:inlet, :wall, :outlet)

symeq(a, b) = isequal(Symbolics.simplify(Num(a) - Num(b)), 0)
errmsg(f) =
try
    f()
    nothing
catch e
    sprint(showerror, e)
end

boundarymap(bcs, v, disc) =
    PDEBase.parse_bcs(bcs, v, Dict(z => d_orders(z, bcs, v) for z in all_ivs(v)), disc)

@parameters t x y k h
@variables u(..)
Dt, Dx, Dy = Differential(t), Differential(x), Differential(y)
heat(κ) = Dt(u(t, x, y)) ~ Dx(κ * Dx(u(t, x, y))) + Dy(κ * Dy(u(t, x, y)))
Dn(Γ) = NormalDerivative(Γ)

@testset "NormalDerivative on box faces equals the Dx form" begin
    rect = Rectangle(0.0 .. 1.0, 0.0 .. 2.0)
    Γl, Γr, Γb, Γt = (boundarypiece(rect, n) for n in (:left, :right, :bottom, :top))
    domains = [(x, y) ∈ rect, t ∈ Interval(0.0, 1.0)]
    dxform = [
        u(0, x, y) ~ 0, -Dx(u(t, 0, y)) ~ h, Dx(u(t, 1, y)) ~ h,
        u(t, x, 0) ~ t * sin(x), k * Dy(u(t, x, 2)) + 2u(t, x, 2) ~ 3,
    ]
    dnform = [
        u(0, x, y) ~ 0, Dn(Γl)(u(t, Γl)) ~ h, Dn(Γr)(u(t, 1, y)) ~ h, u(t, Γb) ~ t * sin(x),
        k * Dn(Γt)(u(t, Γt)) + 2u(t, Γt) ~ 3,
    ]
    traced = [dnform[1:2]; Dn(Γr)(u(t, Γr)) ~ h; dnform[4:5]]
    system(bcs) = PDESystem(
        [heat(k)], bcs, domains, [t, x, y], [u(t, x, y)], [k => 1.0, h => 2.0]; name = :sys
    )
    function boxanalyse(bcs)
        sys = system(bcs)
        v = VariableMap(sys, PieceBoxDisc(t))
        return analyse(sys, v, boundarymap(bcs, v, PieceBoxDisc(t))), v
    end
    reference, v = boxanalyse(dxform)
    s, P = reference.slots, reference.slots.P
    ref = Dict(b.piece => b.Φ for b in reference.natural)
    an = first(boxanalyse(dnform))
    @test an isa StrongFormAnalysis{BoxFace}
    @test length(an.natural) == 3 && all(b -> symeq(b.Φ, ref[b.piece]), an.natural)
    @test only(an.essential).piece == only(reference.essential).piece
    @test symeq(only(an.essential).g, only(reference.essential).g)
    # the outward normal is -e_x on the left face and +e_x on the right: both give -k h
    @test symeq(ref[BoxFace(Symbolics.unwrap(x), :lo)], -P[1] * P[2])
    @test symeq(ref[BoxFace(Symbolics.unwrap(x), :hi)], -P[1] * P[2])
    @test symeq(ref[BoxFace(Symbolics.unwrap(y), :hi)], 2s.U[1] - 3)

    # a mesh discretization of the same box keys by piece, with the faces' numeric normals
    sys = system(traced)
    mesh = PieceMeshDisc(t)
    mv = VariableMap(sys, mesh)
    an = analyse(sys, mv, boundarymap(traced, mv, mesh))
    @test an isa StrongFormAnalysis{typeof(Γl)}
    face = Dict(Γl => BoxFace(Symbolics.unwrap(x), :lo), Γr => BoxFace(Symbolics.unwrap(x), :hi), Γt => BoxFace(Symbolics.unwrap(y), :hi))
    @test length(an.natural) == 3 && all(b -> symeq(b.Φ, ref[face[b.piece]]), an.natural)
    @test only(an.essential).piece == Γb && symeq(only(an.essential).g, only(reference.essential).g)

    ctx = SlotContext(system(dxform), v)
    @test symeq(slot_condition(dnform[2], BoxFace(Symbolics.unwrap(x), :lo), ctx, s), -s.G[1] - P[2])
    @test symeq(slot_condition(dnform[2], Γl, ctx, s), -s.G[1] - P[2])
    # a box-face piece restricts its coordinate like its face
    @test symeq(restrict_to_face(s.X[1] + s.X[2], Γr, ctx, s), 1 + s.X[2])
    # Dx(u) ~ h on the left face has the opposite sign of the normal derivative
    @test symeq(slot_condition(Dx(u(t, Γl)) ~ h, Γl, ctx, s), s.G[1] - P[2])
    @test occursin(
        "takes the normal derivative on Γ_right",
        errmsg(() -> slot_condition(Dn(Γr)(u(t, Γr)) ~ 0, Γl, ctx, s))
    )
end

@testset "mesh pieces end to end: symbolic normal" begin
    Ω = MeshDomain(PieceMesh())
    inlet, wall, outlet = (boundarypiece(Ω, n) for n in (:inlet, :wall, :outlet))
    domains = [(x, y) ∈ Ω, t ∈ Interval(0.0, 1.0)]
    function meshanalyse(eq, bcs)
        @named sys = PDESystem([eq], bcs, domains, [t, x, y], [u(t, x, y)], [k => 1.0, h => 2.0])
        disc = PieceMeshDisc(t)
        try
            SciMLBase.symbolic_discretize(sys, disc)
        catch e
            e isa StopAfterAnalysis || rethrow()
        end
        bmap, v = disc.captured[]
        return analyse(sys, v, bmap), v, bmap, sys
    end
    bcs = [
        u(0, x, y) ~ 0, u(t, inlet) ~ sin(t) * x, Dn(wall)(u(t, wall)) ~ h,
        k * Dn(outlet)(u(t, outlet)) + 2u(t, outlet) ~ 3,
    ]
    an, v, bmap, msys = meshanalyse(heat(k), bcs)
    @test only(bmap[operation(Symbolics.unwrap(u(t, x, y)))][wall]) isa PieceBoundary
    s = an.slots
    U, X, T, N, P = s.U[1], s.X, s.T, s.N, s.P
    @test an isa StrongFormAnalysis{typeof(inlet)}
    ess = only(an.essential)
    # essential data keeps its coordinates: there is no face to restrict them to
    @test ess.piece == inlet && ess.field == 1 && symeq(ess.g, sin(T) * X[1])
    nat = Dict(b.piece => b for b in an.natural)
    @test length(nat) == 2 && all(b -> b.eq == 1, an.natural)
    # Neumann against k∇u with the symbolic normal: -k h
    @test symeq(nat[wall].Φ, -P[1] * P[2])
    # Robin: facet bilinear plus facet source
    @test symeq(nat[outlet].Φ, 2U - 3)
    ctx = SlotContext(msys, v)
    @test symeq(slot_condition(bcs[3], wall, ctx, s), N[1] * s.G[1] + N[2] * s.G[2] - P[2])
    @test isequal(restrict_to_face(X[1], wall, ctx, s), X[1])

    # k(u): a nonlinear facet flux
    an, _, _ = meshanalyse(heat(1 + u(t, x, y)^2), [u(0, x, y) ~ 0, Dn(wall)(u(t, wall)) ~ h])
    @test symeq(only(an.natural).Φ, -(1 + an.slots.U[1]^2) * an.slots.P[2])

    # anisotropic flux: the normal derivative is not the conormal condition
    aniso = Dt(u(t, x, y)) ~ Dx(2Dx(u(t, x, y)) + Dy(u(t, x, y))) + Dy(Dy(u(t, x, y)))
    @test occursin(
        "write the conormal condition",
        errmsg(() -> meshanalyse(aniso, [u(0, x, y) ~ 0, Dn(wall)(u(t, wall)) ~ h]))
    )
    # a raw Cartesian derivative on a mesh piece
    msg = errmsg(() -> meshanalyse(heat(k), [u(0, x, y) ~ 0, Dx(u(t, wall)) ~ h]))
    @test occursin("write the normal derivative", msg)
    # piece mismatch, by the classifier in parse_bcs and by the slot step itself
    @test occursin(
        "not on Γ_wall", errmsg(() -> meshanalyse(heat(k), [u(0, x, y) ~ 0, Dn(wall)(u(t, inlet)) ~ h]))
    )
    @test occursin(
        "takes the normal derivative on Γ_wall",
        errmsg(() -> slot_condition(Dn(wall)(u(t, inlet)) ~ 0, inlet, ctx, s))
    )
    @test occursin(
        "still contains the trace u(t, Γ_inlet)",
        errmsg(() -> slot_condition(u(t, inlet) ~ 0, wall, ctx, s))
    )
    # periodic between mesh pieces needs an inferred map
    @test occursin(
        "not supported yet", errmsg(() -> meshanalyse(heat(k), [u(0, x, y) ~ 0, u(t, inlet) ~ u(t, outlet)]))
    )
    # over-specification on a piece
    @test occursin(
        "impose both an essential and a natural condition on field 1",
        errmsg(() -> meshanalyse(heat(k), [u(0, x, y) ~ 0, u(t, wall) ~ 0, Dn(wall)(u(t, wall)) ~ 1]))
    )
end

@testset "a map mixing box faces and pieces has a union key type" begin
    rect = Rectangle(0.0 .. 1.0, 0.0 .. 2.0)
    notch = BoundaryPiece(rect, :notch)
    box = [u(0, x, y) ~ 0, u(t, 0, y) ~ u(t, 1, y), u(t, x, 0) ~ 1]
    traced = [u(0, x, y) ~ 0, Dn(notch)(u(t, notch)) ~ h]
    @named sys = PDESystem(
        [heat(k)], [box; traced[2]], [(x, y) ∈ rect, t ∈ Interval(0.0, 1.0)], [t, x, y],
        [u(t, x, y)], [k => 1.0, h => 2.0]
    )
    v = VariableMap(sys, PieceBoxDisc(t))
    uop = operation(Symbolics.unwrap(u(t, x, y)))
    bmap = Dict(
        uop => Dict{Any, Any}(
            boundarymap(box, v, PieceBoxDisc(t))[uop]...,
            notch => boundarymap(traced, v, PieceMeshDisc(t))[uop][notch]
        )
    )
    s = StateSlots(1, 2, 2)
    ess, nat, per = boundary_conditions(bmap, SlotContext(sys, v), s, [[-s.P[1] * s.G[1], -s.P[1] * s.G[2]]])
    K = Union{BoxFace, typeof(notch)}
    @test eltype(ess) == EssentialBC{K} && eltype(per) == PeriodicBC{K}
    @test only(ess).piece == BoxFace(Symbolics.unwrap(y), :lo)
    @test only(nat).piece == notch && symeq(only(nat).Φ, -s.P[1] * s.P[2])
    p = only(per)
    @test p.master == BoxFace(Symbolics.unwrap(x), :lo) && p.mirror == BoxFace(Symbolics.unwrap(x), :hi)
    # conditions convert between key types, so outputs keyed by one key type join a mixed set
    @test convert(EssentialBC{K}, EssentialBC(notch, 1, Num(0))) isa EssentialBC{K}
    @test convert(NaturalBC{K}, NaturalBC(notch, 1, Num(0))).piece == notch
    @test convert(PeriodicBC{K}, PeriodicBC(p.master, p.mirror, 1, 1.0)) isa PeriodicBC{K, Float64}
    @test_throws MethodError convert(EssentialBC{BoxFace}, EssentialBC(notch, 1, Num(0)))
end
