using PDEBase, ModelingToolkit, Symbolics, SciMLBase, SymbolicUtils, Test
using DomainSets: Interval, Rectangle, (..)
using PDEBase: BoxFace, FacetNormal, PieceBoundary, boxface, normal
using PDEBase.StrongForm: analyse, StateSlots, SlotContext, slot_condition, check_whitelist, NaturalBC

# A mesh backend with three spatial dimensions that stops `symbolic_discretize` at the boundary map.
struct NormalMesh end
PDEBase.meshdim(::NormalMesh) = 3
PDEBase.facetset_names(::NormalMesh) = (:inlet, :wall)

struct NormalMeshDisc <: PDEBase.AbstractMeshDiscretization
    t::Any
    captured::Base.RefValue{Any}
end
PDEBase.get_time(d::NormalMeshDisc) = d.t
struct StopAtMap <: Exception end
function PDEBase.check_boundarymap(bmap, v::VariableMap, d::NormalMeshDisc)
    d.captured[] = (bmap, v)
    throw(StopAtMap())
end

struct NormalBoxDisc <: PDEBase.AbstractBoxDiscretization
    t::Any
end
PDEBase.get_time(d::NormalBoxDisc) = d.t

symeq(a, b) = isequal(Symbolics.simplify(Symbolics.expand(Num(a) - Num(b))), 0)
function errmsg(f)
    try
        f()
    catch e
        return sprint(showerror, e)
    end
    return nothing
end
boundarymap(bcs, v, disc) =
    PDEBase.parse_bcs(bcs, v, Dict(z => d_orders(z, bcs, v) for z in all_ivs(v)), disc)

@parameters t x y z γ k α
@variables E1(..) E2(..) E3(..) u(..)
Dt, Dx, Dy, Dz = Differential(t), Differential(x), Differential(y), Differential(z)
Dxyz = (Dx, Dy, Dz)

# A vector-valued field E = (E1, E2, E3) on a mesh with a curl-curl equation: component i has the
# conormal flux -Σ_c N_c (G[i, c] - G[c, i]).
Ω = MeshDomain(NormalMesh())
Γ, Γin = boundarypiece(Ω, :wall), boundarypiece(Ω, :inlet)
Evol = (E1(t, x, y, z), E2(t, x, y, z), E3(t, x, y, z))
Etr(P) = (E1(t, P), E2(t, P), E3(t, P))
curlcurl = [Dt(Evol[i]) ~ sum(Dxyz[c](Dxyz[c](Evol[i]) - Dxyz[i](Evol[c])) for c in 1:3) for i in 1:3]
ics = [E1(0, x, y, z) ~ 0, E2(0, x, y, z) ~ 0, E3(0, x, y, z) ~ 0]
meshdomains = [(x, y, z) ∈ Ω, t ∈ Interval(0.0, 1.0)]
n = normal(Γ)
# (n × curl E)_i + γ (n × (n × E))_i, written in components
absorbing(i, Eg = Etr(Γ), nn = n) =
    sum(nn[c] * (Dxyz[c](Eg[i]) - Dxyz[i](Eg[c])) for c in 1:3) +
    γ * (nn[i] * sum(nn[d] * Eg[d] for d in 1:3) - Eg[i]) ~ 0
meshsystem(bcs) = PDESystem(curlcurl, bcs, meshdomains, [t, x, y, z], collect(Evol), [γ]; name = :sys)

function meshanalysis(bcs)
    sys = meshsystem(bcs)
    disc = NormalMeshDisc(t, Ref{Any}(nothing))
    @test_throws StopAtMap SciMLBase.symbolic_discretize(sys, disc)
    bmap, v = disc.captured[]
    return sys, v, bmap
end

@testset "normal(Γ) as a symbolic object" begin
    @test n isa Vector{Num} && length(n) == 3
    @test repr(n[1]) == "n_wall(1)"
    @test operation(unwrap(n[2])) isa FacetNormal
    @test operation(unwrap(n[2])).piece == Γ
    @test isequal(n, normal(boundarypiece(Ω, :wall)))
    @test hash(n[1]) == hash(normal(Γ)[1])
    @test !isequal(n[1], n[2]) && !isequal(n[1], normal(Γin)[1])
    @test FacetNormal(Γ) == FacetNormal(boundarypiece(Ω, :wall))
    @test FacetNormal(Γ) != FacetNormal(Γin)
    @test repr(FacetNormal(Γ)) == "n_wall"
    @test isequal(Symbolics.scalarize(n), n)
    @test isequal(Symbolics.derivative(3 * n[1] + n[2], n[1]), 3)
    @test isequal(Symbolics.simplify(n[1] - n[1]), 0)
    @test length(Symbolics.get_variables(n[1] * u(t, Γ))) == 2
    @test sum(n .* n) isa Num && isequal(sum(n .* [1, 0, 0]), n[1])
    @test [n[2] * 1 - n[3] * 0, n[3] * 0 - n[1] * 1, n[1] * 0 - n[2] * 0] isa Vector{Num}
    @test sum(n .* [1, 2, 3]) isa Num
    @test occursin("component 4", errmsg(() -> FacetNormal(Γ)(4)))
    @test occursin("component 0", errmsg(() -> FacetNormal(Γ)(0)))
    @test !Base.isexported(PDEBase, :normal)
    rules = Dict(unwrap(n[1]) => 7)
    @test isequal(substitute_traces(n[1] * Dx(u(t, Γ)), rules), 7 * Dx(u(t, Γ)))
    @test occursin("remains", errmsg(() -> assert_no_traces(2 * n[1])))
    @test d_orders(x, [n[1] * Dx(E1(t, Γ)) + E1(t, Γ) ~ 0]) == [1]
    @test d_orders(x, [n[1] * E1(t, Γ) ~ 0]) == Int[]
    @test PDEBase.differential_order(unwrap(n[1] * Dx(E1(t, Γ))), unwrap(x)) == Set([1])
    @test PDEBase.count_differentials(unwrap(n[1] * E1(t, Γ)), unwrap(x)) == 0
    @test PDEBase.bckind(absorbing(1)) == PDEBase.Robin()
    @test PDEBase.bckind(n[1] * Dx(E1(t, Γ)) ~ 0) == PDEBase.Natural()
end

@testset "classification of conditions with a normal" begin
    sys, v, _ = meshanalysis([absorbing(1), ics...])
    @test classify_condition(v, absorbing(1)).locations == [PDEBase.Piece(Γ)]
    msg = errmsg(() -> classify_condition(v, n[1] * E1(t, Γin) ~ 0))
    @test occursin("uses the normal of $(Γ) but is on $(Γin)", msg)
    @test occursin("has no call of a dependent variable", errmsg(() -> classify_condition(v, n[1] ~ 0)))
    bmap = boundarymap([absorbing(1), ics...], v, NormalMeshDisc(t, Ref{Any}(nothing)))
    b = only(reduce(vcat, [get(bykey, Γ, []) for bykey in values(bmap)]))
    @test b isa PieceBoundary && b.kind == PDEBase.Robin() && b.order == 1
end

@testset "slots on a mesh piece: symbolic normal" begin
    sys, v, _ = meshanalysis([absorbing(1), ics...])
    ctx = SlotContext(sys, v)
    s = StateSlots(3, 3, 1)
    N, G, U = s.N, s.G, s.U
    B = slot_condition(absorbing(2), Γ, ctx, s)
    expected = sum(N[c] * (G[(2 - 1) * 3 + c] - G[(c - 1) * 3 + 2]) for c in 1:3) +
        s.P[1] * (N[2] * sum(N[d] * U[d] for d in 1:3) - U[2])
    @test symeq(B, expected)
    # a Cartesian derivative of a trace becomes the gradient slot when the normal is present
    @test symeq(slot_condition(n[3] * Dy(E2(t, Γ)) ~ 0, Γ, ctx, s), N[3] * G[(2 - 1) * 3 + 2])
    # normal-free pieces are unchanged
    @test symeq(slot_condition(NormalDerivative(Γ)(E1(t, Γ)) ~ 0, Γ, ctx, s), sum(N[d] * G[d] for d in 1:3))
end

@testset "negative cases" begin
    sys, v, _ = meshanalysis([absorbing(1), ics...])
    ctx = SlotContext(sys, v)
    s = StateSlots(3, 3, 1)
    cartesian = errmsg(() -> slot_condition(Dx(E1(t, Γ)) ~ 0, Γ, ctx, s))
    @test occursin("Cartesian derivative of a trace on $(Γ)", cartesian)
    @test occursin("include PDEBase.normal($(Γ))", cartesian)
    @test occursin(
        "is the normal of $(Γin), but the condition",
        errmsg(() -> slot_condition(n[1] * Dx(E1(t, Γ)) + normal(Γin)[1] ~ 0, Γ, ctx, s))
    )
    @test occursin(
        "derivative of an expression containing the facet normal",
        errmsg(() -> slot_condition(Dx(n[1]) * E1(t, Γ) ~ 0, Γ, ctx, s))
    )
    @test occursin(
        "derivative of an expression containing the facet normal",
        errmsg(() -> slot_condition(Dx(n[1] * E1(t, Γ)) ~ 0, Γ, ctx, s))
    )
    @test occursin("order 2", errmsg(() -> slot_condition(n[1] * Dx(Dx(E1(t, Γ))) ~ 0, Γ, ctx, s)))
    @test occursin("time derivative of a trace", errmsg(() -> slot_condition(n[1] * Dt(E1(t, Γ)) ~ 0, Γ, ctx, s)))
    # the normal is not part of the volume equations
    @test occursin("not supported by the strong-form analysis", errmsg(() -> check_whitelist(n[1] * Evol[1], [t, x, y, z], collect(Evol))))
    # a condition with the normal and no derivative of a trace
    for bad in (n[1] * E1(t, Γ) ~ 0, sum(n[c] * Etr(Γ)[c] for c in 1:3) ~ 0)
        sysbad, vbad, bmapbad = meshanalysis([bad, ics...])
        msg = errmsg(() -> analyse(sysbad, vbad, bmapbad))
        @test occursin("contains the facet normal but no derivative of a trace", msg)
    end
    # a normal left in a condition is not a lowered one
    @test_throws ArgumentError assert_no_traces(absorbing(1))
end

@testset "end to end: a three-component field on a mesh piece" begin
    sys, v, bmap = meshanalysis([absorbing(1), absorbing(2), absorbing(3), ics...])
    an = analyse(sys, v, bmap)
    s = an.slots
    @test an.dim == 3 && an.nfields == 3 && isempty(an.essential)
    @test sort([b.eq for b in an.natural]) == [1, 2, 3]
    for b in an.natural
        @test b.piece == Γ
        # the curl-type part is the conormal of its own equation, so only the impedance term is left
        @test symeq(b.Φ, s.P[1] * (s.N[b.eq] * sum(s.N[d] * s.U[d] for d in 1:3) - s.U[b.eq]))
    end
    # the curl-type term alone is a natural condition with zero data
    curl_only = [sum(n[c] * (Dxyz[c](Etr(Γ)[i]) - Dxyz[i](Etr(Γ)[c])) for c in 1:3) ~ 0 for i in 1:3]
    sys2, v2, bmap2 = meshanalysis([curl_only..., ics...])
    an2 = analyse(sys2, v2, bmap2)
    @test sort([b.eq for b in an2.natural]) == [1, 2, 3] && all(b -> symeq(b.Φ, 0), an2.natural)
    # an inhomogeneous right-hand side is data of the natural condition
    withdata = absorbing(1)
    sys3, v3, bmap3 = meshanalysis([withdata.lhs ~ γ, absorbing(2), absorbing(3), ics...])
    an3 = analyse(sys3, v3, bmap3)
    nat1 = only(filter(b -> b.eq == 1, an3.natural))
    @test symeq(nat1.Φ, s.P[1] * (s.N[1] * sum(s.N[d] * s.U[d] for d in 1:3) - s.U[1]) - s.P[1])
end

@testset "boxes: the normal is numeric and the lowering matches the endpoint form" begin
    rect = Rectangle(0.0 .. 1.0, 0.0 .. 2.0)
    Γl, Γr, Γb, Γt = (boundarypiece(rect, nm) for nm in (:left, :right, :bottom, :top))
    heat = Dt(u(t, x, y)) ~ Dx(k * Dx(u(t, x, y))) + Dy(k * Dy(u(t, x, y)))
    domains = [(x, y) ∈ rect, t ∈ Interval(0.0, 1.0)]
    impedance(Γ) = let nn = normal(Γ)
        k * (nn[1] * Dx(u(t, Γ)) + nn[2] * Dy(u(t, Γ))) + α * u(t, Γ) ~ 1
    end
    bcs = [u(0, x, y) ~ 0, impedance(Γl), impedance(Γr), impedance(Γb), u(t, Γt) ~ 0]
    sys = PDESystem([heat], bcs, domains, [t, x, y], [u(t, x, y)], [k => 1.0, α => 2.0]; name = :sys)
    lowered = to_box_bcs(sys; time = t).bcs
    expected = [
        u(0, x, y) ~ 0,
        -k * Dx(u(t, 0.0, y)) + α * u(t, 0.0, y) ~ 1,
        k * Dx(u(t, 1.0, y)) + α * u(t, 1.0, y) ~ 1,
        -k * Dy(u(t, x, 0.0)) + α * u(t, x, 0.0) ~ 1,
        u(t, x, 2.0) ~ 0,
    ]
    for (got, want) in zip(lowered, expected)
        @test isequal(Symbolics.simplify(got.lhs - got.rhs), Symbolics.simplify(want.lhs - want.rhs))
    end
    # analysis with the numeric normal agrees with the Dx form
    function boxanalysis(bcs)
        sys = PDESystem([heat], bcs, domains, [t, x, y], [u(t, x, y)], [k => 1.0, α => 2.0]; name = :sys)
        v = VariableMap(sys, NormalBoxDisc(t))
        return analyse(sys, v, boundarymap(bcs, v, NormalBoxDisc(t)))
    end
    an, ref = boxanalysis(bcs), boxanalysis([expected[1:4]; u(t, Γt) ~ 0])
    refmap = Dict(b.piece => b.Φ for b in ref.natural)
    @test length(an.natural) == 3 && all(b -> symeq(b.Φ, refmap[b.piece]), an.natural)
    @test !any(b -> any(sym -> any(m -> isequal(sym, unwrap(m)), an.slots.N), Symbolics.get_variables(b.Φ)), an.natural)
    # a piece that is not a face of the box has no numeric normal to lower to
    inlet = BoundaryPiece(rect, :inlet)
    badbcs = [u(0, x, y) ~ 0, normal(inlet)[1] * Dx(u(t, inlet)) ~ 0]
    badsys = PDESystem([heat], badbcs, domains, [t, x, y], [u(t, x, y)], [k => 1.0, α => 2.0]; name = :sys)
    @test occursin("is not a face of the box", errmsg(() -> to_box_bcs(badsys; time = t)))
end

@testset "the component index follows indvars, not the declared box axes" begin
    # u(t, y, x) makes indvars (y, x), while the box axes are (x, y) as declared by (x, y) ∈ rect
    rect = Rectangle(0.0 .. 1.0, 0.0 .. 2.0)
    Γl, Γt = boundarypiece(rect, :left), boundarypiece(rect, :top)
    heat = Dt(u(t, y, x)) ~ Dx(k * Dx(u(t, y, x))) + Dy(k * Dy(u(t, y, x)))
    domains = [(x, y) ∈ rect, t ∈ Interval(0.0, 1.0)]
    impedance(Γ) = let nn = normal(Γ)
        # position 1 of indvars is y, position 2 is x
        k * (nn[1] * Dy(u(t, Γ)) + nn[2] * Dx(u(t, Γ))) + α * u(t, Γ) ~ 1
    end
    bcs = [u(0, y, x) ~ 0, impedance(Γl), impedance(Γt), u(t, 0.0, x) ~ 0, u(t, y, 1.0) ~ 0]
    sys = PDESystem([heat], bcs, domains, [t, x, y], [u(t, y, x)], [k => 1.0, α => 2.0]; name = :sys)
    v = VariableMap(sys, NormalBoxDisc(t))
    @test isequal(indvars(v), [unwrap(y), unwrap(x)])
    @test boxface(v, Γl) == BoxFace(unwrap(x), :lo)
    @test boxface(v, Γt) == BoxFace(unwrap(y), :hi)
    lowered = to_box_bcs(sys; time = t).bcs
    # the left face is x = 0 with outward normal -e_x, which is component 2: -Dx; the top face is
    # y = 2 with +e_y, component 1: +Dy
    left = -k * Dx(u(t, y, 0.0)) + α * u(t, y, 0.0) ~ 1
    top = k * Dy(u(t, 2.0, x)) + α * u(t, 2.0, x) ~ 1
    for (got, want) in zip(lowered[2:3], (left, top))
        @test isequal(Symbolics.simplify(got.lhs - got.rhs), Symbolics.simplify(want.lhs - want.rhs))
    end
    # and the analysis sees the same numeric normal
    function swapped(bcs)
        sys = PDESystem([heat], bcs, domains, [t, x, y], [u(t, y, x)], [k => 1.0, α => 2.0]; name = :sys)
        v = VariableMap(sys, NormalBoxDisc(t))
        return analyse(sys, v, boundarymap(bcs, v, NormalBoxDisc(t)))
    end
    an = swapped(bcs)
    ref = swapped([bcs[1]; left; top; bcs[4:5]])
    refmap = Dict(b.piece => b.Φ for b in ref.natural)
    @test length(an.natural) == 2 && all(b -> symeq(b.Φ, refmap[b.piece]), an.natural)
    # the opposite assignment of components would put -Dy on the left face, which is wrong
    wrong = k * (normal(Γl)[2] * Dy(u(t, Γl)) + normal(Γl)[1] * Dx(u(t, Γl))) + α * u(t, Γl) ~ 1
    wronglowered = to_box_bcs(PDESystem([heat], [bcs[1]; wrong; bcs[3:5]], domains, [t, x, y], [u(t, y, x)], [k => 1.0, α => 2.0]; name = :sys); time = t).bcs[2]
    @test !isequal(Symbolics.simplify(wronglowered.lhs - wronglowered.rhs), Symbolics.simplify(left.lhs - left.rhs))
end
