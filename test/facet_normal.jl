using PDEBase, ModelingToolkit, Symbolics, SciMLBase, SymbolicUtils, Test
using DomainSets: Interval, Rectangle, (..)
using PDEBase: BoxFace, FacetNormal, PieceBoundary, boxface, normal
using PDEBase.StrongForm: analyse, StateSlots, SlotContext, slot_condition, check_whitelist, NaturalBC, match_flux, route_natural

cross(a, b) = [a[2] * b[3] - a[3] * b[2], a[3] * b[1] - a[1] * b[3], a[1] * b[2] - a[2] * b[1]]

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
# -(n × curl E)_i + γ (n × (n × E))_i, written in components
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
    @test occursin("calls on 2 locations", errmsg(() -> classify_condition(v, n[1] * E1(t, Γ) ~ E1(t, Γin))))
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

@testset "the component formula is the cross-product form" begin
    sys, v, _ = meshanalysis([absorbing(1), ics...])
    ctx = SlotContext(sys, v)
    s = StateSlots(3, 3, 1)
    curlE(Gm) = [Gm[3, 2] - Gm[2, 3], Gm[1, 3] - Gm[3, 1], Gm[2, 1] - Gm[1, 2]]  # Gm[f, c] = ∂_c E_f
    for (Nv, Gm, Ev, γv) in (
            ([1.0, 0, 0], [0 0 0; 1 0 0; 0 0 0], [0.0, 0, 0], 0.0),  # E = (0, x, 0), n = e_x
            ([0.6, 0, 0.8], [0.3 -1.2 0.7; 2.0 0.5 -0.4; -0.9 1.1 0.2], [0.5, -1.5, 2.0], 1.7),
        )
        vals = Dict{Any, Float64}(s.N[d] => Nv[d] for d in 1:3)
        for f in 1:3, c in 1:3
            vals[s.G[(f - 1) * 3 + c]] = Gm[f, c]
        end
        for f in 1:3
            vals[s.U[f]] = Ev[f]
        end
        vals[s.P[1]] = γv
        expected = -cross(Nv, curlE(Gm)) + γv * cross(Nv, cross(Nv, Ev))
        for i in 1:3
            B = slot_condition(absorbing(i), Γ, ctx, s)
            @test Symbolics.value(substitute(B, vals)) ≈ expected[i]
        end
    end
    # E = (0, x, 0), n = e_x: n × curl E = (0, -1, 0)
    @test cross([1.0, 0, 0], curlE([0 0 0; 1 0 0; 0 0 0])) == [0.0, -1.0, 0.0]
end

@testset "a scale that holds a normal component is not cancelled" begin
    s = StateSlots(1, 2, 0)
    N, G = s.N, s.G
    Φ = G[1] * N[1] + G[2] * N[2]
    @test match_flux(Φ, Φ, s) !== nothing
    @test match_flux(Φ, N[1] * Φ, s) === nothing
    # a scale that is provably positive for every real normal is cancelled
    @test symeq(match_flux(Φ, (1 + N[1]^2) * (Φ - 3), s), 3)
    @test symeq(last(route_natural((1 + N[1]^2) * (Φ - 3), [[G[1], G[2]]], N, s)), 3)
    @test symeq(match_flux(Φ, -(1 + N[1]^2) * (Φ - 3), s), 3)
    @test symeq(match_flux(Φ, (1 + (N[1] + N[2])^2) * (Φ - 3), s), 3)
    @test symeq(last(route_natural((1 + (N[1] + N[2])^2) * (Φ - 3), [[G[1], G[2]]], N, s)), 3)
    s1 = StateSlots(1, 1, 0)
    @test symeq(match_flux(s1.G[1] * s1.N[1], s1.N[1]^2 * (s1.G[1] * s1.N[1] - 3), s1), 3)
    @test symeq(match_flux(Φ, (N[1]^2 + 2 * N[2]^2) * (Φ - 3), s), 3)
    @test symeq(match_flux(Φ, (N[1]^2 + 2 * N[2]^2 - 1 // 2) * (Φ - 3), s), 3)
    @test match_flux(Φ, (N[1]^2 + 2 * N[2]^2 - 1) * (Φ - 3), s) === nothing
    @test symeq(match_flux(Φ, (3 // 2 - N[1]^2) * (Φ - 3), s), 3)
    @test symeq(last(route_natural((3 // 2 - N[1]^2) * (Φ - 3), [[G[1], G[2]]], N, s)), 3)
    @test match_flux(Φ, (1 - N[1]^2) * (Φ - 3), s) === nothing
    # exact bounds: no integer overflow, and no scale divided by something that can vanish
    Φ2 = N[2] * G[1]
    @test match_flux(Φ2, (-(2^62) - (2^62 + 1) * N[1]^2 + 2^62 * N[2]^2) * Φ2, s) === nothing
    @test match_flux(Φ, ((N[1]^4 - 1) / (N[1]^2 - 1)) * Φ, s) === nothing
    @test symeq(match_flux(Φ, Φ - 1 / (2 + N[1]^2), s), 1 / (2 + N[1]^2))
    @test symeq(match_flux(Φ, (Φ - 3) / (2 + N[1]^2), s), 3)
    # a denominator that is at least 1 for every real state and normal
    s1u = StateSlots(1, 2, 0)
    Φu = s1u.G[1] * s1u.N[1] + s1u.G[2] * s1u.N[2]
    dpos = 1 + s1u.N[1]^2 + s1u.U[1]^2
    @test symeq(match_flux(Φu, Φu - 1 / dpos, s1u), 1 / dpos)
    @test symeq(last(route_natural(Φu - 1 / dpos, [[s1u.G[1], s1u.G[2]]], s1u.N, s1u)), 1 / dpos)
    @test match_flux(Φu, Φu - 1 / (1 - s1u.N[1]^2 - s1u.U[1]^2), s1u) === nothing
    @test symeq(last(route_natural((Φ - 3) / (2 + N[1]^2), [[G[1], G[2]]], N, s)), 3)
    @test symeq(last(route_natural((N[1]^2 + 2 * N[2]^2) * (Φ - 3), [[G[1], G[2]]], N, s)), 3)
    # the unit normal has length one
    @test symeq(match_flux(Φ, (N[1]^2 + N[2]^2) * (Φ - 3), s), 3)
    @test match_flux(Φ, (N[1]^2 + N[1]^2) * (Φ - 3), s) === nothing
    @test match_flux(Φ, (N[1] + N[2])^2 * (Φ - 3), s) === nothing
    @test match_flux(Φ, (N[1]^2) * (Φ - 3), s) === nothing
    @test match_flux(Φ, (1 - N[1]^2) * (Φ - 3), s) === nothing
    @test match_flux(Φ, (1 + N[1]) * (Φ - 3), s) === nothing
    # scalar heat on a mesh: n₁ ∂ₙu ~ 0 vanishes on a facet with n = (0, 1), so it is not the natural condition
    heat3 = Dt(u(t, x, y, z)) ~ sum(Dxyz[c](Dxyz[c](u(t, x, y, z))) for c in 1:3)
    ic = u(0, x, y, z) ~ 0
    for bc in (n[1] * NormalDerivative(Γ)(u(t, Γ)) ~ 0, NormalDerivative(Γ)(u(t, Γ)) * n[1] ~ 0)
        sysh = PDESystem([heat3], [bc, ic], meshdomains, [t, x, y, z], [u(t, x, y, z)], [γ]; name = :sys)
        disc = NormalMeshDisc(t, Ref{Any}(nothing))
        @test_throws StopAtMap SciMLBase.symbolic_discretize(sysh, disc)
        bmaph, vh = disc.captured[]
        @test occursin("not a conormal condition", errmsg(() -> analyse(sysh, vh, bmaph)))
    end
    # without the normal factor the same condition is the natural condition with zero data
    sysh = PDESystem([heat3], [NormalDerivative(Γ)(u(t, Γ)) ~ 0, ic], meshdomains, [t, x, y, z], [u(t, x, y, z)], [γ]; name = :sys)
    disc = NormalMeshDisc(t, Ref{Any}(nothing))
    @test_throws StopAtMap SciMLBase.symbolic_discretize(sysh, disc)
    bmaph, vh = disc.captured[]
    anh = analyse(sysh, vh, bmaph)
    @test only(anh.natural).eq == 1 && symeq(only(anh.natural).Φ, 0)
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

@testset "a normal in a condition with calls on two locations is rejected" begin
    rect = Rectangle(0.0 .. 1.0, 0.0 .. 2.0)
    Γl, Γr, Γt = (boundarypiece(rect, nm) for nm in (:left, :right, :top))
    heat = Dt(u(t, x, y)) ~ Dx(k * Dx(u(t, x, y))) + Dy(k * Dy(u(t, x, y)))
    domains = [(x, y) ∈ rect, t ∈ Interval(0.0, 1.0)]
    bcs = [u(0, x, y) ~ 0, normal(Γt)[2] * u(t, Γl) ~ u(t, Γr), u(t, x, 0.0) ~ 0, u(t, x, 2.0) ~ 0]
    sys = PDESystem([heat], bcs, domains, [t, x, y], [u(t, x, y)], [k => 1.0]; name = :sys)
    v = VariableMap(sys, NormalBoxDisc(t))
    @test occursin("calls on 2 locations", errmsg(() -> classify_condition(v, bcs[2])))
    @test occursin("calls on 2 locations", errmsg(() -> to_box_bcs(sys; time = t)))
end

@testset "box lowering rejects the normal under a derivative and in an essential condition" begin
    rect = Rectangle(0.0 .. 1.0, 0.0 .. 2.0)
    Γl = boundarypiece(rect, :left)
    nl = normal(Γl)
    heat = Dt(u(t, x, y)) ~ Dx(k * Dx(u(t, x, y))) + Dy(k * Dy(u(t, x, y)))
    domains = [(x, y) ∈ rect, t ∈ Interval(0.0, 1.0)]
    function lowering(bc)
        sys = PDESystem([heat], [u(0, x, y) ~ 0, bc, u(t, x, 0.0) ~ 0, u(t, x, 2.0) ~ 0, u(t, 1.0, y) ~ 0], domains, [t, x, y], [u(t, x, y)], [k => 1.0]; name = :sys)
        return errmsg(() -> to_box_bcs(sys; time = t))
    end
    @test occursin("derivative of an expression containing the facet normal", lowering(Dx(nl[1]) + Dx(u(t, Γl)) ~ 1))
    @test occursin("contains the facet normal but no derivative of a trace", lowering(nl[1] * u(t, Γl) ~ 1))
    @test occursin("derivative of an expression containing the facet normal", lowering(NormalDerivative(Γl)(nl[1] * u(t, Γl)) ~ 0))
    @test occursin("contains the facet normal but no derivative of a trace", lowering(nl[1] * Dt(u(t, Γl)) ~ 1))
    # a derivative whose coefficient is a zero normal component, or that cancels across the sides, does not count
    zero_coefficient = nl[2] * Dx(u(t, Γl)) + nl[1] * u(t, Γl) ~ 1
    cancelled = nl[1] * Dx(u(t, Γl)) + u(t, Γl) ~ nl[1] * Dx(u(t, Γl))
    for bc in (zero_coefficient, cancelled)
        @test occursin("contains the facet normal but no derivative of a trace", lowering(bc))
    end
    sysbc = PDESystem([heat], [u(0, x, y) ~ 0, zero_coefficient, u(t, x, 0.0) ~ 0, u(t, x, 2.0) ~ 0, u(t, 1.0, y) ~ 0], domains, [t, x, y], [u(t, x, y)], [k => 1.0]; name = :sys)
    vbc = VariableMap(sysbc, NormalBoxDisc(t))
    @test occursin(
        "contains the facet normal but no derivative of a trace",
        errmsg(() -> analyse(sysbc, vbc, boundarymap(sysbc.bcs, vbc, NormalBoxDisc(t))))
    )
    @test occursin(
        "contains the facet normal but no derivative of a trace",
        errmsg(() -> slot_condition(nl[1] * u(t, Γl) ~ 1, Γl, SlotContext(sysbc, vbc), StateSlots(1, 2, 0)))
    )
    # the zero natural condition at a degenerate face is accepted, a nonzero gradient-free residual is not
    degenerate = x * nl[1] * Dx(u(t, Γl)) ~ 0
    sysdeg = PDESystem([Dt(u(t, x, y)) ~ Dx(x * Dx(u(t, x, y))) + Dy(Dy(u(t, x, y)))], [u(0, x, y) ~ 0, degenerate, u(t, x, 0.0) ~ 0, u(t, x, 2.0) ~ 0, u(t, 1.0, y) ~ 0], domains, [t, x, y], [u(t, x, y)], [k => 1.0]; name = :sys)
    @test isequal(Symbolics.simplify(to_box_bcs(sysdeg; time = t).bcs[2].lhs), Symbolics.simplify(-x * Dx(u(t, 0.0, y))))
    vdeg = VariableMap(sysdeg, NormalBoxDisc(t))
    andeg = analyse(sysdeg, vdeg, boundarymap(sysdeg.bcs, vdeg, NormalBoxDisc(t)))
    @test only(andeg.natural).eq == 1 && symeq(only(andeg.natural).Φ, 0)
    @test occursin("no derivative of a trace", lowering(x * nl[1] * u(t, Γl) ~ 1))
    # the classifier reads the differentiation variable's normal too
    Γr = boundarypiece(rect, :right)
    @test occursin("calls on 2 locations", errmsg(() -> classify_condition(vdeg, Differential(nl[1])(u(t, Γl)) ~ u(t, Γr))))
    @test occursin("uses the normal of", errmsg(() -> classify_condition(vdeg, Differential(normal(Γr)[1])(u(t, Γl)) ~ 0)))
    @test occursin("no derivative of a trace", lowering(nl[2] * u(t, Γl) ~ 0))
    # a derivative that vanished with the numeric normal imposes nothing
    @test occursin("no derivative of a trace", lowering(nl[2] * Dx(u(t, Γl)) ~ 0))
    @test occursin("no derivative of a trace", lowering(nl[1] * Dx(u(t, Γl)) ~ -Dx(u(t, Γl))))
    # a mixed derivative of a field that lacks one of the variables is zero
    @variables w2(..)
    sysmix = PDESystem([heat, Dt(w2(t, y)) ~ 0], [u(0, x, y) ~ 0, Dx(Dy(w2(t, y))) + nl[1] * u(t, Γl) ~ 1, u(t, x, 0.0) ~ 0, u(t, x, 2.0) ~ 0, u(t, 1.0, y) ~ 0], domains, [t, x, y], [u(t, x, y), w2(t, y)], [k => 1.0]; name = :sys)
    @test occursin("no derivative of a trace", errmsg(() -> to_box_bcs(sysmix; time = t)))
    @variables w3(..)
    sysx2 = PDESystem([heat, Dt(w3(t, x)) ~ 0], [u(0, x, y) ~ 0, Dx(Dy(w3(t, 0.0))) + nl[1] * u(t, Γl) ~ 1, u(t, x, 0.0) ~ 0, u(t, x, 2.0) ~ 0, u(t, 1.0, y) ~ 0], domains, [t, x, y], [u(t, x, y), w3(t, x)], [k => 1.0]; name = :sys)
    @test occursin("no derivative of a trace", errmsg(() -> to_box_bcs(sysx2; time = t)))
    # the two reviewer examples: a parameter derivative is zero, a coefficient does not hide a derivative
    @variables w4(..)
    sysk = PDESystem([heat, Dt(w4(t, y)) ~ 0], [u(0, x, y) ~ 0, Dy(Differential(k)(w4(t, y))) + nl[1] * u(t, Γl) ~ 1, u(t, x, 0.0) ~ 0, u(t, x, 2.0) ~ 0, u(t, 1.0, y) ~ 0], domains, [t, x, y], [u(t, x, y), w4(t, y)], [k => 1.0]; name = :sys)
    @test occursin("no derivative of a trace", errmsg(() -> to_box_bcs(sysk; time = t)))
    sysc = PDESystem([heat, Dt(w4(t, y)) ~ 0], [u(0, x, y) ~ 0, Dx(Dy(x * w4(t, y))) + nl[1] * u(t, Γl) ~ 1, u(t, x, 0.0) ~ 0, u(t, x, 2.0) ~ 0, u(t, 1.0, y) ~ 0], domains, [t, x, y], [u(t, x, y), w4(t, y)], [k => 1.0]; name = :sys)
    gotc = to_box_bcs(sysc; time = t).bcs[2]
    @test symeq(gotc.lhs - gotc.rhs, Dx(Dy(x * w4(t, y))) - u(t, 0.0, y) - 1)
    # mixed partials cancel across spellings, in both orders
    for (inner, outer) in ((Dy(Dx(u(t, 0.0, y))), Dx(Dy(u(t, Γl)))), (Dx(Dy(u(t, 0.0, y))), Dy(Dx(u(t, Γl)))))
        @test occursin("no derivative of a trace", lowering(nl[1] * outer + u(t, Γl) ~ -inner))
    end
    # a derivative of data is no derivative of the solution, and a placeholder name is no collision
    @parameters a(..) _dg1
    @test occursin("no derivative of a trace", lowering(Dx(a(x)) + nl[1] * u(t, Γl) ~ 1))
    @test occursin("no derivative of a trace", lowering(x * Dx(u(t, Γl)) + nl[1] * _dg1 * u(t, Γl) ~ 1))
    @test occursin("no derivative of a trace", lowering(x * Dx(u(t, Γl)) + nl[1] * Symbolics.variable(Symbol("##derivative#1")) * u(t, Γl) ~ 1))
    # a spatial derivative inside a time derivative counts, in either order
    @test to_box_bcs(PDESystem([heat], [u(0, x, y) ~ 0, nl[1] * Dt(Dx(u(t, Γl))) ~ 0, u(t, x, 0.0) ~ 0, u(t, x, 2.0) ~ 0, u(t, 1.0, y) ~ 0], domains, [t, x, y], [u(t, x, y)], [k => 1.0]; name = :sys); time = t) isa PDESystem
    # Dx of a time-only field is zero and does not count
    @variables w(..)
    sysw = PDESystem([heat, Dt(w(t)) ~ 0], [u(0, x, y) ~ 0, Dx(w(t)) + nl[1] * u(t, Γl) ~ 1, u(t, x, 0.0) ~ 0, u(t, x, 2.0) ~ 0, u(t, 1.0, y) ~ 0], domains, [t, x, y], [u(t, x, y), w(t)], [k => 1.0]; name = :sys)
    @test occursin("no derivative of a trace", errmsg(() -> to_box_bcs(sysw; time = t)))
    vw = VariableMap(sysw, NormalBoxDisc(t))
    sw = StateSlots(2, 2, 0)
    @test occursin("no derivative of a trace", errmsg(() -> slot_condition(nl[1] * Dx(w(t)) + u(t, Γl) ~ 0, Γl, SlotContext(sysw, vw), sw)))
    # the normal derivative of w(t, y) has no x part
    @variables w5(..)
    sysnd = PDESystem([heat, Dt(w5(t, y)) ~ 0], [u(0, x, y) ~ 0, u(t, x, 0.0) ~ 0], domains, [t, x, y], [u(t, x, y), w5(t, y)], [k => 1.0]; name = :sys)
    vnd = VariableMap(sysnd, NormalBoxDisc(t))
    @test symeq(slot_condition(NormalDerivative(Γl)(w5(t, Γl)) ~ 0, Γl, SlotContext(sysnd, vnd), StateSlots(2, 2, 0)), 0)
    # operators that cancel between the sides without a normal are unchanged
    @test symeq(slot_condition(Dt(u(t, Γl)) + Dx(u(t, Γl)) ~ Dt(u(t, Γl)), Γl, SlotContext(sysnd, vnd), StateSlots(2, 2, 0)), StateSlots(2, 2, 0).G[1])
    @test occursin("derivative of an expression containing the facet normal", errmsg(() -> slot_condition(Dx(nl[1]) + nl[1] * Dx(u(t, Γl)) ~ Dx(nl[1]), Γl, SlotContext(sysnd, vnd), StateSlots(2, 2, 0))))
    cancelled_nd = NormalDerivative(Γl)(nl[1] * u(t, Γl))
    @test occursin("derivative of an expression containing the facet normal", errmsg(() -> slot_condition(cancelled_nd + Dx(u(t, Γl)) ~ cancelled_nd, Γl, SlotContext(sysnd, vnd), StateSlots(2, 2, 0))))
    # division by a zero normal component is no condition
    nonfinite = -Dx(u(t, Γl)) ~ 1 / nl[2]
    @test_throws ArgumentError to_box_bcs(PDESystem([heat], [u(0, x, y) ~ 0, nonfinite, u(t, x, 0.0) ~ 0, u(t, x, 2.0) ~ 0, u(t, 1.0, y) ~ 0], domains, [t, x, y], [u(t, x, y)], [k => 1.0]; name = :sys); time = t)
    @test occursin("divides by a component of the facet normal", errmsg(() -> slot_condition(nonfinite, Γl, SlotContext(sysnd, vnd), StateSlots(2, 2, 0))))
    # the denominator vanishes only once the face coordinate and the zero normal component are both substituted
    endpoint_zero = nl[1] * Dx(u(t, Γl)) ~ 1 / (x + nl[2])
    @test_throws ArgumentError to_box_bcs(PDESystem([heat], [u(0, x, y) ~ 0, endpoint_zero, u(t, x, 0.0) ~ 0, u(t, x, 2.0) ~ 0, u(t, 1.0, y) ~ 0], domains, [t, x, y], [u(t, x, y)], [k => 1.0]; name = :sys); time = t)
    @test occursin("divides by a component of the facet normal", errmsg(() -> slot_condition(endpoint_zero, Γl, SlotContext(sysnd, vnd), StateSlots(2, 2, 0))))
    s0 = StateSlots(2, 2, 0)
    @test occursin("divides by a component of the facet normal", errmsg(() -> route_natural(s0.G[1] * s0.N[1] + s0.G[2] * s0.N[2] - 1 / s0.N[2], [[s0.G[1], s0.G[2]]], [-1, 0], s0)))
    # a normal of another piece cannot cancel between the sides
    @test occursin("is the normal of", errmsg(() -> slot_condition(NormalDerivative(Γl)(u(t, Γl)) + normal(boundarypiece(rect, :right))[1] ~ normal(boundarypiece(rect, :right))[1], Γl, SlotContext(sysnd, vnd), StateSlots(2, 2, 0))))
    # a derivative in a parameter is no spatial derivative, and the normal cannot be the variable of one
    @test occursin("contains the facet normal but no derivative of a trace", lowering(Differential(k)(u(t, Γl)) + nl[1] * u(t, Γl) ~ 0))
    @test occursin("derivative of an expression containing the facet normal", lowering(Differential(nl[1])(Dx(u(t, Γl))) ~ 0))
    @test occursin("derivative of an expression containing the facet normal", lowering(Differential(nl[1])(u(t, 0.0, y)) ~ 0))
    @test_throws ArgumentError assert_no_traces(Differential(nl[1])(u(t, 0.0, y)))
    @test occursin("contains the facet normal but no derivative of a trace", lowering(x * Dx(u(t, Γl)) + nl[1] * u(t, Γl) ~ 1))
    # the endpoint can remove the derivative: x = 0 on the left face
    sysx = PDESystem([heat], [u(0, x, y) ~ 0, x * Dx(u(t, Γl)) + nl[1] * u(t, Γl) ~ 1, u(t, x, 0.0) ~ 0, u(t, x, 2.0) ~ 0, u(t, 1.0, y) ~ 0], domains, [t, x, y], [u(t, x, y)], [k => 1.0]; name = :sys)
    vx = VariableMap(sysx, NormalBoxDisc(t))
    @test occursin(
        "contains the facet normal but no derivative of a trace",
        errmsg(() -> analyse(sysx, vx, boundarymap(sysx.bcs, vx, NormalBoxDisc(t))))
    )
    # a derivative spelled in endpoint form still counts as the derivative of the Robin condition
    robin = nl[1] * Dx(u(t, 0.0, y)) + u(t, Γl) ~ 0
    sysrobin = PDESystem([heat], [u(0, x, y) ~ 0, robin, u(t, x, 0.0) ~ 0, u(t, x, 2.0) ~ 0, u(t, 1.0, y) ~ 0], domains, [t, x, y], [u(t, x, y)], [k => 1.0]; name = :sys)
    got = to_box_bcs(sysrobin; time = t).bcs[2]
    @test symeq(got.lhs - got.rhs, -Dx(u(t, 0.0, y)) + u(t, 0.0, y))
    # analyse rejects the essential condition with the normal on a box face as well
    sysess = PDESystem([heat], [u(0, x, y) ~ 0, nl[1] * u(t, Γl) ~ 1, u(t, x, 0.0) ~ 0, u(t, x, 2.0) ~ 0, u(t, 1.0, y) ~ 0], domains, [t, x, y], [u(t, x, y)], [k => 1.0]; name = :sys)
    vess = VariableMap(sysess, NormalBoxDisc(t))
    bmapess = boundarymap(sysess.bcs, vess, NormalBoxDisc(t))
    @test occursin("contains the facet normal but no derivative of a trace", errmsg(() -> analyse(sysess, vess, bmapess)))
end
