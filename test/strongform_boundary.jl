using PDEBase, ModelingToolkit, Symbolics, Test
using DomainSets: Interval
using PDEBase.StrongForm: analyse, Weighted, StateSlots, SlotContext, BoxFace, EssentialBC, NaturalBC, PeriodicBC,
    match_flux, slot_condition, route_natural, essential_condition, periodic_condition,
    boundary_conditions, pair_equations

struct TimeDisc <: PDEBase.AbstractEquationSystemDiscretization
    t::Any
end
PDEBase.get_time(d::TimeDisc) = d.t

bcdata(y) = y^3
Symbolics.@register_symbolic bcdata(y)
customcoef(y) = abs(y - 20)
Symbolics.@register_symbolic customcoef(y)

symeq(a, b) = isequal(Symbolics.simplify(Num(a) - Num(b)), 0)

# Variable map, slot context and parse_bcs boundary map of a system.
function boxsetup(eqs, bcs, domains, ivs, dvs, ps, t)
    @named sys = PDESystem(eqs, bcs, domains, ivs, dvs, ps)
    v = VariableMap(sys, TimeDisc(t))
    orders = Dict(x => d_orders(x, bcs) for x in all_ivs(v))
    return v, PDEBase.parse_bcs(bcs, v, orders)
end

@testset "match_flux: natural-condition matching table (test_misc.out)" begin
    s = StateSlots(1, 2, 1)
    G, N, U, k = s.G, s.N, s.U[1], s.P[1]
    @variables h
    dn = N[1] * G[1] + N[2] * G[2]
    @test symeq(match_flux(-k * dn, dn - h, s), -k * h)
    @test symeq(match_flux(-k * dn, k * dn + 2U - h, s), 2U - h)
    @test symeq(match_flux(-(1 + U^2) * dn, dn - h, s), -(1 + U^2) * h)
    K = [2.0 0.5; 0.5 1.0]
    conorm = N[1] * (K[1, 1] * G[1] + K[1, 2] * G[2]) + N[2] * (K[2, 1] * G[1] + K[2, 2] * G[2])
    @test symeq(match_flux(-conorm, conorm - h, s), -h)
    @test match_flux(-conorm, dn - h, s) === nothing
    # not affine in the gradient
    @test match_flux(-dn, dn^2 - h, s) === nothing
    # only a symbolic proof is accepted (O3); the cross-multiplied, expanded check proves these
    @test symeq(match_flux(-(1 + U)^2 * dn, (1 + 2U + U^2) * dn - h, s), -h)
    @test symeq(match_flux(-sqrt(U) * (1 + U)^2 * dn, sqrt(U) * (1 + 2U + U^2) * dn - h, s), -h)
    # Codex counterexamples to sampling: non-smooth, registered, composed and underflowing
    @test match_flux(G[1] + U * G[2], G[1] + abs(U) * G[2] - 1, s) === nothing
    @test match_flux(G[1] + (2 - U) * G[2], G[1] + abs(U - 2) * G[2] - 1, s) === nothing
    @test match_flux(G[1] + (2 - U) * G[2], G[1] + ifelse(U > 20, U, 2 - U) * G[2] - 1, s) ===
        nothing
    @test match_flux(G[1] + (20 - U) * G[2], G[1] + customcoef(U) * G[2] - 1, s) === nothing
    @test match_flux(G[1] + (20 - U) * G[2], G[1] + sqrt((U - 20)^2) * G[2] - 1, s) === nothing
    @test match_flux(G[1] + exp(-1.0e6 * U^2) * G[2], G[1] - 1, s) === nothing
    # tiny coefficients must not underflow into a false proof
    @test match_flux(1.0e-200 * G[1], 1.0e-200 * (G[1] + G[2]) - 1, s) === nothing
    # numeric agreement only sharpens the rejection message
    @test_throws "agree numerically with the conormal of equation 1" route_natural(
        G[1] + tan(U) * G[2] - h, [[G[1] + sin(U) / cos(U) * G[2], 0]], [1, 0], s
    )
end

@testset "routing over several equations' fluxes" begin
    s = StateSlots(2, 2, 0)
    G = s.G
    @variables h
    lap(f) = [-G[2f - 1], -G[2f]]
    n = [1, 0]
    # Neumann on u routes to the u equation only
    e, Φ = route_natural(G[1] - h, [lap(1), lap(2)], n, s)
    @test e == 1 && symeq(Φ, -h)
    # m3: two equations with proportional conormals is an error naming both
    err = try
        route_natural(G[1] - h, [lap(1), 2 .* lap(1)], n, s; label = "Dx(u) ~ h")
    catch ex
        sprint(showerror, ex)
    end
    @test occursin("equations 1 and 2", err)
    # Dn of u against the coupled flux ∇u + χ∇w: write the conormal condition
    χ = 3
    @test_throws "not a conormal condition for any equation's flux; write the conormal condition" route_natural(
        G[1] - h, [lap(1) .+ χ .* lap(2)], n, s
    )
    # raw Dy on an x-face: write the normal derivative
    @test_throws "not a conormal condition for any equation's flux; write the normal derivative" route_natural(
        G[2] - h, [lap(1), lap(2)], n, s
    )
end

@testset "essential extraction" begin
    s = StateSlots(2, 1, 0)
    U, X, T = s.U, s.X, s.T
    f, g = essential_condition(2U[1] + 1 - sin(T), s)
    @test f == 1 && symeq(g, (sin(T) - 1) / 2)
    f, g = essential_condition(U[2] - X[1] * T, s; field = 2)
    @test f == 2 && symeq(g, X[1] * T)
    # m1: two-field essential condition, and u + w ~ 1, rejected cleanly
    @test_throws "Dirichlet data depends on another field" essential_condition(U[1] - U[2], s)
    @test_throws "Dirichlet data depends on another field" essential_condition(U[1] - U[2], s; field = 1)
    @test_throws "Dirichlet data depends on another field" essential_condition(U[1] + U[2] - 1, s; field = 1)
    @test_throws "Dirichlet data depends on another field" essential_condition(U[1] * U[2] - 1, s; field = 1)
    @test_throws "Dirichlet data must be linear in one field" essential_condition(U[1]^2 - 1, s)
    @test_throws "Dirichlet data must be linear in one field" essential_condition(T - 1, s)
end

@testset "box faces from parse_bcs, 2D heat with k" begin
    @parameters t x y k h
    @variables u(..)
    Dt, Dx, Dy = Differential(t), Differential(x), Differential(y)
    bcs = [
        u(0, x, y) ~ 0,
        u(t, 0, y) ~ t * sin(y),
        Dx(u(t, 1, y)) ~ h,
        Dy(u(t, x, 0)) ~ 0,
        k * Dy(u(t, x, 1)) + 2u(t, x, 1) ~ 3,
    ]
    domains = [t ∈ Interval(0.0, 1.0), x ∈ Interval(0.0, 1.0), y ∈ Interval(0.0, 1.0)]
    v, bmap = boxsetup(
        [Dt(u(t, x, y)) ~ Dx(k * Dx(u(t, x, y))) + Dy(k * Dy(u(t, x, y)))], bcs, domains,
        [t, x, y], [u(t, x, y)], [k => 1.0, h => 2.0], t
    )
    s = StateSlots(1, 2, 2)
    ctx = SlotContext(v)
    @test isequal(ctx.sivs, Symbolics.unwrap.([x, y]))
    G, U, X, T, P = s.G, s.U[1], s.X, s.T, s.P
    flux = [[-P[1] * G[1], -P[1] * G[2]]]
    ess, nat, per = boundary_conditions(bmap, ctx, s, flux)
    @test isempty(per)
    @test length(ess) == 1
    @test ess[1].piece == BoxFace(Symbolics.unwrap(x), :lo) && ess[1].field == 1
    @test symeq(ess[1].g, T * sin(X[2]))
    byface = Dict(b.piece => b for b in nat)
    @test length(byface) == 3 && all(b -> b.eq == 1, nat)
    # Dx(u) ~ h on the x-high face: -k h
    @test symeq(byface[BoxFace(Symbolics.unwrap(x), :hi)].Φ, -P[1] * P[2])
    # Dy(u(t, x, 0)) ~ 0 is a natural (zero-flux) condition, not an essential one (S1)
    @test symeq(byface[BoxFace(Symbolics.unwrap(y), :lo)].Φ, 0)
    # Robin: facet bilinear plus facet source
    @test symeq(byface[BoxFace(Symbolics.unwrap(y), :hi)].Φ, 2U - 3)
    # derivatives of trace-free boundary data are expanded and their coordinates slotted
    Bdata = slot_condition(u(t, 0, y) ~ Dy(sin(y)), BoxFace(Symbolics.unwrap(x), :lo), ctx, s)
    @test symeq(Bdata, U - cos(X[2]))
    # a registered function without a derivative rule errors instead of recursing forever
    @test_throws "register the derivative" slot_condition(
        u(t, 0, y) ~ Dy(bcdata(y)), BoxFace(Symbolics.unwrap(x), :lo), ctx, s
    )
end

@testset "coefficients are restricted to the face" begin
    @parameters t x y
    @variables u(..)
    Dt, Dx, Dy = Differential(t), Differential(x), Differential(y)
    bcs = [u(0, x, y) ~ 0, Dx(u(t, 0, y)) ~ 1, u(t, 1, y) ~ x + y]
    domains = [t ∈ Interval(0.0, 1.0), x ∈ Interval(0.0, 1.0), y ∈ Interval(0.0, 1.0)]
    v, bmap = boxsetup(
        [Dt(u(t, x, y)) ~ Dx(Dx(u(t, x, y)) + x * Dy(u(t, x, y))) + Dy(Dy(u(t, x, y)))],
        bcs, domains, [t, x, y], [u(t, x, y)], [], t
    )
    s = StateSlots(1, 2, 0)
    G, X = s.G, s.X
    # conormal on the low x-face is G[1] + X[1] G[2], which is G[1] at X[1] = 0
    ess, nat, per = boundary_conditions(bmap, SlotContext(v), s, [[-G[1] - X[1] * G[2], -G[2]]])
    @test length(nat) == 1 && symeq(only(nat).Φ, 1)
    @test length(ess) == 1 && symeq(only(ess).g, 1 + X[2])
end

@testset "1D: trace never derivative-expanded, raw Dx equals the Dn form" begin
    @parameters t x h
    @variables u(..)
    Dt, Dx = Differential(t), Differential(x)
    domains = [t ∈ Interval(0.0, 1.0), x ∈ Interval(0.0, 1.0)]
    bcs = [u(0, x) ~ 0, Dx(u(t, 0)) ~ 0, Dx(u(t, 1)) ~ h]
    v, bmap = boxsetup(
        [Dt(u(t, x)) ~ Dx(Dx(u(t, x)))], bcs, domains, [t, x], [u(t, x)], [h => 2.0], t
    )
    s = StateSlots(1, 1, 1)
    ctx = SlotContext(v)
    G, N, P = s.G, s.N, s.P
    flux = [[-G[1]]]
    lo, hi = BoxFace(Symbolics.unwrap(x), :lo), BoxFace(Symbolics.unwrap(x), :hi)
    # S1: the trace becomes a gradient slot instead of expanding to zero
    @test symeq(slot_condition(bcs[2], lo, ctx, s), G[1])
    ess, nat, per = boundary_conditions(bmap, ctx, s, flux)
    @test isempty(ess) && isempty(per)
    byface = Dict(b.piece => b for b in nat)
    @test symeq(byface[lo].Φ, 0)
    Φraw = byface[hi].Φ
    @test !symeq(Φraw, 0)
    # the Dn form N·G - h with N = +1 gives the same matched flux
    Bdn = substitute(N[1] * G[1] - P[1], Dict(N[1] => 1))
    @test symeq(Φraw, match_flux(-G[1], Bdn, s))
    @test symeq(Φraw, -P[1])
    # inhomogeneous Neumann on the low face: sign flips with the normal (rerun_misc: 2V)
    e, Φ = route_natural(slot_condition(Dx(u(t, 0)) ~ 2, lo, ctx, s), flux, [-1], s)
    @test e == 1 && symeq(Φ, 2)
    # a trace on another face is reported by name, before anything is expanded
    @test_throws "still contains the trace u(t, 1)" slot_condition(Dx(u(t, 0)) ~ u(t, 1), lo, ctx, s)
    @test_throws "time derivative of a trace" slot_condition(Dt(u(t, 0)) ~ 0, lo, ctx, s)
    @test_throws "order 2" slot_condition(Dx(Dx(u(t, 0))) ~ 0, lo, ctx, s)
    @test_throws "nested derivative" slot_condition(Dx(Dt(u(t, 0))) ~ 0, lo, ctx, s)
    @variables z
    @test_throws "derivative in z, which is not an independent variable" slot_condition(
        Differential(z)(u(t, 0)) ~ 0, lo, ctx, s
    )
    @test_throws "derivative of an expression containing a trace" slot_condition(
        Dx(u(t, 0)^2) ~ 0, lo, ctx, s
    )
end

@testset "two fields: essential coupling rejected (m1)" begin
    @parameters t x
    @variables u(..) w(..)
    Dt, Dx = Differential(t), Differential(x)
    domains = [t ∈ Interval(0.0, 1.0), x ∈ Interval(0.0, 1.0)]
    eqs = [Dt(u(t, x)) ~ Dx(Dx(u(t, x))), Dt(w(t, x)) ~ Dx(Dx(w(t, x)))]
    s = StateSlots(2, 1, 0)
    flux = [[-s.G[1]], [-s.G[2]]]
    for bc in (u(t, 0) ~ w(t, 0), u(t, 0) + w(t, 0) ~ 1)
        bcs = [u(0, x) ~ 0, w(0, x) ~ 0, bc]
        v, bmap = boxsetup(eqs, bcs, domains, [t, x], [u(t, x), w(t, x)], [], t)
        @test_throws "Dirichlet data depends on another field" boundary_conditions(
            bmap, SlotContext(v), s, flux
        )
    end
end

@testset "periodic" begin
    @parameters t x
    @variables u(..) w(..)
    Dt, Dx = Differential(t), Differential(x)
    domains = [t ∈ Interval(0.0, 1.0), x ∈ Interval(0.0, 2.0)]
    bcs = [u(0, x) ~ 0, u(t, 0) ~ u(t, 2)]
    v, bmap = boxsetup(
        [Dt(u(t, x)) ~ Dx(Dx(u(t, x)))], bcs, domains, [t, x], [u(t, x)], [], t
    )
    s = StateSlots(1, 1, 0)
    ess, nat, per = boundary_conditions(bmap, SlotContext(v), s, [[-s.G[1]]])
    @test isempty(ess) && isempty(nat)
    @test length(per) == 1
    p = only(per)
    @test p.master == BoxFace(Symbolics.unwrap(x), :lo) &&
        p.mirror == BoxFace(Symbolics.unwrap(x), :hi)
    @test p.field == 1 && p.shift == 2.0

    bcs = [u(0, x) ~ 0, w(0, x) ~ 0, u(t, 0) ~ w(t, 2)]
    v, bmap = boxsetup(
        [Dt(u(t, x)) ~ Dx(Dx(u(t, x))), Dt(w(t, x)) ~ Dx(Dx(w(t, x)))], bcs, domains,
        [t, x], [u(t, x), w(t, x)], [], t
    )
    s = StateSlots(2, 1, 0)
    @test_throws "must equate the same field" boundary_conditions(
        bmap, SlotContext(v), s, [[-s.G[1]], [-s.G[2]]]
    )
end

@testset "array parameters and undeclared symbols in boundary conditions" begin
    @parameters t x
    @parameters a[1:2]
    @variables u(..) c
    Dt, Dx = Differential(t), Differential(x)
    domains = [t ∈ Interval(0.0, 1.0), x ∈ Interval(0.0, 1.0)]
    bcs = [u(0, x) ~ 0, u(t, 0) ~ a[1], Dx(u(t, 1)) ~ a[2]]
    v, bmap = boxsetup([Dt(u(t, x)) ~ Dx(Dx(u(t, x)))], bcs, domains, [t, x], [u(t, x)], [a], t)
    ctx = SlotContext(v)
    s = StateSlots(1, 1, length(ctx.params))
    ess, nat, _ = boundary_conditions(bmap, ctx, s, [[-s.G[1]]])
    @test symeq(only(ess).g, s.P[1])
    @test symeq(only(nat).Φ, -s.P[2])
    @test_throws "unresolved symbols" slot_condition(
        u(t, 0) ~ c, BoxFace(Symbolics.unwrap(x), :lo), ctx, s
    )
end

@testset "conormal vanishing identically on a face (weighted r²)" begin
    @parameters t r
    @variables u(..)
    Dt, Dr = Differential(t), Differential(r)
    domains = [t ∈ Interval(0.0, 1.0), r ∈ Interval(0.0, 1.0)]
    eq = Dt(u(t, r)) ~ 1 / r^2 * Dr(r^2 * Dr(u(t, r)))
    function weighted(bcs)
        @named sys = PDESystem([eq], bcs, domains, [t, r], [u(t, r)])
        v = VariableMap(sys, TimeDisc(t))
        bmap = PDEBase.parse_bcs(bcs, v, Dict(x => d_orders(x, bcs) for x in all_ivs(v)))
        return analyse(sys, v, bmap; ibp = Dict(1 => Weighted(r^2)))
    end
    an = weighted([u(0, r) ~ 1, Dr(u(t, 0)) ~ 0, u(t, 1) ~ 0])
    nat = only(an.natural)
    @test nat.piece == BoxFace(Symbolics.unwrap(r), :lo) && nat.eq == 1 && symeq(nat.Φ, 0)
    @test only(an.essential).piece == BoxFace(Symbolics.unwrap(r), :hi)
    @test_throws "the conormal flux of equation 1 vanishes identically" weighted(
        [u(0, r) ~ 1, Dr(u(t, 0)) ~ 1, u(t, 1) ~ 0]
    )
end

@testset "over-specified faces" begin
    @parameters t x
    @variables u(..)
    Dt, Dx = Differential(t), Differential(x)
    domains = [t ∈ Interval(0.0, 1.0), x ∈ Interval(0.0, 1.0)]
    s = StateSlots(1, 1, 0)
    function classify(bcs)
        v, bmap = boxsetup([Dt(u(t, x)) ~ Dx(Dx(u(t, x)))], bcs, domains, [t, x], [u(t, x)], [], t)
        return boundary_conditions(bmap, SlotContext(v), s, [[-s.G[1]]]; mass = [s.Ud[1]])
    end
    @test_throws "both give the natural condition of equation 1" classify(
        [u(0, x) ~ 0, Dx(u(t, 0)) ~ 1, Dx(u(t, 0)) ~ 2]
    )
    @test_throws "impose both an essential and a natural condition on field 1" classify(
        [u(0, x) ~ 0, u(t, 0) ~ 1, Dx(u(t, 0)) ~ 2]
    )
    @test_throws "both fix field 1" classify([u(0, x) ~ 0, u(t, 0) ~ 1, u(t, 0) ~ 2])
    # periodic faces take no other condition on the field
    @test_throws "field 1 is periodic on BoxFace(x, :lo)" classify(
        [u(0, x) ~ 0, u(t, 0) ~ u(t, 1), Dx(u(t, 0)) ~ 1]
    )
    @test_throws "field 1 is periodic on BoxFace(x, :hi)" classify(
        [u(0, x) ~ 0, u(t, 0) ~ u(t, 1), u(t, 1) ~ 3]
    )
    # one condition per face is fine
    ess, nat, _ = classify([u(0, x) ~ 0, u(t, 0) ~ 1, Dx(u(t, 1)) ~ 2])
    @test length(ess) == 1 && length(nat) == 1
end

@testset "pairing (spec §5.4)" begin
    s = StateSlots(3, 2, 0)
    Ud = s.Ud
    @test pair_equations([Ud[2], 2Ud[1], Ud[3]], s) == [2, 1, 3]
    # algebraic row: needs an explicit pairing; positional is not the fallback
    @test_throws "equation 2 has no time derivative" pair_equations([Ud[1], 0, Ud[2]], s)
    @test pair_equations([Ud[1], 0, Ud[2]], s; pairing = Dict(2 => 3)) == [1, 3, 2]
    @test_throws "time derivatives of fields [1, 2]" pair_equations([Ud[1] + Ud[2], Ud[3]], s)
    @test_throws "but there are 3 fields" pair_equations([Ud[1], 0, Ud[2]], s; pairing = Dict(2 => 4))
    @test_throws "names equation 5" pair_equations([Ud[1], Ud[2], Ud[3]], s; pairing = Dict(5 => 1))
    @test_throws "field 1 is paired with equations 1 and 2" pair_equations(
        [Ud[1], Ud[2]], s; pairing = Dict(2 => 1)
    )
end
