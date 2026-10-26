using PDEBase
using PDEBase: parse_bcs
using ModelingToolkit
using Symbolics: unwrap
using Test
import DomainSets
using DomainSets: Interval, Rectangle, (..)

include(joinpath(@__DIR__, "shared", "legacy_parse_bcs.jl"))

struct GenericBCDiscretization <: PDEBase.AbstractEquationSystemDiscretization
    time::Any
end
struct BoxBCDiscretization <: PDEBase.AbstractBoxDiscretization
    time::Any
end
PDEBase.get_time(disc::Union{GenericBCDiscretization, BoxBCDiscretization}) = disc.time

# The steps `symbolic_discretize` runs before `parse_bcs`.
function prepare(pdesys, disc)
    pdesys = PDEBase.lower_boundaries(pdesys, disc)
    pdesys, _ = PDEBase.handle_complex(pdesys)
    PDEBase.cardinalize_eqs!(pdesys)
    pdesys, replaced = PDEBase.make_pdesys_compatible(pdesys)
    v = VariableMap(pdesys, disc; replaced_vars = replaced)
    bcs = PDEBase.get_bcs(pdesys)
    orders = Dict(map(x -> x => d_orders(x, bcs, v), all_ivs(v)))
    return bcs, v, orders
end

function same_boundary(a, b; fields = fieldnames(typeof(a)))
    typeof(a) === typeof(b) || return false
    return all(f -> isequal(getfield(a, f), getfield(b, f)), fields)
end

# Same keys, in the same iteration order when `ordered`, and the same boundaries field by
# field. Maps of two different systems need not share the order of their spatial variables.
function same_map(a, b; ordered = true, kwargs...)
    samekeys(p, q) = ordered ? isequal(collect(keys(p)), collect(keys(q))) : issetequal(keys(p), keys(q))
    samekeys(a, b) || return false
    for u in keys(a)
        samekeys(a[u], b[u]) || return false
        for x in keys(a[u])
            length(a[u][x]) == length(b[u][x]) || return false
            all(same_boundary(p, q; kwargs...) for (p, q) in zip(a[u][x], b[u][x])) || return false
        end
    end
    return true
end

# Runs the old matcher and every new path on the same input and checks they agree, and
# that the shared constructor helper reproduces the old inline field computation.
function check_equivalence(pdesys, t)
    bcs, v, orders = prepare(pdesys, BoxBCDiscretization(t))
    legacy = legacy_parse_bcs(bcs, v, orders)
    @test same_map(parse_bcs(bcs, v, orders), legacy)
    @test same_map(parse_bcs(bcs, v, orders, GenericBCDiscretization(t)), legacy)
    box = parse_bcs(bcs, v, orders, BoxBCDiscretization(t))
    @test same_map(box, legacy)
    for b in flatten_vardict(legacy)
        hasfield(typeof(b), :depvars) || continue
        @test isequal((b.depvars, b.indvars), legacy_boundary_vars(b.eq, v.time, v))
    end
    @test length(flatten_vardict(box)) >= length(bcs)
    return box
end

@parameters t x y x1 x2
@variables u(..) w(..) c1(..) c2(..)
Dt, Dx, Dy = Differential(t), Differential(x), Differential(y)
Dxx, Dyy = Dx^2, Dy^2
Dx1, Dx2 = Differential(x1), Differential(x2)

heat1d(bcs; domain = x ∈ Interval(0.0, 1.0)) = PDESystem(
    [Dt(u(t, x)) ~ Dxx(u(t, x))], bcs, [t ∈ Interval(0.0, 1.0), domain], [t, x], [u(t, x)];
    name = :heat1d
)
heat2d(bcs, domains) = PDESystem(
    [Dt(u(t, x, y)) ~ Dxx(u(t, x, y)) + Dyy(u(t, x, y))], bcs, domains, [t, x, y], [u(t, x, y)];
    name = :heat2d
)

rect = Rectangle(0.0 .. 1.0, 0.0 .. 2.0)
Γ_left, Γ_right, Γ_bottom, Γ_top = (boundarypiece(rect, n) for n in (:left, :right, :bottom, :top))
rect_domains = [(x, y) ∈ rect, t ∈ Interval(0.0, 1.0)]

@testset "parse_bcs box method matches the substitution matcher" begin
    @testset "Dirichlet, Neumann, Robin and time-dependent 1D" begin
        check_equivalence(heat1d([u(0, x) ~ 0, u(t, 0) ~ 0, u(t, 1) ~ 0]), t)
        check_equivalence(heat1d([u(0, x) ~ sin(pi * x), Dx(u(t, 0)) ~ 0, Dx(u(t, 1)) + 2u(t, 1) ~ sin(t)]), t)
        check_equivalence(heat1d([u(0, x) ~ x, 3u(t, 0) - Dx(u(t, 0)) ~ exp(-t), Dt(u(t, 1)) ~ cos(t)]), t)
        check_equivalence(heat1d([u(0.0, x) ~ x^2, u(t, 0.0) ~ t, Dx(u(t, 1.0)) ~ 2]), t)
    end

    @testset "wave equation, second-order initial conditions" begin
        wave = PDESystem(
            [Dt(Dt(u(t, x))) ~ Dxx(u(t, x))],
            [u(0, x) ~ sin(pi * x), Dt(u(0, x)) ~ 0, u(t, 0) ~ 0, u(t, 1) ~ 0],
            [t ∈ Interval(0.0, 1.0), x ∈ Interval(0.0, 1.0)], [t, x], [u(t, x)]; name = :wave
        )
        box = check_equivalence(wave, t)
        uop = operation(unwrap(u(t, x)))
        @test sort([b.order for b in box[uop][unwrap(t)]]) == [0, 1]
    end

    @testset "higher-order boundary conditions" begin
        beam = PDESystem(
            [Dt(u(t, x)) ~ -Dx(Dx(Dx(Dx(u(t, x)))))],
            [u(0, x) ~ 0, u(t, 0) ~ 0, Dxx(u(t, 0)) ~ 0, Dx(u(t, 1)) ~ 1, (Dx^3)(u(t, 1)) ~ 0],
            [t ∈ Interval(0.0, 1.0), x ∈ Interval(0.0, 1.0)], [t, x], [u(t, x)]; name = :beam
        )
        check_equivalence(beam, t)
    end

    @testset "2D, separate intervals and a tuple pairing on a Rectangle" begin
        bcs = [u(0, x, y) ~ 0, u(t, 0, y) ~ 0, u(t, 1, y) ~ y, Dy(u(t, x, 0)) ~ 0, u(t, x, 2) ~ 0]
        check_equivalence(
            heat2d(bcs, [t ∈ Interval(0.0, 1.0), x ∈ Interval(0.0, 1.0), y ∈ Interval(0.0, 2.0)]), t
        )
        check_equivalence(heat2d(bcs, rect_domains), t)
    end

    @testset "steady problem" begin
        poisson = PDESystem(
            [Dxx(u(x, y)) + Dyy(u(x, y)) ~ 1],
            [u(0, y) ~ 0, u(1, y) ~ 0, Dy(u(x, 0)) ~ 0, u(x, 2) ~ x],
            [x ∈ Interval(0.0, 1.0), y ∈ Interval(0.0, 2.0)], [x, y], [u(x, y)]; name = :poisson
        )
        check_equivalence(poisson, nothing)
    end

    @testset "periodic" begin
        box = check_equivalence(heat1d([u(0, x) ~ sin(2pi * x), u(t, 0) ~ u(t, 1)]), t)
        uop = operation(unwrap(u(t, x)))
        @test all(b -> b isa InterfaceBoundary, box[uop][unwrap(x)])
        @test length(box[uop][unwrap(x)]) == 2
    end

    @testset "Brusselator, time-last declaration" begin
        Dt_, Dx_, Dy_ = Differential(t), Differential(x)^2, Differential(y)^2
        eqs = [
            Dt_(u(x, y, t)) ~ 1.0 + w(x, y, t) * u(x, y, t)^2 - 4.4 * u(x, y, t) +
                10.0 * (Dx_(u(x, y, t)) + Dy_(u(x, y, t))),
            Dt_(w(x, y, t)) ~ 3.4 * u(x, y, t) - w(x, y, t) * u(x, y, t)^2 +
                10.0 * (Dx_(w(x, y, t)) + Dy_(w(x, y, t))),
        ]
        bcs = [
            u(x, y, 0) ~ 22(y * (1 - y))^(3 / 2), u(0, y, t) ~ u(1, y, t), u(x, 0, t) ~ u(x, 1, t),
            w(x, y, 0) ~ 27(x * (1 - x))^(3 / 2), w(0, y, t) ~ w(1, y, t), w(x, 0, t) ~ w(x, 1, t),
        ]
        domains = [x ∈ Interval(0.0, 1.0), y ∈ Interval(0.0, 1.0), t ∈ Interval(0.0, 11.5)]
        @named bruss = PDESystem(eqs, bcs, domains, [x, y, t], [u(x, y, t), w(x, y, t)])
        check_equivalence(bruss, t)
    end

    @testset "interface between two variables" begin
        eqs = [Dt(c1(t, x1)) ~ Dx1(Dx1(c1(t, x1))), Dt(c2(t, x2)) ~ Dx2(Dx2(c2(t, x2)))]
        domains = [t ∈ Interval(0.0, 1.0), x1 ∈ Interval(0.0, 0.5), x2 ∈ Interval(0.5, 1.0)]
        bcs = [
            c1(0, x1) ~ -x1 * (x1 - 1) * sin(x1), c2(0, x2) ~ x2 * (x2 - 1) * sin(x2),
            c1(t, 0) ~ 0, c1(t, 0.5) ~ c2(t, 0.5), c2(t, 1) ~ 0,
        ]
        @named iface = PDESystem(eqs, bcs, domains, [t, x1, x2], [c1(t, x1), c2(t, x2)])
        box = check_equivalence(iface, t)
        @test count(b -> b isa InterfaceBoundary, flatten_vardict(box)) == 2

        flux = [bcs; Dx1(c1(t, 0.5)) ~ Dx2(c2(t, 0.5))]
        @named fluxsys = PDESystem(eqs, flux, domains, [t, x1, x2], [c1(t, x1), c2(t, x2)])
        box = check_equivalence(fluxsys, t)
        @test count(b -> b isa HigherOrderInterfaceBoundary, flatten_vardict(box)) == 1
    end

    @testset "two variables on one face" begin
        eqs = [Dt(u(t, x)) ~ Dxx(u(t, x)) + w(t, x), Dt(w(t, x)) ~ Dxx(w(t, x))]
        bcs = [
            u(0, x) ~ 0, w(0, x) ~ 0, u(t, 0) ~ w(t, 0), w(t, 0) ~ 1,
            u(t, 1) + w(t, 1) ~ 0, Dx(w(t, 1)) ~ 0,
        ]
        @named coupled = PDESystem(
            eqs, bcs, [t ∈ Interval(0.0, 1.0), x ∈ Interval(0.0, 1.0)], [t, x], [u(t, x), w(t, x)]
        )
        check_equivalence(coupled, t)
    end

    @testset "complex boundary conditions" begin
        @variables ψ(..) χ(..)
        complex_pair = PDESystem(
            [im * Dt(ψ(t, x)) ~ (1 + 2im) * Dxx(χ(t, x)), im * Dt(χ(t, x)) ~ ψ(t, x) * χ(t, x)],
            [ψ(t, 0) ~ 0, ψ(t, 1) ~ 0, χ(t, 0) ~ 0, χ(t, 1) ~ 0],
            [t ∈ (0, 1), x ∈ (0, 1)], [t, x], [ψ(t, x), χ(t, x)]; name = :complex_pair
        )
        check_equivalence(complex_pair, t)
    end
end

@testset "parse_bcs box method: endpoint tolerance" begin
    for (domain, given, exact) in ((Interval(0.0, π), π, Float64(π)), (Interval(0.0, 0.3), 0.1 + 0.2, 0.3))
        for condition in (c -> u(t, c) ~ 1, c -> Dx(u(t, c)) ~ 0)
            bcs_given = [u(0, x) ~ 0, u(t, 0) ~ 0, condition(given)]
            bcs_exact = [u(0, x) ~ 0, u(t, 0) ~ 0, condition(exact)]
            bcs, v, orders = prepare(heat1d(bcs_given; domain = x ∈ domain), BoxBCDiscretization(t))
            # The substitution matcher needs the endpoint exactly.
            @test_throws AssertionError parse_bcs(bcs, v, orders)
            box = parse_bcs(bcs, v, orders, BoxBCDiscretization(t))
            # The near spelling is rewritten to the endpoint, so `u`, `eq` and `indvars` agree.
            b = only(box[operation(unwrap(u(t, x)))][unwrap(x)][2:end])
            @test b isa UpperBoundary
            @test isempty(get_depvars(PDEBase.substitute_traces(b.eq.lhs, Dict(b.u => 42.0)), v.depvar_ops))
            if b.order == 0
                @test isequal(PDEBase.unwrap_const(substitute(b.eq.lhs, Dict(b.u => 42.0))), 42.0)
            end
            bcs, v, orders = prepare(heat1d(bcs_exact; domain = x ∈ domain), BoxBCDiscretization(t))
            @test same_map(box, legacy_parse_bcs(bcs, v, orders); ordered = false)
        end
    end
end

@testset "parse_bcs box method: traces on a Rectangle" begin
    Dn_right, Dn_bottom = NormalDerivative(Γ_right), NormalDerivative(Γ_bottom)
    traced = [
        u(0, x, y) ~ 0, u(t, Γ_left) ~ 0, Dn_right(u(t, Γ_right)) ~ 0,
        2Dn_bottom(u(t, Γ_bottom)) ~ 1, Dy(u(t, Γ_top)) ~ sin(t),
    ]
    endpoint = [
        u(0, x, y) ~ 0, u(t, 0, y) ~ 0, Dx(u(t, 1, y)) ~ 0,
        -2Dy(u(t, x, 0)) ~ 1, Dy(u(t, x, 2)) ~ sin(t),
    ]
    bcs, v, orders = prepare(heat2d(traced, rect_domains), BoxBCDiscretization(t))
    box = parse_bcs(bcs, v, orders, BoxBCDiscretization(t))
    bcs, v, orders = prepare(heat2d(endpoint, rect_domains), BoxBCDiscretization(t))
    # `eq` is the condition as written and `indvars` holds the trace's piece. A Robin
    # condition is left out: its `order` is that of whichever term of the sum comes first,
    # and the two spellings need not sort their terms alike.
    endpoint_map = legacy_parse_bcs(bcs, v, orders)
    @test same_map(box, endpoint_map; ordered = false, fields = (:u, :x, :depvars, :order))
    @test !same_map(box, endpoint_map; ordered = false)

    # Nested normal derivatives, alone or mixed with the face's coordinate derivative, have
    # the order of the folded endpoint spelling.
    Dn_left, Dn_top = NormalDerivative(Γ_left), NormalDerivative(Γ_top)
    nested = [
        u(0, x, y) ~ 0, Dn_right(Dn_right(u(t, Γ_right))) ~ 0, Dx(Dn_left(u(t, Γ_left))) ~ 0,
        Dn_top(Dy(u(t, Γ_top))) ~ 0, Dn_bottom(Dn_bottom(u(t, Γ_bottom))) ~ 0,
    ]
    nested_endpoint = [
        u(0, x, y) ~ 0, Dxx(u(t, 1, y)) ~ 0, -Dxx(u(t, 0, y)) ~ 0,
        Dyy(u(t, x, 2)) ~ 0, Dyy(u(t, x, 0)) ~ 0,
    ]
    bcs, v, orders = prepare(heat2d(nested, rect_domains), BoxBCDiscretization(t))
    @test orders[unwrap(x)] == [2]
    box = parse_bcs(bcs, v, orders, BoxBCDiscretization(t))
    uop = operation(unwrap(u(t, x, y)))
    @test [b.order for b in box[uop][unwrap(x)]] == [2, 2]
    @test [b.order for b in box[uop][unwrap(y)]] == [2, 2]
    bcs, v, orders = prepare(heat2d(nested_endpoint, rect_domains), BoxBCDiscretization(t))
    @test same_map(box, legacy_parse_bcs(bcs, v, orders); ordered = false, fields = (:u, :x, :depvars, :order))

    mismatched = [u(0, x, y) ~ 0, NormalDerivative(Γ_left)(u(t, Γ_right)) ~ 0]
    bcs, v, orders = prepare(heat2d(mismatched, rect_domains), BoxBCDiscretization(t))
    @test_throws ArgumentError parse_bcs(bcs, v, orders, BoxBCDiscretization(t))
end
