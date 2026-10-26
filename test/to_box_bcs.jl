using PDEBase
using PDEBase: BoxFace, Piece, location, boxface
using ModelingToolkit
using ModelingToolkit: get_bcs, get_eqs, get_domain, get_ivs, get_dvs, get_ps
using SciMLBase
using Symbolics: unwrap
using SymbolicUtils: iscall, operation, arguments, substitute
using Test
import DomainSets
using DomainSets: Interval, Rectangle, (..)
import PDESystemLibrary

struct BoxLoweringDisc <: PDEBase.AbstractBoxDiscretization
    time::Any
    seen::Base.RefValue{Any}
end
BoxLoweringDisc(t) = BoxLoweringDisc(t, Ref{Any}(nothing))
struct MeshLoweringDisc <: PDEBase.AbstractMeshDiscretization
    time::Any
end
PDEBase.get_time(d::Union{BoxLoweringDisc, MeshLoweringDisc}) = d.time

# Stops `symbolic_discretize` at the first hook that sees the processed system.
struct PipelineStop <: Exception end
function PDEBase.interface_errors(sys::PDESystem, ::VariableMap, d::BoxLoweringDisc)
    d.seen[] = sys
    throw(PipelineStop())
end

struct LoweringMesh end
PDEBase.meshdim(::LoweringMesh) = 2
PDEBase.facetset_names(::LoweringMesh) = (:inlet,)

function errmsg(f)
    try
        f()
    catch e
        e isa ArgumentError || rethrow()
        return e.msg
    end
    return nothing
end

lowered_bcs(sys; kwargs...) = get_bcs(to_box_bcs(sys; kwargs...))
bcs_equal(a, b) = length(a) == length(b) && all(isequal.(a, b))

@parameters t x y k α
@variables u(..)
Dt, Dx, Dy = Differential(t), Differential(x), Differential(y)
heat = Dt(u(t, x, y)) ~ Dx(Dx(u(t, x, y))) + Dy(Dy(u(t, x, y)))
rect = Rectangle(0.0 .. 1.0, 0.0 .. 2.0)
Γ_left, Γ_right, Γ_bottom, Γ_top = (boundarypiece(rect, n) for n in (:left, :right, :bottom, :top))
Dn_left, Dn_right, Dn_bottom, Dn_top = NormalDerivative.((Γ_left, Γ_right, Γ_bottom, Γ_top))

# The section 9 lowering table of the FiniteElementMethod design, each row as
# (trace form, hand-written endpoint form).
table_rows = [
    (u(t, Γ_left) ~ sin(t), u(t, 0, y) ~ sin(t)),
    (u(t, Γ_top) ~ x, u(t, x, 2) ~ x),
    (Dn_left(u(t, Γ_left)) ~ 1, -Dx(u(t, 0, y)) ~ 1),
    (Dn_right(u(t, Γ_right)) ~ 0, Dx(u(t, 1, y)) ~ 0),
    (Dn_bottom(u(t, Γ_bottom)) ~ x, -Dy(u(t, x, 0)) ~ x),
    (
        k * Dn_right(u(t, Γ_right)) + α * u(t, Γ_right) ~ α * 0.3,
        k * Dx(u(t, 1, y)) + α * u(t, 1, y) ~ α * 0.3,
    ),
    (k * Dn_left(u(t, Γ_left)) + α * u(t, Γ_left) ~ 0, -k * Dx(u(t, 0, y)) + α * u(t, 0, y) ~ 0),
    (u(t, Γ_left) ~ u(t, Γ_right), u(t, 0, y) ~ u(t, 1, y)),
    (u(t, Γ_bottom) ~ u(t, Γ_top), u(t, x, 0) ~ u(t, x, 2)),
    (Dx(u(t, Γ_right)) ~ 2, Dx(u(t, 1, y)) ~ 2),
    (Dy(u(t, Γ_top)) ~ 0, Dy(u(t, x, 2)) ~ 0),
]

box_domains = [(x, y) ∈ rect, t ∈ Interval(0.0, 1.0)]
split_domains = [t ∈ Interval(0.0, 1.0), x ∈ Interval(0.0, 1.0), y ∈ Interval(0.0, 2.0)]
heat_system(bcs, domains; kwargs...) =
    PDESystem([heat], [u(0, x, y) ~ 0; bcs], domains, [t, x, y], [u(t, x, y)]; name = :heat, kwargs...)

@testset "section 9 lowering table" begin
    for domains in (box_domains, split_domains), (trace, endpoint) in table_rows
        lowered = lowered_bcs(heat_system([trace], domains))
        @test isequal(lowered[1], u(0, x, y) ~ 0)
        @test isequal(lowered[2], endpoint)
        @test isequal(lowered_bcs(heat_system([trace], domains); time = t)[2], endpoint)
    end
    # All rows at once, as one system would carry them.
    for domains in (box_domains, split_domains)
        @test bcs_equal(
            lowered_bcs(heat_system(first.(table_rows), domains)), [u(0, x, y) ~ 0; last.(table_rows)]
        )
    end

    @testset "endpoint values" begin
        # The endpoint comes from the domain, as a Float64, and is isequal to the integer
        # endpoint a user writes.
        lowered = only(filter(eq -> !isequal(eq, u(0, x, y) ~ 0), lowered_bcs(heat_system([u(t, Γ_left) ~ 0], box_domains))))
        @test isequal(lowered, u(t, 0.0, y) ~ 0)
        @test isequal(lowered, u(t, 0, y) ~ 0)
    end

    @testset "1D interval" begin
        @variables w(..)
        line = Interval(0.0, π)
        Γ_lo, Γ_hi = boundarypieces(line)
        @test (Γ_lo.name, Γ_hi.name) == (:x1_lo, :x1_hi)
        @test isequal(Γ_lo.geometry, DomainSets.Point(0.0))
        @test boundarypiece(line, :x1_hi) == Γ_hi
        Dn_lo, Dn_hi = NormalDerivative(Γ_lo), NormalDerivative(Γ_hi)
        @named line_sys = PDESystem(
            [Dt(w(t, x)) ~ Dx(Dx(w(t, x)))],
            [
                w(0, x) ~ cos(x), w(t, Γ_lo) ~ exp(-t), Dn_hi(w(t, Γ_hi)) ~ 0,
                k * Dn_lo(w(t, Γ_lo)) + α * w(t, Γ_lo) ~ 1, w(t, Γ_lo) ~ w(t, Γ_hi),
                Dx(Dx(w(t, Γ_hi))) ~ 0,
            ],
            [t ∈ Interval(0.0, 1.0), x ∈ line], [t, x], [w(t, x)]
        )
        # The endpoint is the domain's own value: `Interval(0.0, π)` holds the Float64 π,
        # which is not `isequal` to the irrational `π`.
        Lπ = Float64(π)
        @test !isequal(Lπ, π)
        @test bcs_equal(
            lowered_bcs(line_sys),
            [
                w(0, x) ~ cos(x), w(t, 0) ~ exp(-t), Dx(w(t, Lπ)) ~ 0,
                -k * Dx(w(t, 0)) + α * w(t, 0) ~ 1, w(t, 0) ~ w(t, Lπ), Dx(Dx(w(t, Lπ))) ~ 0,
            ]
        )
    end

    @testset "argument orders and steady systems" begin
        @variables s(..) b(..)
        @named steady = PDESystem(
            [Dx(Dx(s(x, y))) + Dy(Dy(s(x, y))) ~ 0],
            [s(Γ_left) ~ 0, Dn_top(s(Γ_top)) ~ 1], [(x, y) ∈ rect], [x, y], [s(x, y)]
        )
        @test bcs_equal(lowered_bcs(steady), [s(0, y) ~ 0, Dy(s(x, 2)) ~ 1])
        # Time after the spatial arguments, as in the Brusselator.
        @named late = PDESystem(
            [Dt(b(x, y, t)) ~ Dx(Dx(b(x, y, t)))],
            [b(x, y, 0) ~ 0, b(Γ_left, t) ~ b(Γ_right, t), Dn_bottom(b(Γ_bottom, t)) ~ t],
            [(x, y) ∈ rect, t ∈ Interval(0.0, 1.0)], [x, y, t], [b(x, y, t)]
        )
        @test bcs_equal(
            lowered_bcs(late), [b(x, y, 0) ~ 0, b(0, y, t) ~ b(1, y, t), -Dy(b(x, 0, t)) ~ t]
        )
        # A space-time box: the spatial face is an endpoint of the x axis only.
        @variables p(..)
        st = Rectangle(0.0 .. 1.0, 0.0 .. 3.0)
        Γ_x = boundarypiece(Interval(0.0, 3.0), :x1_hi)
        @named spacetime = PDESystem(
            [Dt(p(t, x)) ~ Dx(Dx(p(t, x)))], [p(0, x) ~ 0, NormalDerivative(Γ_x)(p(t, Γ_x)) ~ 0],
            [(t, x) ∈ st], [t, x], [p(t, x)]
        )
        @test bcs_equal(lowered_bcs(spacetime), [p(0, x) ~ 0, Dx(p(t, 3)) ~ 0])
    end
end

@testset "box axes follow the domain declaration" begin
    @variables q(..)
    I = Interval(0.0, 1.0)
    # `q(t, y)` comes first, so `indvars` starts with `y`; the box axes do not.
    coupled = [Dt(q(t, y)) ~ Dy(Dy(q(t, y))), Dt(u(t, x, y)) ~ Dx(Dx(u(t, x, y))) + q(t, y)]
    for (Ix, Iy) in ((I, I), (I, Interval(0.0, 2.0)))
        Γ = boundarypiece(Rectangle(Ix, Iy), :left)
        sys = PDESystem(
            coupled, [u(0, x, y) ~ 0, q(0, y) ~ 0, u(t, Γ) ~ 1], [t ∈ I, x ∈ Ix, y ∈ Iy],
            [t, x, y], [q(t, y), u(t, x, y)]; name = :coupled
        )
        @test isequal(lowered_bcs(sys)[3], u(t, 0.0, y) ~ 1)
        Γ_top = boundarypiece(Rectangle(Ix, Iy), :top)
        sys = PDESystem(
            coupled, [u(0, x, y) ~ 0, q(0, y) ~ 0, u(t, Γ_top) ~ 1], [t ∈ I, x ∈ Ix, y ∈ Iy],
            [t, x, y], [q(t, y), u(t, x, y)]; name = :coupled
        )
        @test isequal(lowered_bcs(sys)[3], u(t, x, Iy.right) ~ 1)
    end
    # A space-time pairing, with the variable declared in another order than the pairing.
    @variables r(..)
    Γ = boundarypiece(Rectangle(I, I), :left)
    @named spacetime = PDESystem(
        [Dt(r(t, y, x)) ~ Dx(Dx(r(t, y, x)))], [r(0, y, x) ~ 0, r(t, Γ) ~ 1],
        [(t, x, y) ∈ DomainSets.ProductDomain(I, I, I)], [t, x, y], [r(t, y, x)]
    )
    @test isequal(lowered_bcs(spacetime)[2], r(t, y, 0.0) ~ 1)
end

@testset "time inference" begin
    @variables w(..) m(..)
    # Time between the spatial arguments: the piece replaces x and y, which are not adjacent.
    @named middle = PDESystem(
        [Dt(m(x, t, y)) ~ Dx(Dx(m(x, t, y)))],
        [m(x, 0, y) ~ 0, m(Γ_left, t) ~ 1, Dn_top(m(Γ_top, t)) ~ 0], box_domains, [t, x, y], [m(x, t, y)]
    )
    expected = [m(x, 0, y) ~ 0, m(0, t, y) ~ 1, Dy(m(x, t, 2)) ~ 0]
    @test bcs_equal(lowered_bcs(middle), expected)
    @test bcs_equal(lowered_bcs(middle; time = t), expected)
    # Two layouts in one system.
    @named mixed = PDESystem(
        [Dt(m(x, t, y)) ~ Dx(Dx(m(x, t, y))), Dt(w(t, x, y)) ~ Dy(Dy(w(t, x, y)))],
        [m(x, 0, y) ~ 0, w(0, x, y) ~ 0, m(Γ_left, t) ~ w(t, Γ_left), w(t, Γ_right) ~ 1],
        box_domains, [t, x, y], [m(x, t, y), w(t, x, y)]
    )
    expected = [m(x, 0, y) ~ 0, w(0, x, y) ~ 0, m(0, t, y) ~ w(t, 0, y), w(t, 1, y) ~ 1]
    @test bcs_equal(lowered_bcs(mixed), expected)
    @test bcs_equal(lowered_bcs(mixed; time = t), expected)
    # A given time is used without inference: on equal intervals `u(Γ)` cannot tell `t`
    # from `x`, and with `time = t` the classifier reports the missing argument.
    line = Interval(0.0, 1.0)
    Γ_lo = boundarypiece(line, :x1_lo)
    @named ambiguous = PDESystem(
        [Dt(w(t, x)) ~ Dx(Dx(w(t, x)))], [w(0, x) ~ 0, w(Γ_lo) ~ 0], [t ∈ line, x ∈ line], [t, x], [w(t, x)]
    )
    @test occursin("pass it as to_box_bcs(pdesys; time = t)", errmsg(() -> to_box_bcs(ambiguous)))
    # A normal derivative of an endpoint-form operand carries the only piece; it decides.
    @named dn_only = PDESystem(
        [Dt(w(t, x)) ~ Dx(Dx(w(t, x)))], [w(0, x) ~ 0, NormalDerivative(Γ_lo)(w(t, 0)) ~ 0],
        [t ∈ line, x ∈ line], [t, x], [w(t, x)]
    )
    @test bcs_equal(lowered_bcs(dn_only), [w(0, x) ~ 0, -Dx(w(t, 0)) ~ 0])
    @test occursin("trace is missing the time argument", errmsg(() -> to_box_bcs(ambiguous; time = t)))
end

@testset "normal derivative of an expression" begin
    @variables w(..)
    two = heat_system([NormalDerivative(Γ_left)(2 * u(t, Γ_left)) ~ 0], box_domains)
    @test isequal(lowered_bcs(two)[2], -Dx(2 * u(t, 0.0, y)) ~ 0)
    @named sum_sys = PDESystem(
        [heat, Dt(w(t, x, y)) ~ Dx(Dx(w(t, x, y)))],
        [u(0, x, y) ~ 0, w(0, x, y) ~ 0, Dn_right(u(t, Γ_right) + k * w(t, Γ_right)) ~ 1],
        box_domains, [t, x, y], [u(t, x, y), w(t, x, y)]
    )
    @test isequal(lowered_bcs(sum_sys)[3], Dx(u(t, 1.0, y) + k * w(t, 1.0, y)) ~ 1)
end

@testset "Pair conditions" begin
    lowered = lowered_bcs(heat_system([Dn_right(u(t, Γ_right)) => 0], box_domains))
    @test lowered[2] isa Pair
    @test isequal(lowered[2].first, Dx(u(t, 1, y)))
    @test isequal(lowered[2].second, 0)
end

@testset "rejections" begin
    mesh_sys = PDESystem(
        [heat], [u(0, x, y) ~ 0], [(x, y) ∈ MeshDomain(LoweringMesh()), t ∈ Interval(0.0, 1.0)],
        [t, x, y], [u(t, x, y)]; name = :mesh
    )
    mesh_msg = "MethodOfLines needs an interval or box domain; this system uses a MeshDomain, use FEMDiscretization"
    @test errmsg(() -> to_box_bcs(mesh_sys)) == mesh_msg
    @test errmsg(() -> PDEBase.lower_boundaries(mesh_sys, BoxLoweringDisc(t))) == mesh_msg
    # A mesh discretization leaves the system alone.
    @test PDEBase.lower_boundaries(mesh_sys, MeshLoweringDisc(t)) === mesh_sys

    disk_sys = PDESystem(
        [heat], [u(0, x, y) ~ 0], [(x, y) ∈ DomainSets.UnitDisk(), t ∈ Interval(0.0, 1.0)],
        [t, x, y], [u(t, x, y)]; name = :disk
    )
    msg = errmsg(() -> to_box_bcs(disk_sys))
    @test startswith(msg, "MethodOfLines needs an interval or box domain; this system uses a ")
    @test occursin("Ball", msg)
    @test endswith(msg, ", use FEMDiscretization")

    # The trace rule: `u(Γ)` on a time-dependent variable.
    missing_time = heat_system([u(Γ_left) ~ 0], box_domains)
    @test occursin("trace is missing the time argument", errmsg(() -> to_box_bcs(missing_time; time = t)))
    @test occursin(
        "trace is missing the time argument",
        errmsg(() -> PDEBase.lower_boundaries(missing_time, BoxLoweringDisc(t)))
    )
    msg = errmsg(() -> to_box_bcs(missing_time))
    @test occursin("trace is missing the time argument", msg)
    # A number in the place of time is classified, not mistaken for a missing argument.
    msg = errmsg(() -> to_box_bcs(heat_system([u(0.5, Γ_left) ~ 0], box_domains)))
    @test occursin("0.5 in the place of t", msg)
    # A trace on a face whose coordinate the variable does not depend on.
    @variables q(..)
    @named partial = PDESystem(
        [heat, Dt(q(t, y)) ~ Dy(Dy(q(t, y)))], [u(0, x, y) ~ 0, q(0, y) ~ 0, q(t, Γ_left) ~ 1],
        box_domains, [t, x, y], [u(t, x, y), q(t, y)]
    )
    @test occursin("does not depend on x", errmsg(() -> to_box_bcs(partial)))
    # Array and complex dependent variables.
    line = Interval(0.0, 1.0)
    Γ_lo = boundarypiece(line, :x1_lo)
    @variables a(..)[1:2] ψ(..)::Complex
    @named arr = PDESystem(
        [Dt(a(t, x)[i]) ~ Dx(Dx(a(t, x)[i])) for i in 1:2], [a(0, x)[1] ~ 0, a(t, Γ_lo)[1] ~ 1],
        [t ∈ line, x ∈ line], [t, x], [a(t, x)]
    )
    @test occursin("traces on array or complex dependent variables are not supported yet", errmsg(() -> to_box_bcs(arr)))
    @named cplx = PDESystem(
        [Dt(ψ(t, x)) ~ im * Dx(Dx(ψ(t, x)))], [ψ(0, x) ~ 0, ψ(t, Γ_lo) ~ 0],
        [t ∈ line, x ∈ line], [t, x], [ψ(t, x)]
    )
    @test occursin("trace of a complex dependent variable", errmsg(() -> to_box_bcs(cplx)))
    # A periodic pair that is not the two faces of one axis.
    msg = errmsg(() -> to_box_bcs(heat_system([u(t, Γ_left) ~ u(t, Γ_bottom)], box_domains)))
    @test occursin("relates the locations", msg)
    # An operator on one piece applied to a trace on another.
    msg = errmsg(() -> to_box_bcs(heat_system([Dn_left(u(t, Γ_right)) ~ 0], box_domains)))
    @test occursin("write the trace on the operator's piece", msg)
    # A piece of the box that is not one of its faces.
    Γ_odd = BoundaryPiece(rect, :inlet)
    msg = errmsg(() -> to_box_bcs(heat_system([u(t, Γ_odd) ~ 0], box_domains)))
    @test occursin("is not a face of the box", msg)
end

@testset "mixed styles and passthrough" begin
    mixed = [
        u(t, 0, y) ~ sin(t), Dn_right(u(t, Γ_right)) ~ 0, Dy(u(t, x, 0)) ~ 0,
        k * Dn_top(u(t, Γ_top)) + α * u(t, Γ_top) ~ α,
    ]
    expected = [
        u(0, x, y) ~ 0, u(t, 0, y) ~ sin(t), Dx(u(t, 1, y)) ~ 0, Dy(u(t, x, 0)) ~ 0,
        k * Dy(u(t, x, 2)) + α * u(t, x, 2) ~ α,
    ]
    sys = heat_system(mixed, box_domains; metadata = :meta, description = "mixed")
    lowered = to_box_bcs(sys)
    @test bcs_equal(get_bcs(lowered), expected)
    @test get_bcs(lowered)[2] === get_bcs(sys)[2]
    @test isequal(get_eqs(lowered), get_eqs(sys))
    @test get_domain(lowered) === get_domain(sys)
    @test isequal(get_ivs(lowered), get_ivs(sys))
    @test isequal(get_dvs(lowered), get_dvs(sys))
    @test nameof(lowered) === nameof(sys)
    @test ModelingToolkit.get_metadata(lowered) === :meta
    @test ModelingToolkit.get_description(lowered) == "mixed"

    @named with_ps = PDESystem(
        [Dt(u(t, x, y)) ~ k * Dx(Dx(u(t, x, y)))], [u(0, x, y) ~ 0, u(t, Γ_left) ~ 0], box_domains,
        [t, x, y], [u(t, x, y)], [k => 2.0]; initial_conditions = Dict(k => 2.0)
    )
    lowered = to_box_bcs(with_ps)
    @test isequal(get_ps(lowered), get_ps(with_ps))
    @test isequal(lowered.initial_conditions, with_ps.initial_conditions)

    endpoint = heat_system([u(t, 0, y) ~ 0, Dx(u(t, 1, y)) ~ 0], box_domains)
    @test to_box_bcs(endpoint) === endpoint
    @test PDEBase.lower_boundaries(endpoint, BoxLoweringDisc(t)) === endpoint
end

@testset "automatic lowering for box discretizations" begin
    function seen_by_pipeline(sys)
        disc = BoxLoweringDisc(t)
        @test_throws PipelineStop SciMLBase.symbolic_discretize(sys, disc)
        return disc.seen[]
    end
    trace_sys = heat_system(first.(table_rows[[1, 4, 6, 10]]), box_domains)
    endpoint_sys = heat_system(last.(table_rows[[1, 4, 6, 10]]), box_domains)
    @test bcs_equal(get_bcs(seen_by_pipeline(trace_sys)), get_bcs(seen_by_pipeline(endpoint_sys)))
    @test PDEBase.lower_boundaries(trace_sys, MeshLoweringDisc(t)) === trace_sys
end

### Library round trip (FiniteElementMethod design section 10.5) ###

# The independent variable named `t` is time in every library system.
library_time(sys) = only(filter(iv -> nameof(unwrap(iv)) === :t, get_ivs(sys)))

# `sys` with its spatial interval pairings combined into one box pairing in two and more
# dimensions; one-dimensional systems keep their interval.
function boxed(sys, t)
    space = filter(p -> !isequal(p.variables, t), get_domain(sys))
    length(space) == 1 && return sys
    coords = Tuple(p.variables for p in space)
    domains = [coords ∈ DomainSets.ProductDomain((p.domain for p in space)...); filter(p -> isequal(p.variables, t), get_domain(sys))]
    return PDESystem(
        get_eqs(sys), get_bcs(sys), domains, get_ivs(sys), get_dvs(sys), get_ps(sys); name = nameof(sys)
    )
end

# The trace spelling of an endpoint-form condition: endpoint calls become traces on the
# face's piece, and a first derivative across the face becomes the normal derivative.
function to_trace_form(v, pieces, eq)
    rules = Dict{Any, Any}()
    function visit(ex)
        ex = unwrap(ex)
        iscall(ex) || return
        op = operation(ex)
        if haskey(v.args, op)
            face = location(v, classify_call(v, ex)[2])
            if face isa BoxFace && haskey(pieces, face)
                kept = [a for (a, d) in zip(arguments(ex), v.args[op]) if isequal(d, v.time)]
                spatial = findfirst(d -> !isequal(d, v.time), v.args[op])
                rules[ex] = op(insert!(kept, spatial, pieces[face])...)
            end
            return
        elseif op isa Differential && op.order == 1
            call = unwrap(only(arguments(ex)))
            if iscall(call) && haskey(v.args, operation(call))
                face = location(v, classify_call(v, call)[2])
                if face isa BoxFace && isequal(face.iv, op.x) && haskey(pieces, face)
                    visit(call)
                    Dn = NormalDerivative(pieces[face])(rules[call])
                    rules[ex] = face.side === :lo ? -Dn : Dn
                    return
                end
            end
        end
        return foreach(visit, arguments(ex))
    end
    visit(eq.lhs)
    visit(eq.rhs)
    return substitute_traces(eq, rules)
end

@testset "PDESystemLibrary round trip" begin
    library = PDESystemLibrary.all_systems
    @test length(library) == 427
    ntraces = nnormal = 0
    for sys in library
        t_ = library_time(sys)
        box = boxed(sys, t_)
        v = VariableMap(box, BoxLoweringDisc(t_))
        pieces = Dict(boxface(v, p) => p for p in boundarypieces(spatial_domain(v)))
        traced = map(get_bcs(box)) do bc
            classify_condition(v, bc).kind === :initial ? bc : to_trace_form(v, pieces, bc)
        end
        ntraces += count(bc -> PDEBase._has_trace(bc), traced)
        nnormal += count(bc -> occursin("Dn_", string(bc)), traced)
        # Every boundary condition of the library is rewritten.
        @test count(bc -> classify_condition(v, bc).kind !== :initial, get_bcs(box)) ==
            count(bc -> PDEBase._has_trace(bc), traced)
        # The classifier maps both spellings to the same faces.
        @test all(zip(get_bcs(box), traced)) do (bc, tr)
            isequal(classify_condition(v, bc).locations, classify_condition(v, tr).locations)
        end
        trace_sys = PDESystem(
            get_eqs(box), traced, get_domain(box), get_ivs(box), get_dvs(box), get_ps(box); name = nameof(sys)
        )
        @test bcs_equal(get_bcs(to_box_bcs(trace_sys)), get_bcs(sys))
    end
    @test ntraces >= length(library)
    # The Neumann conditions of `pdesys` and `sph` go through the normal derivative.
    @test nnormal >= 3
end
