# Stage d (spec §5.3): put each boundary condition into slots, then classify it as essential,
# natural (matched against an equation's conormal flux) or periodic.

"""
    BoxFace(iv, side)

The `:lo` or `:hi` face of an analytic box along independent variable `iv`. Boundary key of
box conditions; stands in for the classifier's `BoxFace` (spec §4.2) until that lands.
"""
struct BoxFace
    iv::Any
    side::Symbol
end
Base.:(==)(a::BoxFace, b::BoxFace) = isequal(a.iv, b.iv) && a.side == b.side
Base.hash(a::BoxFace, h::UInt) = hash(a.side, hash(a.iv, h))
Base.show(io::IO, f::BoxFace) = print(io, "BoxFace(", f.iv, ", :", f.side, ")")

boundary_face(b::LowerBoundary) = BoxFace(b.x, :lo)
boundary_face(b::UpperBoundary) = BoxFace(b.x, :hi)

"Outward unit normal of `face` in the coordinate order `sivs`: `-e_iv` on `:lo`, `+e_iv` on `:hi`."
face_normal(face::BoxFace, sivs) = [isequal(x, face.iv) ? (face.side === :lo ? -1 : 1) : 0 for x in sivs]

_depends(ex, syms) = any(s -> any(isequal(unwrap(s)), Symbolics.get_variables(ex)), syms)
_iszero(ex) = isequal(Symbolics.simplify(Num(ex)), 0)

"Is the call `ex` of field `f` the trace of that field on `face`?"
function _on_face(ctx::SlotContext, ex, f, face::BoxFace)
    declared = ctx.args[ctx.ops[f]]
    args = arguments(ex)
    length(args) == length(declared) || return false
    lo, hi = ctx.intervals[face.iv]
    endpoint = face.side === :lo ? lo : hi
    return all(zip(declared, args)) do (iv, arg)
        isequal(iv, face.iv) || return isequal(iv, arg)
        val = unwrap_const(arg)
        return val isa Number && abs(val - endpoint) <= 8 * eps(Float64) * abs(hi - lo)
    end
end

# To be replaced by PDEBase.substitute_traces (spec §4.2) once PDEBase.boundary lands: plain
# `substitute` does not descend into `Differential`, so `Dx(u(t, 0))` would keep its trace
# (report 6 S1). Traces not on `face` are left in place for `_assert_no_trace` to report.
function _slot_substitute(ctx::SlotContext, ex, face, slots)
    ex = unwrap(ex)
    iscall(ex) || return _slot_leaf(ctx, ex, slots)
    j = findfirst(p -> isequal(p, ex), ctx.params)
    j === nothing || return unwrap(slots.P[j])
    op = operation(ex)
    op isa SymbolicUtils.Operator && return _slot_operator(ctx, op, ex, face, slots)
    f = _field(ctx, ex)
    if f !== nothing
        return _on_face(ctx, ex, f, face) ? unwrap(slots.U[f]) : ex
    end
    args = map(a -> _slot_substitute(ctx, a, face, slots), arguments(ex))
    return maketerm(typeof(ex), op, args, metadata(ex))
end

# A `NormalDerivative(Γ)` method goes here when Symbolics gains the operator: its trace of
# field `f` becomes `N · G_f`.
function _slot_operator(ctx::SlotContext, op::Differential, ex, face, slots)
    inner = unwrap(only(arguments(ex)))
    f = _field(ctx, inner)
    if f === nothing
        _has_trace(ctx, inner) && error(
            iscall(inner) && operation(inner) isa Differential ?
                "nested derivative in a boundary condition is not supported: $(wrap(ex))" :
                "derivative of an expression containing a trace is not supported in a boundary condition; differentiate the trace itself: $(wrap(ex))"
        )
        # Trace-free data such as `Dy(sin(y))` is safe to expand; only traces meet the §3.3 trap.
        expanded = unwrap(Symbolics.expand_derivatives(wrap(ex)))
        _has_differential(expanded) &&
            error("unresolved symbols in boundary condition: $(wrap(expanded)) has no derivative rule; register the derivative")
        return _slot_substitute(ctx, expanded, face, slots)
    end
    _on_face(ctx, inner, f, face) || return ex
    d = findfirst(x -> isequal(x, op.x), ctx.sivs)
    d === nothing &&
        error("time derivative of a trace is not supported in a boundary condition: $(wrap(ex))")
    op.order == 1 ||
        error("derivative of order $(op.order) in a boundary condition is not supported: $(wrap(ex))")
    return unwrap(slots.G[(f - 1) * length(ctx.sivs) + d])
end
function _slot_operator(ctx::SlotContext, op, ex, face, slots)
    return error("operator $op is not supported in a boundary condition: $(wrap(ex))")
end

_has_differential(ex) =
    iscall(ex) && (operation(ex) isa Differential || any(a -> _has_differential(unwrap(a)), arguments(ex)))

_has_trace(ctx::SlotContext, ex) =
    iscall(ex) && (_field(ctx, ex) !== nothing || any(a -> _has_trace(ctx, unwrap(a)), arguments(ex)))

function _assert_no_trace(ctx::SlotContext, ex, eq, face)
    ex = unwrap(ex)
    iscall(ex) || return nothing
    _field(ctx, ex) === nothing ||
        error("boundary condition $eq on $face still contains the trace $(wrap(ex)) after slot substitution; every trace must lie on $face")
    foreach(a -> _assert_no_trace(ctx, a, eq, face), arguments(ex))
    return nothing
end

"""
    slot_condition(eq, face, ctx::SlotContext, slots) -> Num

The residual `eq.lhs - eq.rhs` of a boundary condition on `face` in slot form: each trace
`u_f` becomes `U[f]`, `Dx_d(u_f)` becomes `G[(f-1)*dim+d]` and the face coordinate becomes its
endpoint value. Errors if a trace or any other unslotted symbol remains. Only trace-free
subexpressions are derivative-expanded (spec §3.3 trap).
"""
function slot_condition(eq::Equation, face, ctx::SlotContext, slots)
    B = _slot_substitute(ctx, unwrap(eq.lhs - eq.rhs), face, slots)
    _assert_no_trace(ctx, B, eq, face)
    slotset = Set(unwrap.(vcat(slots.U, slots.G, slots.X, [slots.T], slots.N, slots.P)))
    leftover = unique(_leftovers!(Any[], B, slotset, false))
    isempty(leftover) ||
        error("unresolved symbols $(Num.(leftover)) in boundary condition $eq on $face; declare them as parameters or register the derivative")
    return restrict_to_face(Num(B), face, ctx, slots)
end

"""
    restrict_to_face(ex, face::BoxFace, ctx::SlotContext, slots)

`ex` with the face's own coordinate slot `X[d]` set to the face's endpoint, so coefficients are
compared on the face rather than in the volume.
"""
function restrict_to_face(ex, face::BoxFace, ctx::SlotContext, slots)
    d = findfirst(x -> isequal(x, face.iv), ctx.sivs)
    lo, hi = ctx.intervals[face.iv]
    return substitute(ex, Dict(slots.X[d] => (face.side === :lo ? lo : hi)))
end

# Proportionality at random slot values (open question O3). Sampling cannot certify identity, so
# this never accepts a match; it only sharpens the rejection message. Values are signed and
# log-uniform in [0.1, 10]; a sample that throws (`sqrt` of a negative value) is redrawn.
function _numerically_proportional(aB, aΦ; rtol = 1.0e-8, nsamples = 8, maxattempts = 40)
    syms = union(Symbolics.get_variables.(vcat(aB, aΦ))...)
    ev(ex, vals) = Float64(Symbolics.value(substitute(Num(ex), vals; fold = Val(true))))
    accepted = 0
    for attempt in 1:maxattempts
        signof() = attempt == 1 ? -1 : attempt == 2 ? 1 : rand((-1, 1))
        vals = Dict(s => signof() * exp10(2 * rand() - 1) for s in syms)
        b, a = try
            ev.(aB, (vals,)), ev.(aΦ, (vals,))
        catch
            continue
        end
        aa = sum(abs2, a)
        λ = aa == 0 ? 0.0 : sum(b .* a) / aa
        sqrt(sum(abs2, b .- λ .* a)) <= rtol * (sqrt(sum(abs2, b)) + abs(λ) * sqrt(aa)) ||
            return false
        accepted += 1
        accepted == nsamples && return true
    end
    return accepted > 0
end

"""
    match_flux(Φ, B, slots) -> Union{Num, Nothing}

Match the boundary residual `B` against the conormal flux `Φ = G·N` (spec §5.3). Both must be
affine in the gradient slots, `Φ = a_Φ·g + Φ₀` and `B = a_B·g + B₀`; if `a_B = λ a_Φ` for a scalar
`λ` the matched flux `Φ₀ - B₀/λ` is returned, otherwise `nothing`. Proportionality must be proven
by `simplify`; a numeric agreement is never accepted (open question O3).
"""
function match_flux(Φ, B, slots)
    g = slots.G
    aΦ = [Symbolics.derivative(Φ, gi) for gi in g]
    aB = [Symbolics.derivative(B, gi) for gi in g]
    any(a -> _depends(a, g), vcat(aΦ, aB)) && return nothing
    i = findfirst(!_iszero, aΦ)
    (i === nothing || all(_iszero, aB)) && return nothing
    λ = Symbolics.simplify(Symbolics.simplify_fractions(aB[i] / aΦ[i]))
    _iszero(λ) && return nothing
    # Expanded and fraction-simplified so polynomial identities such as (1 + U)^2 = 1 + 2U + U^2
    # are proven; no coefficient-by-coefficient products, which can underflow to zero.
    all(j -> _iszero(Symbolics.simplify_fractions(Symbolics.expand(aB[j] - λ * aΦ[j]))), eachindex(g)) ||
        return nothing
    nograd = Dict(gi => 0 for gi in g)
    return Symbolics.simplify(substitute(Φ, nograd) - substitute(B, nograd) / λ)
end

# Is `B` a normal derivative of one field, i.e. its gradient coefficients ∝ the normal?
function _is_normal_derivative(B, n, slots)
    dim = length(n)
    aB = [Symbolics.derivative(B, gi) for gi in slots.G]
    blocks = [aB[((f - 1) * dim + 1):(f * dim)] for f in eachindex(slots.U)]
    nonzero = filter(blk -> !all(_iszero, blk), blocks)
    return length(nonzero) == 1 && _numerically_proportional(only(nonzero), Num.(n))
end

_gradcoefs(ex, slots) = [Symbolics.derivative(ex, gi) for gi in slots.G]
_conormal_flux(G, n) = sum(G[d] * n[d] for d in eachindex(n))

"""
    route_natural(B, flux, n, slots; label = B) -> (eq::Int, Φ::Num)

Route the natural condition `B` with outward normal `n` to the one equation whose conormal flux
`flux[e]·n` it matches (spec §5.3), returning that equation and the matched flux. Errors if no
equation or more than one matches.
"""
function route_natural(B, flux, n, slots; label = B)
    nsub = Dict(slots.N[d] => n[d] for d in eachindex(n))
    Bn = substitute(B, nsub)
    matches = Tuple{Int, Num}[]
    for (e, G) in enumerate(flux)
        Φbc = match_flux(substitute(_conormal_flux(G, n), nsub), Bn, slots)
        Φbc === nothing || push!(matches, (e, Φbc))
    end
    length(matches) == 1 && return only(matches)
    if length(matches) > 1
        error("boundary condition $label is a conormal condition for the fluxes of equations $(matches[1][1]) and $(matches[2][1]); the routing is ambiguous")
    end
    hint = _is_normal_derivative(Bn, n, slots) ? "write the conormal condition" :
        "write the normal derivative"
    aB = _gradcoefs(Bn, slots)
    near = findfirst(flux) do G
        aΦ = _gradcoefs(substitute(_conormal_flux(G, n), nsub), slots)
        !any(a -> _depends(a, slots.G), vcat(aB, aΦ)) && !all(_iszero, aΦ) &&
            _numerically_proportional(aB, aΦ)
    end
    if near !== nothing
        hint *= ". Its gradient coefficients agree numerically with the conormal of equation $near at sample points, but proportionality could not be proven symbolically; write the condition as a multiple of that conormal flux"
    end
    return error("boundary condition $label is not a conormal condition for any equation's flux; $hint")
end

"""
    essential_condition(B, slots; field = nothing, label = B) -> (field::Int, g::Num)

Solve the gradient-free boundary residual `B = 0` for `U[field] = g(X, T, P)`. `B` must be linear
in exactly one field and `g` free of the others; `field`, if given, is the field to solve for.
"""
function essential_condition(B, slots; field = nothing, label = B)
    U = slots.U
    present = [f for f in eachindex(U) if _depends(B, [U[f]])]
    isempty(present) && error("Dirichlet data must be linear in one field: $label contains no field")
    f = field === nothing ? first(present) : field
    others = [U[j] for j in eachindex(U) if j != f]
    (field === nothing && length(present) > 1 || f ∉ present) &&
        error("Dirichlet data depends on another field: $label")
    a = Symbolics.derivative(B, U[f])
    _depends(a, [U[f]]) && error("Dirichlet data must be linear in one field: $label")
    g = Symbolics.simplify(-substitute(B, Dict(U[f] => 0)) / a)
    (_depends(a, others) || _depends(g, others)) &&
        error("Dirichlet data depends on another field: $label")
    return f, g
end

"""
    periodic_condition(b::InterfaceBoundary, ctx::SlotContext) -> PeriodicBC{BoxFace}

A condition `u(.., x_lo, ..) ~ u(.., x_hi, ..)` equating one field on the two faces of one axis,
as `parse_bcs` records it; the shift is the axis length.
"""
function periodic_condition(b::InterfaceBoundary, ctx::SlotContext)
    f = _field(ctx, unwrap(b.u))
    (f === nothing || !isequal(operation(unwrap(b.u)), operation(unwrap(b.u2)))) &&
        error("periodic condition $(b.eq) must equate the same field on two faces")
    isequal(b.x, b.x2) ||
        error("periodic pairing needs a translation; general maps are unsupported: $(b.eq)")
    residual = Num(unwrap(b.eq.lhs - b.eq.rhs))
    pair = Num(unwrap(b.u)) - Num(unwrap(b.u2))
    (_iszero(residual - pair) || _iszero(residual + pair)) ||
        error("periodic condition $(b.eq) must have the form u(.., x_lo, ..) ~ u(.., x_hi, ..)")
    lo, hi = ctx.intervals[b.x]
    return PeriodicBC(BoxFace(b.x, :lo), BoxFace(b.x, :hi), f, hi - lo)
end

# Routing on a face where some equation's conormal vanishes identically (e.g. the weighted
# r²∂ᵣ flux at r = 0): a condition that matches only that equation's volume conormal imposes a
# zero facet flux if it is homogeneous, and cannot be imposed otherwise.
function _route_on_face(B, flux, faceflux, face, n, slots, label)
    vanishing = [all(_iszero, _gradcoefs(_conormal_flux(G, n), slots)) for G in faceflux]
    if any(vanishing) &&
            all(G -> match_flux(_conormal_flux(G, n), B, slots) === nothing, faceflux)
        degenerate = filter(findall(vanishing)) do e
            match_flux(_conormal_flux(flux[e], n), B, slots) !== nothing
        end
        if length(degenerate) == 1
            e = only(degenerate)
            _iszero(substitute(B, Dict(g => 0 for g in slots.G))) ||
                error("the conormal flux of equation $e vanishes identically on $face; the condition $(label) cannot be imposed there")
            return e, Num(0)
        end
    end
    return route_natural(B, faceflux, n, slots; label)
end

_gradient_fields(B, slots, dim) =
    [f for f in eachindex(slots.U) if _depends(B, slots.G[((f - 1) * dim + 1):(f * dim)])]

"""
    boundary_conditions(boundarymap, ctx::SlotContext, slots, flux)
        -> (essential, natural, periodic)

Classify every spatial boundary condition of a `parse_bcs` boundary map. A gradient-free
condition is essential for the field `parse_bcs` filed it under; any other is routed to the
equation whose conormal flux it matches. `flux[e]` is equation `e`'s flux vector in slots.
Errors if two conditions fix the same field, or give the same equation's natural condition, on
one face, or if one face has both an essential and a natural condition on the same field.
"""
function boundary_conditions(boundarymap, ctx::SlotContext, slots, flux)
    essential = EssentialBC{BoxFace}[]
    natural = NaturalBC{BoxFace}[]
    periodic = PeriodicBC{BoxFace}[]
    dim = length(ctx.sivs)
    fixed = Dict{Tuple{BoxFace, Int}, String}()        # (face, field) => essential condition
    fluxed = Dict{Tuple{BoxFace, Int}, String}()       # (face, eq) => natural condition
    graded = Dict{Tuple{BoxFace, Int}, String}()       # (face, field) => natural condition
    for (op, bykey) in boundarymap, (x, bs) in bykey
        ctx.t !== nothing && isequal(x, ctx.t) && continue
        for b in bs
            if b isa InterfaceBoundary
                isupper(b) || push!(periodic, periodic_condition(b, ctx))
                continue
            end
            b isa Union{LowerBoundary, UpperBoundary} ||
                error("boundary condition $(b.eq) relates derivatives on two faces; only u(.., x_lo, ..) ~ u(.., x_hi, ..) periodic pairs are supported")
            face = boundary_face(b)
            B = slot_condition(b.eq, face, ctx, slots)
            label = "$(b.eq) on $face"
            if _depends(B, slots.G)
                faceflux = [restrict_to_face.(G, (face,), (ctx,), (slots,)) for G in flux]
                e, Φ = _route_on_face(B, flux, faceflux, face, face_normal(face, ctx.sivs), slots, label)
                haskey(fluxed, (face, e)) &&
                    error("boundary conditions $(fluxed[(face, e)]) and $label both give the natural condition of equation $e")
                fluxed[(face, e)] = label
                for f in _gradient_fields(B, slots, dim)
                    haskey(fixed, (face, f)) &&
                        error("boundary conditions $(fixed[(face, f)]) and $label impose both an essential and a natural condition on field $f")
                    graded[(face, f)] = label
                end
                push!(natural, NaturalBC(face, e, Φ))
            else
                f = findfirst(o -> isequal(o, op), ctx.ops)
                f, g = essential_condition(B, slots; field = f, label)
                haskey(fixed, (face, f)) &&
                    error("boundary conditions $(fixed[(face, f)]) and $label both fix field $f")
                haskey(graded, (face, f)) &&
                    error("boundary conditions $(graded[(face, f)]) and $label impose both an essential and a natural condition on field $f")
                fixed[(face, f)] = label
                push!(essential, EssentialBC(face, f, g))
            end
        end
    end
    return essential, natural, periodic
end
