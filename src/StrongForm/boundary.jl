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
    if d === nothing
        ctx.t !== nothing && isequal(op.x, ctx.t) &&
            error("time derivative of a trace is not supported in a boundary condition: $(wrap(ex))")
        error("derivative in $(op.x), which is not an independent variable of the system, in a boundary condition: $(wrap(ex))")
    end
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
slot_condition(eq::Equation, face, ctx::SlotContext, slots) =
    restrict_to_face(_slot_unrestricted(eq, face, ctx, slots), face, ctx, slots)

function _slot_unrestricted(eq::Equation, face, ctx::SlotContext, slots)
    B = _slot_substitute(ctx, unwrap(eq.lhs - eq.rhs), face, slots)
    _assert_no_trace(ctx, B, eq, face)
    slotset = Set(unwrap.(vcat(slots.U, slots.G, slots.X, [slots.T], slots.N, slots.P)))
    leftover = unique(_leftovers!(Any[], B, slotset, false))
    isempty(leftover) ||
        error("unresolved symbols $(Num.(leftover)) in boundary condition $eq on $face; declare them as parameters or register the derivative")
    return Num(B)
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

# A condition written with a coefficient that vanishes on the face, such as r²∂ᵣu ~ 0 at r = 0:
# it is homogeneous there, so it is the zero natural condition of the one equation whose volume
# conormal it matches.
function _route_vanishing(Bvolume, flux, n, slots, label)
    matches = filter(e -> match_flux(_conormal_flux(flux[e], n), Bvolume, slots) !== nothing, eachindex(flux))
    length(matches) == 1 && return only(matches)
    isempty(matches) &&
        error("boundary condition $label is not a conormal condition for any equation's flux; write the normal derivative")
    return error("boundary condition $label is a conormal condition for the fluxes of equations $(matches[1]) and $(matches[2]); the routing is ambiguous")
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

# Each equation's own field: the unique field whose time derivative its mass contains, then by
# forced completion while exactly one equation and one field are left (a steady single-field
# problem, or the last algebraic row). `nothing` where the pairing stays ambiguous.
function _own_fields(mass, slots)
    owner = Union{Int, Nothing}[]
    for m in mass
        owned = [f for f in eachindex(slots.Ud) if _depends(m, [slots.Ud[f]])]
        push!(owner, length(owned) == 1 ? only(owned) : nothing)
    end
    while true
        open_eqs = findall(isnothing, owner)
        free = setdiff(eachindex(slots.U), filter(!isnothing, owner))
        (length(open_eqs) == 1 && length(free) == 1) || break
        owner[only(open_eqs)] = only(free)
    end
    return owner
end

# Conflicts between the classified conditions, each given with its label. An equation's own
# field is decided by `_own_fields`; an equation without one is not checked against essential
# or periodic conditions, since the discretizer owns its pairing.
function _check_conflicts(essential, natural, periodic, labels, owner)
    conflict(i, j, what) = error("boundary conditions $(labels[i]) and $(labels[j]) $what")
    ess = [(b.piece, b.field) for b in essential]
    nat = [(b.piece, b.eq) for b in natural]
    offset = length(essential)
    for i in eachindex(ess), j in (i + 1):length(ess)
        ess[i] == ess[j] && conflict(i, j, "both fix field $(ess[i][2])")
    end
    for i in eachindex(nat), j in (i + 1):length(nat)
        nat[i] == nat[j] &&
            conflict(offset + i, offset + j, "both give the natural condition of equation $(nat[i][2])")
    end
    for i in eachindex(ess), j in eachindex(nat)
        ess[i][1] == nat[j][1] && owner[nat[j][2]] == ess[i][2] &&
            conflict(i, offset + j, "impose both an essential and a natural condition on field $(ess[i][2])")
    end
    poffset = offset + length(natural)
    for (k, p) in enumerate(periodic), face in (p.master, p.mirror)
        for i in eachindex(ess)
            ess[i] == (face, p.field) &&
                conflict(poffset + k, i, "conflict: field $(p.field) is periodic on $face")
        end
        for j in eachindex(nat)
            nat[j][1] == face && owner[nat[j][2]] == p.field &&
                conflict(poffset + k, offset + j, "conflict: field $(p.field) is periodic on $face")
        end
    end
    return nothing
end

"""
    boundary_conditions(boundarymap, ctx::SlotContext, slots, flux; mass = nothing)
        -> (essential, natural, periodic)

Classify every spatial boundary condition of a `parse_bcs` boundary map. A gradient-free
condition is essential for the field `parse_bcs` filed it under; any other is routed to the
equation whose conormal flux it matches. `flux[e]` is equation `e`'s flux vector and `mass[e]`
its time-derivative part, in slots.

Errors if two conditions fix the same field, or give the same equation's natural condition, on
one face. With `mass`, an equation's own field is the one whose time derivative alone it
contains, or the one left over when a single equation and a single field remain unpaired (a
steady single-field problem); where the pairing stays ambiguous these checks are skipped.
Otherwise it is also an error to give a face both an essential condition on a field and a
natural condition for that field's equation, or to give either kind of condition on a face of a
field's periodic pair.
"""
function boundary_conditions(boundarymap, ctx::SlotContext, slots, flux; mass = nothing)
    essential = EssentialBC{BoxFace}[]
    natural = NaturalBC{BoxFace}[]
    periodic = PeriodicBC{BoxFace}[]
    esslabels, natlabels, perlabels = String[], String[], String[]
    for (op, bykey) in boundarymap, (x, bs) in bykey
        ctx.t !== nothing && isequal(x, ctx.t) && continue
        for b in bs
            if b isa InterfaceBoundary
                isupper(b) && continue
                p = periodic_condition(b, ctx)
                any(q -> (q.master, q.mirror, q.field) == (p.master, p.mirror, p.field), periodic) &&
                    continue
                push!(periodic, p)
                push!(perlabels, string(b.eq))
                continue
            end
            b isa Union{LowerBoundary, UpperBoundary} ||
                error("boundary condition $(b.eq) relates derivatives on two faces; only u(.., x_lo, ..) ~ u(.., x_hi, ..) periodic pairs are supported")
            face = boundary_face(b)
            Bvolume = _slot_unrestricted(b.eq, face, ctx, slots)
            B = restrict_to_face(Bvolume, face, ctx, slots)
            label = "$(b.eq) on $face"
            if _depends(Bvolume, slots.G) && _iszero(B)
                e = _route_vanishing(Bvolume, flux, face_normal(face, ctx.sivs), slots, label)
                push!(natural, NaturalBC(face, e, Num(0)))
                push!(natlabels, label)
            elseif !_depends(B, vcat(slots.G, slots.U)) && !_iszero(B)
                error("the condition $label reduces to a nonzero constant on the face and cannot be imposed")
            elseif _depends(B, slots.G)
                faceflux = [restrict_to_face.(G, (face,), (ctx,), (slots,)) for G in flux]
                e, Φ = _route_on_face(B, flux, faceflux, face, face_normal(face, ctx.sivs), slots, label)
                push!(natural, NaturalBC(face, e, Φ))
                push!(natlabels, label)
            else
                f = findfirst(o -> isequal(o, op), ctx.ops)
                f, g = essential_condition(B, slots; field = f, label)
                push!(essential, EssentialBC(face, f, g))
                push!(esslabels, label)
            end
        end
    end
    owner = mass === nothing ? fill(nothing, length(flux)) : _own_fields(mass, slots)
    _check_conflicts(essential, natural, periodic, vcat(esslabels, natlabels, perlabels), owner)
    foreach(b -> check_evaluable(b.g), essential)
    foreach(b -> check_evaluable(b.Φ), natural)
    return essential, natural, periodic
end
