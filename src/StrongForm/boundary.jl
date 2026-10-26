# Stage d (spec §5.3): put each boundary condition into slots, then classify it as essential,
# natural (matched against an equation's conormal flux) or periodic.

# A boundary key is the classifier's `BoxFace` for an endpoint face and the `BoundaryPiece`
# for a trace on a piece. `_location` is where the key lies by the classifier, so a piece that
# is a face of the box is that `BoxFace`, with its numeric normal and endpoint.
boundary_face(b::LowerBoundary) = BoxFace(b.x, :lo)
boundary_face(b::UpperBoundary) = BoxFace(b.x, :hi)
boundary_face(b::PieceBoundary) = b.piece

_location(ctx::SlotContext, face::BoxFace) = face
_location(ctx::SlotContext, piece::BoundaryPiece) = location(ctx.v, [Piece(piece)])

"Outward unit normal of `face` in the coordinate order `sivs`: `-e_iv` on `:lo`, `+e_iv` on `:hi`."
face_normal(face::BoxFace, sivs) = [isequal(x, face.iv) ? (face.side === :lo ? -1 : 1) : 0 for x in sivs]

# The normal of a location in slots: numeric on a box face; the symbolic slots `N` on a piece,
# whose normal the backend supplies per quadrature point (spec §3.4).
_normal(face::BoxFace, ctx::SlotContext, slots) = face_normal(face, ctx.sivs)
_normal(::Piece, ctx::SlotContext, slots) = slots.N

_depends(ex, syms) = any(s -> any(isequal(unwrap(s)), Symbolics.get_variables(ex)), syms)
_iszero(ex) = isequal(Symbolics.simplify(Num(ex)), 0)

# How messages name a location or key, the same however `BoxFace` is imported.
_where(loc::BoxFace) = "BoxFace($(loc.iv), :$(loc.side))"
_where(loc::Piece) = string(loc.piece)
_where(piece::BoundaryPiece) = string(piece)

# Is the call `ex` of field `f` a trace on `loc`? A call of a field that does not have the
# coordinates `loc` fixes, such as `w(t)` or `q(t, y)` on an x-face, is at its declared
# arguments everywhere on `loc` and counts as on it.
function _on_location(ctx::SlotContext, ex, f, loc)
    callloc = location(ctx.v, classify_call(ctx.v, ex)[2])
    isequal(callloc, loc) && return true
    callloc isa Interior || return false
    fixed = loc isa BoxFace ? Any[loc.iv] : ctx.sivs
    return !any(a -> any(isequal(a), fixed), ctx.args[ctx.ops[f]])
end

function _check_on_location(ctx::SlotContext, ex, f, loc, label)
    _on_location(ctx, ex, f, loc) ||
        error("boundary condition $label still contains the trace $(wrap(ex)) after slot substitution; every trace must lie on $(_where(loc))")
    return nothing
end

# The substitution rules that put the boundary expression `ex` on `loc` into slots, for
# `substitute_traces`: traces become `U`, their first derivatives `G`, their normal derivatives
# `N·G`, coordinates `X`, time `T` and parameters `P`. Only trace-free subexpressions are
# derivative-expanded (spec §3.3 trap).
function _slot_rules!(rules, ctx::SlotContext, ex, loc, slots, label)
    ex = unwrap(ex)
    haskey(rules, ex) && return rules
    if !iscall(ex)
        slot = _slot_leaf(ctx, ex, slots)
        slot === ex || (rules[ex] = slot)
        return rules
    end
    j = findfirst(p -> isequal(p, ex), ctx.params)
    if j !== nothing
        rules[ex] = unwrap(slots.P[j])
        return rules
    end
    op = operation(ex)
    if op isa SymbolicUtils.Operator
        rules[ex] = _slot_operator(ctx, op, ex, loc, slots, label)
        return rules
    end
    f = _field(ctx, ex)
    if f !== nothing
        _check_on_location(ctx, ex, f, loc, label)
        rules[ex] = unwrap(slots.U[f])
        return rules
    end
    foreach(a -> _slot_rules!(rules, ctx, a, loc, slots, label), arguments(ex))
    return rules
end

function _slot_expression(ctx::SlotContext, ex, loc, slots, label)
    rules = _slot_rules!(Dict{Any, Any}(), ctx, ex, loc, slots, label)
    return unwrap(assert_no_traces(substitute_traces(ex, rules)))
end

function _slot_operator(ctx::SlotContext, op::Differential, ex, loc, slots, label)
    inner = unwrap(only(arguments(ex)))
    f = _field(ctx, inner)
    if f === nothing
        _has_trace(ctx, inner) && error(
            iscall(inner) && operation(inner) isa Union{Differential, NormalDerivative} ?
                "nested derivative in a boundary condition is not supported: $(wrap(ex))" :
                "derivative of an expression containing a trace is not supported in a boundary condition; differentiate the trace itself: $(wrap(ex))"
        )
        # Trace-free data such as `Dy(sin(y))` is safe to expand; only traces meet the §3.3 trap.
        expanded = unwrap(Symbolics.expand_derivatives(wrap(ex)))
        _has_differential(expanded) &&
            error("unresolved symbols in boundary condition: $(wrap(expanded)) has no derivative rule; register the derivative")
        return _slot_expression(ctx, expanded, loc, slots, label)
    end
    _check_on_location(ctx, inner, f, loc, label)
    d = findfirst(x -> isequal(x, op.x), ctx.sivs)
    if d === nothing
        ctx.t !== nothing && isequal(op.x, ctx.t) &&
            error("time derivative of a trace is not supported in a boundary condition: $(wrap(ex))")
        error("derivative in $(op.x), which is not an independent variable of the system, in a boundary condition: $(wrap(ex))")
    end
    loc isa Piece && error(
        "$(wrap(ex)) is a Cartesian derivative of a trace on $(loc.piece), which equals the normal derivative only on a planar piece perpendicular to $(op.x); write the normal derivative NormalDerivative($(loc.piece))(…)"
    )
    op.order == 1 ||
        error("derivative of order $(op.order) in a boundary condition is not supported: $(wrap(ex))")
    return unwrap(slots.G[(f - 1) * length(ctx.sivs) + d])
end

# `NormalDerivative(Γ)(u_f(t, Γ))` is `N·G_f` (spec §5.3), with the numeric normal on a box face.
function _slot_operator(ctx::SlotContext, op::NormalDerivative, ex, loc, slots, label)
    inner = unwrap(only(arguments(ex)))
    f = _field(ctx, inner)
    f === nothing && error(
        _has_trace(ctx, inner) ?
            "normal derivative of an expression containing a trace is not supported in a boundary condition; differentiate the trace itself: $(wrap(ex))" :
            "normal derivative of $(wrap(inner)), which is not a trace, in a boundary condition: $(wrap(ex))"
    )
    op.piece isa BoundaryPiece ||
        error("the normal derivative $(wrap(ex)) is not taken on a boundary piece")
    oploc = location(ctx.v, [Piece(op.piece)])
    isequal(oploc, loc) ||
        error("$(wrap(ex)) takes the normal derivative on $(op.piece), but the condition $label is on $(_where(loc)); write the condition on the operator's piece")
    _check_on_location(ctx, inner, f, loc, label)
    dim = length(ctx.sivs)
    n = _normal(loc, ctx, slots)
    return unwrap(sum(n[d] * slots.G[(f - 1) * dim + d] for d in 1:dim))
end

function _slot_operator(ctx::SlotContext, op, ex, loc, slots, label)
    return error("operator $op is not supported in a boundary condition: $(wrap(ex))")
end

_has_differential(ex) =
    iscall(ex) && (operation(ex) isa Differential || any(a -> _has_differential(unwrap(a)), arguments(ex)))

_has_trace(ctx::SlotContext, ex) =
    iscall(ex) && (_field(ctx, ex) !== nothing || any(a -> _has_trace(ctx, unwrap(a)), arguments(ex)))

"""
    slot_condition(eq, face, ctx::SlotContext, slots) -> Num

The residual `eq.lhs - eq.rhs` of a boundary condition on `face`, a `PDEBase.BoxFace` or a
`BoundaryPiece`, in slot form: each trace `u_f` becomes `U[f]`, `Dx_d(u_f)` becomes
`G[(f-1)*dim+d]` and `NormalDerivative(Γ)(u_f)` becomes `Σ_d N[d] G[(f-1)*dim+d]`. On a box face
the normal is the face's `±e_iv` and the face coordinate becomes its endpoint value; on a piece
the normal stays the symbolic slots `N`, and a Cartesian derivative of a trace is an error.
Traces are replaced with `PDEBase.substitute_traces`. Errors if a trace or any other unslotted
symbol remains. Only trace-free subexpressions are derivative-expanded (spec §3.3 trap).
"""
slot_condition(eq::Equation, face, ctx::SlotContext, slots) =
    restrict_to_face(_slot_unrestricted(eq, face, ctx, slots), face, ctx, slots)

function _slot_unrestricted(eq::Equation, face, ctx::SlotContext, slots)
    loc = _location(ctx, face)
    B = _slot_expression(ctx, unwrap(eq.lhs - eq.rhs), loc, slots, "$eq on $(_where(loc))")
    slotset = Set(unwrap.(vcat(slots.U, slots.G, slots.X, [slots.T], slots.N, slots.P)))
    leftover = unique(_leftovers!(Any[], B, slotset, false))
    isempty(leftover) ||
        error("unresolved symbols $(Num.(leftover)) in boundary condition $eq on $(_where(loc)); declare them as parameters or register the derivative")
    return Num(B)
end

"""
    restrict_to_face(ex, face, ctx::SlotContext, slots)

`ex` with a box face's own coordinate slot `X[d]` set to the face's endpoint, so coefficients
are compared on the face rather than in the volume. The identity on a `BoundaryPiece`, whose
coordinates vary over the piece, unless the piece is a face of the box.
"""
function restrict_to_face(ex, face::BoxFace, ctx::SlotContext, slots)
    d = findfirst(x -> isequal(x, face.iv), ctx.sivs)
    lo, hi = ctx.intervals[face.iv]
    return substitute(ex, Dict(slots.X[d] => (face.side === :lo ? lo : hi)))
end
restrict_to_face(ex, ::Piece, ctx::SlotContext, slots) = ex
restrict_to_face(ex, piece::BoundaryPiece, ctx::SlotContext, slots) =
    restrict_to_face(ex, _location(ctx, piece), ctx, slots)

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

# A condition written with a coefficient that vanishes on the face, such as r²∂ᵣu ~ 0 at r = 0,
# reads 0 ~ 0 there. It is the zero natural condition of the one equation whose volume conormal
# it matches only if that equation's whole facet flux vanishes on the face too; otherwise it
# would silently impose a zero flux the user never wrote.
function _route_vanishing(Bvolume, flux, faceflux, n, slots, label)
    matches = filter(e -> match_flux(_conormal_flux(flux[e], n), Bvolume, slots) !== nothing, eachindex(flux))
    isempty(matches) &&
        error("boundary condition $label is not a conormal condition for any equation's flux; write the normal derivative")
    length(matches) == 1 ||
        error("boundary condition $label is a conormal condition for the fluxes of equations $(matches[1]) and $(matches[2]); the routing is ambiguous")
    e = only(matches)
    _iszero(_conormal_flux(faceflux[e], n)) ||
        error("the condition $label reduces to 0 ~ 0 on the face and imposes nothing there, while equation $e has a nonzero flux on that face; write the condition without the vanishing factor")
    return e
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
        error("boundary condition $(b.eq) relates values on two faces of different axes, $(b.x) and $(b.x2); only u(.., x_lo, ..) ~ u(.., x_hi, ..) periodic pairs are supported")
    residual = Num(unwrap(b.eq.lhs - b.eq.rhs))
    pair = Num(unwrap(b.u)) - Num(unwrap(b.u2))
    (_iszero(residual - pair) || _iszero(residual + pair)) ||
        error("periodic condition $(b.eq) must have the form u(.., x_lo, ..) ~ u(.., x_hi, ..)")
    lo, hi = ctx.intervals[b.x]
    return PeriodicBC(BoxFace(b.x, :lo), BoxFace(b.x, :hi), f, hi - lo)
end

# Routing on a face where some equation's conormal gradient coefficients vanish (e.g. the
# weighted r²∂ᵣ flux at r = 0). If its whole facet flux vanishes, a homogeneous condition that
# matches its volume conormal is the zero natural condition, and other data cannot be imposed.
# If a gradient-free facet flux remains, no natural condition can prescribe it.
function _route_on_face(B, flux, faceflux, face, n, slots, label)
    vanishing = [all(_iszero, _gradcoefs(_conormal_flux(G, n), slots)) for G in faceflux]
    if any(vanishing) &&
            all(G -> match_flux(_conormal_flux(G, n), B, slots) === nothing, faceflux)
        degenerate = filter(findall(vanishing)) do e
            match_flux(_conormal_flux(flux[e], n), B, slots) !== nothing
        end
        if length(degenerate) == 1
            e = only(degenerate)
            Φ₀ = _conormal_flux(faceflux[e], n)
            _iszero(Φ₀) ||
                error("the conormal gradient coefficients of equation $e vanish on $face; its facet flux there is $Φ₀ and cannot be prescribed by $label")
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
                conflict(poffset + k, i, "conflict: field $(p.field) is periodic on $(_where(face))")
        end
        for j in eachindex(nat)
            nat[j][1] == face && owner[nat[j][2]] == p.field &&
                conflict(poffset + k, offset + j, "conflict: field $(p.field) is periodic on $(_where(face))")
        end
    end
    return nothing
end

"""
    boundary_conditions(boundarymap, ctx::SlotContext, slots, flux; mass = nothing, bcs = nothing)
        -> (essential, natural, periodic)

Classify every spatial boundary condition of a `parse_bcs` boundary map: the endpoint faces of
a box discretization's map, keyed by their `PDEBase.BoxFace`, and the `PieceBoundary` entries
of a mesh discretization's map, keyed by their `BoundaryPiece`; a mixed map gives a union key
type. Each output vector lists its conditions in the order of `bcs`, the system's boundary
conditions; without `bcs` the order is that of the boundary map. Which condition wins at a
corner is the discretizer's choice. A gradient-free condition is essential for the field
`parse_bcs` filed it under; any other is routed to the equation whose conormal flux it matches,
with the normal `±e_iv` on a box face (or on a piece that is a face of the box) and the
symbolic normal slots `N` on any other piece. `flux[e]` is equation `e`'s flux vector and
`mass[e]` its time-derivative part, in slots.

Errors if two conditions fix the same field, or give the same equation's natural condition, on
one face. With `mass`, an equation's own field is the one whose time derivative alone it
contains, or the one left over when a single equation and a single field remain unpaired (a
steady single-field problem); where the pairing stays ambiguous these checks are skipped.
Otherwise it is also an error to give a face both an essential condition on a field and a
natural condition for that field's equation, or to give either kind of condition on a face of a
field's periodic pair.
"""
function boundary_conditions(
        boundarymap, ctx::SlotContext, slots, flux; mass = nothing, bcs = nothing
    )
    rank(eq) = bcs === nothing ? 0 : something(findfirst(c -> isequal(c, eq), bcs), length(bcs) + 1)
    located = Tuple{Int, Equation, Any, Union{Int, Nothing}}[]    # (rank, condition, key, field)
    found = Tuple{Int, Any, String}[]    # (rank, classified condition, label)
    for (op, bykey) in boundarymap, (x, bs) in bykey
        ctx.t !== nothing && isequal(x, ctx.t) && continue
        for b in bs
            if b isa InterfaceBoundary
                isupper(b) || push!(found, (rank(b.eq), periodic_condition(b, ctx), string(b.eq)))
                continue
            end
            b isa Union{LowerBoundary, UpperBoundary, PieceBoundary} ||
                error("boundary condition $(b.eq) relates derivatives on two faces; only u(.., x_lo, ..) ~ u(.., x_hi, ..) periodic pairs are supported")
            push!(located, (rank(b.eq), b.eq, boundary_face(b), findfirst(o -> isequal(o, op), ctx.ops)))
        end
    end
    K = _key_type([key for (_, _, key, _) in located], !isempty(found))
    for (r, eq, key, field) in located
        loc = _location(ctx, key)
        n = _normal(loc, ctx, slots)
        Bvolume = _slot_unrestricted(eq, key, ctx, slots)
        B = restrict_to_face(Bvolume, key, ctx, slots)
        label = "$eq on $(_where(loc))"
        faceflux = [restrict_to_face.(G, (key,), (ctx,), (slots,)) for G in flux]
        condition = if _depends(Bvolume, slots.G) && _iszero(B)
            e = _route_vanishing(Bvolume, flux, faceflux, n, slots, label)
            NaturalBC{K}(key, e, Num(0))
        elseif !_depends(B, vcat(slots.G, slots.U)) && !_iszero(B)
            error("the condition $label reduces to a nonzero constant on the face and cannot be imposed")
        elseif _depends(B, slots.G)
            e, Φ = _route_on_face(B, flux, faceflux, _where(loc), n, slots, label)
            NaturalBC{K}(key, e, Φ)
        else
            f, g = essential_condition(B, slots; field, label)
            EssentialBC{K}(key, f, g)
        end
        push!(found, (r, condition, label))
    end
    sort!(found; by = first)
    essential = EssentialBC{K}[c for (_, c, _) in found if c isa EssentialBC]
    natural = NaturalBC{K}[c for (_, c, _) in found if c isa NaturalBC]
    periodic = PeriodicBC{K}[]
    labels = [[l for (_, c, l) in found if c isa T] for T in (EssentialBC, NaturalBC)]
    perlabels = String[]
    for (_, p, l) in found
        p isa PeriodicBC || continue
        any(q -> (q.master, q.mirror, q.field) == (p.master, p.mirror, p.field), periodic) && continue
        push!(periodic, p)
        push!(perlabels, l)
    end
    owner = mass === nothing ? fill(nothing, length(flux)) : _own_fields(mass, slots)
    _check_conflicts(essential, natural, periodic, vcat(labels..., perlabels), owner)
    foreach(b -> check_evaluable(b.g), essential)
    foreach(b -> check_evaluable(b.Φ), natural)
    return essential, natural, periodic
end

# The key type of a set of conditions: the one key type, or the union of a mixed set.
function _key_type(keys, hasperiodic)
    types = unique(vcat(Any[typeof(key) for key in keys], hasperiodic ? Any[BoxFace] : Any[]))
    return isempty(types) ? BoxFace : Union{types...}
end
