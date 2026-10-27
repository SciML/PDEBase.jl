### ARGUMENT CLASSIFICATION ###

# One classification step for every call of a dependent variable: `u(t, x, y)`,
# `u(t, 0, y)` and the trace `u(t, Γ)` go through the same code, dispatching on the type of
# each argument.

"""
    ArgumentClass

Supertype of the classes [`classify_argument`](@ref) assigns to one argument of a call of
a dependent variable: [`Free`](@ref PDEBase.Free), [`BoxFace`](@ref PDEBase.BoxFace) or
[`Piece`](@ref PDEBase.Piece).
"""
abstract type ArgumentClass end

"""
    Free(iv) <: ArgumentClass

The argument is the declared independent variable `iv` in its own position, as `x` in
`u(t, x, y)`.
"""
struct Free <: ArgumentClass
    iv::Any
end

"""
    BoxFace(iv, side::Symbol) <: ArgumentClass

The argument is an endpoint of the interval of `iv`, the lower one when `side` is `:lo`
and the upper one when it is `:hi`, as `0` in `u(t, 0, y)` on `x ∈ Interval(0, 1)`. As a
location, it is the face of the box on which `iv` is fixed at that endpoint. A
[`Piece`](@ref PDEBase.Piece) that is a face of the box normalises to the same `BoxFace`
(see [`location`](@ref PDEBase.location)).
"""
struct BoxFace <: ArgumentClass
    iv::Any
    side::Symbol
end

"""
    Piece(piece::BoundaryPiece) <: ArgumentClass

The argument is a [`BoundaryPiece`](@ref) standing for all the spatial arguments, as `Γ`
in the trace `u(t, Γ)`.
"""
struct Piece <: ArgumentClass
    piece::BoundaryPiece
end

"""
    Interior()

The [`location`](@ref PDEBase.location) of a call whose arguments are all
[`Free`](@ref PDEBase.Free): the call is in the interior of the domain.
"""
struct Interior end

Base.isequal(a::Free, b::Free) = isequal(a.iv, b.iv)
Base.isequal(a::BoxFace, b::BoxFace) = isequal(a.iv, b.iv) && a.side === b.side
Base.isequal(a::Piece, b::Piece) = isequal(a.piece, b.piece)
Base.:(==)(a::T, b::T) where {T <: ArgumentClass} = isequal(a, b)
Base.hash(a::Free, h::UInt) = hash(a.iv, hash(Free, h))
Base.hash(a::BoxFace, h::UInt) = hash(a.side, hash(a.iv, hash(BoxFace, h)))
Base.hash(a::Piece, h::UInt) = hash(a.piece, hash(Piece, h))

_declared_call(v, u_op) = Num(u_op(v.args[u_op]...))

# The interval of `iv` converted to its float type, for the endpoint test.
function _float_interval(v, iv)
    lo, hi = v.intervals[iv]
    T = float(promote_type(typeof(lo), typeof(hi)))
    return T, convert(T, lo), convert(T, hi)
end

# Whether `val` is the endpoint `e` of an interval of width `width`, to within `8ϵ` relative
# to `max(width, |e|)`. Infinite endpoints match only exactly, and an infinite width falls
# back to the magnitude of the finite endpoint.
function _on_endpoint(val, e, width, ϵ)
    isfinite(e) || return val == e
    isfinite(val) || return false
    return abs(val - e) <= 8 * ϵ * (isfinite(width) ? max(width, abs(e)) : abs(e))
end

"""
    classify_argument(v::VariableMap, u_op, position::Int, arg) -> ArgumentClass

Classify argument `arg` of a call of the dependent variable `u_op`, given in the place of
the `position`-th argument of its declared argument list `v.args[u_op]`.

  - A [`BoundaryPiece`](@ref) whose parent is the spatial domain (`isequal` to
    [`spatial_domain(v)`](@ref spatial_domain)) is a [`Piece`](@ref PDEBase.Piece).
  - A number or numeric constant is converted to the float type of the interval of the
    declared variable `iv`. It is a [`BoxFace`](@ref PDEBase.BoxFace) when it lies within
    `8eps·max(width, |endpoint|)` of an endpoint, where `eps` is that of the coarser of the
    interval's and the number's float types, so `u(t, π)` and `u(t, Float32(π))` on
    `Interval(0.0, π)` and `u(t, 0.1 + 0.2)` on `Interval(0.0, 0.3)` are on the upper face.
    An infinite endpoint matches only an infinite number of the same sign. An exact endpoint
    always wins; when a number is within the tolerance of both endpoints of a narrow
    interval, the nearer endpoint wins, and a number equally close to both is an error. Any
    other number is an error naming the nearest endpoint and the distance to it.
  - The declared variable `iv` itself is [`Free`](@ref PDEBase.Free).
  - Anything else is an error that shows the declared argument list.
"""
function classify_argument(v::VariableMap, u_op, position::Int, arg)
    arg = unwrap_const(safe_unwrap(arg))
    declared = v.args[u_op]
    iv = declared[position]
    if arg isa BoundaryPiece
        Ω = spatial_domain(v)
        isequal(parent(arg), Ω) || throw(
            ArgumentError(
                "$(arg) in $(_declared_call(v, u_op)) is a piece of the boundary of $(parent(arg)), " *
                    "not of the spatial domain $(Ω)"
            )
        )
        return Piece(arg)
    elseif arg isa Number
        haskey(v.intervals, iv) || throw(
            ArgumentError(
                "$(iv) has no interval, so the number $(arg) in the place of $(iv) in " *
                    "$(_declared_call(v, u_op)) cannot be an endpoint; write a trace on a boundary piece instead"
            )
        )
        T, lo, hi = _float_interval(v, iv)
        val = convert(T, arg)
        ϵ = arg isa AbstractFloat ? max(eps(T), T(eps(typeof(arg)))) : eps(T)
        val == lo && return BoxFace(iv, :lo)
        val == hi && return BoxFace(iv, :hi)
        onlo, onhi = _on_endpoint(val, lo, hi - lo, ϵ), _on_endpoint(val, hi, hi - lo, ϵ)
        if onlo && onhi
            # The tolerances overlap, as on Interval(1e16, 1e16 + 2): the nearer endpoint wins.
            abs(val - lo) < abs(val - hi) && return BoxFace(iv, :lo)
            abs(val - hi) < abs(val - lo) && return BoxFace(iv, :hi)
            throw(
                ArgumentError(
                    "$(arg) in the place of $(iv) in $(_declared_call(v, u_op)) is equally close to both " *
                        "endpoints of $(iv) ∈ [$(lo), $(hi)], which is too narrow to tell its endpoints " *
                        "apart at the precision of $(arg); write an endpoint exactly"
                )
            )
        end
        onlo && return BoxFace(iv, :lo)
        onhi && return BoxFace(iv, :hi)
        near, side = abs(val - lo) <= abs(val - hi) ? (lo, "lower") : (hi, "upper")
        throw(
            ArgumentError(
                "$(arg) in the place of $(iv) in $(_declared_call(v, u_op)) is not on the boundary " *
                    "of $(iv) ∈ [$(lo), $(hi)]: the nearest endpoint is the $(side) endpoint $(near), " *
                    "at distance $(abs(val - near))"
            )
        )
    elseif isequal(arg, iv)
        return Free(iv)
    end
    throw(
        ArgumentError(
            "argument $(arg) in position $(position) of a call of $(u_op) is neither $(iv), an " *
                "endpoint of its interval, nor a boundary piece; the declared argument list is " *
                "$(_declared_call(v, u_op))"
        )
    )
end

_ispiece(arg) = unwrap_const(safe_unwrap(arg)) isa BoundaryPiece

"""
    classify_call(v::VariableMap, call) -> (u_op, Vector{ArgumentClass})

Classify every argument of `call`, a call of a dependent variable of `v`, with
[`classify_argument`](@ref), and return the operation of the call with the classes in
argument order.

A call without a [`BoundaryPiece`](@ref) must have the declared number of arguments, with
each declared variable in its own position or replaced by a number. A call with a piece
must be the declared argument list with every spatial argument replaced by the one piece
and nothing else changed: on `u(t, x, y)` the trace is `u(t, Γ)`, on Brusselator's
`u(x, y, t)` it is `u(Γ, t)`, and on a steady `u(x, y)` it is `u(Γ)`. A trace `u(Γ)` of a
time-dependent variable is an error, "trace is missing the time argument". A trace has
one class fewer than the declaration has arguments for each spatial argument beyond the
first.
"""
function classify_call(v::VariableMap, call)
    call = safe_unwrap(call)
    u_op = operation(call)
    haskey(v.args, u_op) ||
        throw(ArgumentError("$(call) is not a call of a dependent variable; the dependent variables are $(v.ū)"))
    args = arguments(call)
    declared = v.args[u_op]
    npieces = count(_ispiece, args)
    if npieces == 0
        length(args) == length(declared) || throw(
            ArgumentError(
                "$(call) has $(length(args)) arguments, but $(u_op) is declared as $(_declared_call(v, u_op))"
            )
        )
        return u_op, ArgumentClass[classify_argument(v, u_op, i, a) for (i, a) in enumerate(args)]
    end
    npieces == 1 || throw(ArgumentError("$(call) has more than one boundary piece argument"))
    # Declared positions of the trace's arguments, with `0` where the piece stands in for all spatial ones.
    spatial = findall(a -> !isequal(a, v.time), declared)
    isempty(spatial) && throw(
        ArgumentError("$(call) has a boundary piece, but $(_declared_call(v, u_op)) has no spatial arguments")
    )
    layout = [i == first(spatial) ? 0 : i for i in eachindex(declared) if !(i in spatial) || i == first(spatial)]
    piece = unwrap_const(safe_unwrap(args[findfirst(_ispiece, args)]))
    expected = Num(u_op((i == 0 ? piece : declared[i] for i in layout)...))
    if length(args) != length(layout)
        if length(args) == length(layout) - 1 && v.time !== nothing && any(a -> isequal(a, v.time), declared)
            throw(ArgumentError("$(call): trace is missing the time argument; write `$(expected)`"))
        end
        throw(
            ArgumentError(
                "$(call) has $(length(args)) arguments, but a trace of $(_declared_call(v, u_op)) " *
                    "is written $(expected)"
            )
        )
    end
    classes = map(zip(layout, args)) do (i, a)
        (i == 0) == _ispiece(a) || throw(
            ArgumentError(
                "$(call) does not have the arguments of $(_declared_call(v, u_op)) in their declared " *
                    "positions; the trace is written $(expected)"
            )
        )
        classify_argument(v, u_op, i == 0 ? first(spatial) : i, a)
    end
    return u_op, Vector{ArgumentClass}(classes)
end

"""
    location(classes)
    location(v::VariableMap, classes)

Return where a call with argument classes `classes` (from [`classify_call`](@ref)) lies:
[`Interior()`](@ref PDEBase.Interior) when every class is [`Free`](@ref PDEBase.Free),
the one [`BoxFace`](@ref PDEBase.BoxFace), or the one [`Piece`](@ref PDEBase.Piece). Two
faces in one call (an edge or a corner) or a face with a piece are an error.

With `v`, a piece that is a face of the box `spatial_domain(v)` is normalised to its
`BoxFace` with [`boxface`](@ref PDEBase.boxface), so `u(t, 0, y)` and the trace
`u(t, Γ)` on the left face of `(x, y) ∈ Rectangle(...)` have equal locations.
"""
function location(classes)
    bound = filter(c -> !(c isa Free), classes)
    isempty(bound) && return Interior()
    length(bound) == 1 && return only(bound)
    any(c -> c isa Piece, bound) && throw(
        ArgumentError(
            "a trace on $(only(c.piece for c in bound if c isa Piece)) also fixes a non-spatial " *
                "argument: $(bound); write an initial condition over the whole domain, as " *
                "`u(0, x, y)`, and keep the time argument free in a trace, as `u(t, Γ)`"
        )
    )
    throw(
        ArgumentError(
            "a call lies on more than one face of the domain, an edge or a corner: $(bound); " *
                "write the condition on one face"
        )
    )
end

function location(v::VariableMap, classes)
    loc = location(classes)
    loc isa Piece || return loc
    face = boxface(v, loc.piece)
    return face === nothing ? loc : face
end

"""
    boxface(v::VariableMap, piece::BoundaryPiece) -> Union{BoxFace, Nothing}

Return the [`BoxFace`](@ref PDEBase.BoxFace) that `piece` is, when it is a face of the box
`spatial_domain(v)` (a `DomainSets.HyperRectangle`, with pieces from
[`boundarypieces`](@ref)), and `nothing` otherwise. The axes of the box follow the domain
declaration, as in `spatial_domain`: the spatial coordinates of a tuple pairing in tuple
order, or the spatial interval pairings in the order of the domain list, never the order
of `indvars(v)`.

A piece is matched to a face by its name, as `boundarypieces` names the faces, so a
user-built piece with a face's name and a custom geometry, such as part of that face,
normalises to the whole face. A face name beyond the dimension of the box, as `:x5_lo` on
a 3D box, is an `ArgumentError` that lists the valid names.
"""
function boxface(v::VariableMap, piece::BoundaryPiece)
    Ω = spatial_domain(v)
    Ω isa DomainSets.HyperRectangle && isequal(parent(piece), Ω) || return nothing
    face = try
        _box_face(piece)
    catch e
        e isa ArgumentError || rethrow()
        return nothing
    end
    axis, side = face
    axes = _box_axes(v)
    axis <= length(axes) || throw(
        ArgumentError(
            "$(piece) names axis $(axis) of a $(length(axes))-dimensional box; valid face names are " *
                join(map(repr, _box_face_names(length(axes))), ", ")
        )
    )
    return BoxFace(axes[axis], side)
end

"""
    normal(face::BoxFace) -> (iv, sign)
    normal(piece::Piece, x)

Return the outward normal of a location. For a [`BoxFace`](@ref PDEBase.BoxFace) it is
`±e_iv`, returned as the pair `(iv, -1)` on the lower face and `(iv, 1)` on the upper
face. For a [`Piece`](@ref PDEBase.Piece) it is `DomainSets.normal(parent, x)` at the
point `x` when DomainSets knows the normal of the parent domain; otherwise the
discretization backend supplies it, and this is an error.
"""
normal(face::BoxFace) = (face.iv, face.side === :lo ? -1 : 1)

function normal(p::Piece, x)
    Ω = parent(p.piece)
    hasmethod(DomainSets.normal, Tuple{typeof(Ω), typeof(x)}) && return DomainSets.normal(Ω, x)
    throw(
        ArgumentError(
            "the normal of $(p.piece) is not known to PDEBase or DomainSets for a $(typeof(Ω)); " *
                "the discretization backend supplies it"
        )
    )
end

"""
    ConditionClass

What [`classify_condition`](@ref) found a condition to be.

# Fields
- `kind`: `:initial`, `:boundary` or `:periodic`.
- `locations`: The distinct locations of the condition's calls. One `BoxFace(t, :lo)` for
  an initial condition, one [`BoxFace`](@ref PDEBase.BoxFace) or
  [`Piece`](@ref PDEBase.Piece) for a boundary condition, and the lower then the upper
  face of one axis for a periodic condition.
"""
struct ConditionClass
    kind::Symbol
    locations::Vector{Any}
end

# An interior call of a variable declared with the arguments `args`. It is compatible with,
# and adds no location to, a condition on box faces that fix variables it does not have, as
# `q(t, y)` in `u(t, 0, y) ~ q(t, y)`.
struct _FreeCall
    args::Any
end

_compatible(c::_FreeCall, bound) =
    !isempty(bound) && all(l -> l isa BoxFace && !any(a -> isequal(a, l.iv), c.args), bound)

# The locations in `locs`, with the free calls compatible with the others dropped and the
# rest replaced by `Interior()`.
function _resolve_free_calls(locs)
    bound = filter(l -> !(l isa _FreeCall), locs)
    any(l -> l isa _FreeCall && !_compatible(l, bound), locs) && push!(bound, Interior())
    return unique(bound)
end

# Locations of the dependent-variable calls in `ex`, checking that each normal derivative's
# piece is the location of every trace it is applied to. A call of a variable with no
# spatial arguments, such as an ODE variable `w(t)`, has no spatial location and
# contributes none unless it fixes the time, as `w(0)`. Other interior calls are recorded
# as `_FreeCall`s, resolved by `_resolve_free_calls`.
function _call_locations!(locs, v, ex)
    ex = safe_unwrap(ex)
    iscall(ex) || return locs
    op = operation(ex)
    if haskey(v.args, op)
        loc = location(v, classify_call(v, ex)[2])
        if loc isa Interior
            all(a -> isequal(a, v.time), v.args[op]) && return locs
            loc = _FreeCall(v.args[op])
        end
        push!(locs, loc)
        return locs
    elseif op isa NormalDerivative
        opface = location(v, [Piece(op.piece)])
        inner = _call_locations!(Any[], v, only(arguments(ex)))
        all(l -> l isa _FreeCall, inner) &&
            throw(ArgumentError("$(ex) takes a normal derivative of an expression with no trace on $(op.piece)"))
        arglocs = _resolve_free_calls(push!(inner, opface))
        all(loc -> isequal(loc, opface), arglocs) || throw(
            ArgumentError(
                "$(ex) takes the normal derivative on $(op.piece) of an expression located at " *
                    "$(filter(l -> !isequal(l, opface), arglocs)), which is not on $(op.piece); write the trace on the operator's piece"
            )
        )
        push!(locs, opface)
        return locs
    end
    for arg in arguments(ex)
        _call_locations!(locs, v, arg)
    end
    return locs
end

_istimeface(v, loc, side) = loc isa BoxFace && v.time !== nothing && isequal(loc.iv, v.time) && loc.side === side

"""
    classify_condition(v::VariableMap, eq) -> ConditionClass

Classify the condition `eq` (an `Equation` or a `Pair`) by the set of
[`location`](@ref PDEBase.location)s of its dependent-variable calls, and return a
[`ConditionClass`](@ref PDEBase.ConditionClass):

| Locations of the calls | Result |
|:--- |:--- |
| all on the lower time face, as `u(0, x)` | `:initial` |
| any on the upper time face | error, not an initial condition |
| one spatial face or one piece | `:boundary` |
| the two opposite faces of one spatial axis | `:periodic` |
| two pieces that are not opposite faces of a box | error, not supported yet |
| interior calls mixed with traces | error |
| two faces in one call (an edge or a corner) | error |

A trace written with a piece that is a face of the box classifies like its endpoint form.
A call of a variable with no spatial arguments, such as `w(t)`, has no location, so
`u(t, 0) ~ w(t)` is a boundary condition; `w(0)` is on the lower time face. Likewise an
interior call of a variable that does not have the variables its condition's faces fix,
as `q(t, y)` in `u(t, 0, y) ~ q(t, y)`, adds no location, while `u(t, 0, y) ~ u(t, x, y)`
mixes an interior call with a trace. A
`NormalDerivative(Γa)` applied to an expression whose traces are not all on `Γa` is an
error.
"""
function classify_condition(v::VariableMap, eq)
    lhs, rhs = _get_lhs_rhs(eq)
    locs = _resolve_free_calls(_call_locations!(_call_locations!(Any[], v, lhs), v, rhs))
    isempty(locs) && throw(ArgumentError("condition $(eq) has no call of a dependent variable"))
    any(loc -> loc isa Interior, locs) && throw(
        ArgumentError(
            "condition $(eq) is not on the boundary: " *
                (length(locs) == 1 ? "every call is in the interior" : "it mixes interior calls with traces")
        )
    )
    any(loc -> _istimeface(v, loc, :hi), locs) &&
        throw(ArgumentError("condition $(eq) is on the final time $(v.intervals[v.time][2]), not an initial condition"))
    if length(locs) == 1
        loc = only(locs)
        return ConditionClass(_istimeface(v, loc, :lo) ? :initial : :boundary, locs)
    elseif length(locs) == 2 && all(loc -> loc isa BoxFace && !_istimeface(v, loc, :lo), locs)
        lo, hi = sort(locs; by = loc -> loc.side === :hi)
        isequal(lo.iv, hi.iv) && lo.side !== hi.side && return ConditionClass(:periodic, [lo, hi])
    elseif length(locs) == 2 && all(loc -> loc isa Piece, locs)
        throw(
            ArgumentError(
                "condition $(eq) relates the pieces $(locs[1].piece) and $(locs[2].piece), which are not " *
                    "opposite faces of a box; periodic conditions between such pieces need the map " *
                    "between them, inferred as a translation (marked INFERRED in the spec), which is not supported yet"
            )
        )
    end
    throw(
        ArgumentError(
            "condition $(eq) relates the locations $(locs); a condition must be on one face or piece, " *
                "or relate the two opposite faces of one axis"
        )
    )
end

### SUBSTITUTION THROUGH OPERATORS ###

_trace_substitute_filter(ex) =
    SymbolicUtils.default_substitute_filter(ex) ||
    (iscall(ex) && operation(ex) isa Union{Differential, NormalDerivative})

"""
    substitute_traces(ex, rules)

`substitute(ex, rules)` that also substitutes inside `Differential` and
[`NormalDerivative`](@ref) applications, which plain `substitute` leaves alone. `ex` may be
an expression or an `Equation`. Use it wherever traces are rewritten, so that
`Dx(u(t, Γ))`, `Dn(u(t, Γ))` and `k * Dn(u(t, Γ)) + α * u(t, Γ)` all have their trace
replaced, and check the result with [`assert_no_traces`](@ref).
"""
substitute_traces(ex, rules) = substitute(ex, rules; filterer = _trace_substitute_filter)
substitute_traces(eq::Equation, rules) = substitute_traces(eq.lhs, rules) ~ substitute_traces(eq.rhs, rules)

function _find_trace(ex)
    ex = safe_unwrap(ex)
    iscall(ex) || return nothing
    operation(ex) isa NormalDerivative && return ex
    any(_ispiece, arguments(ex)) && return ex
    for arg in arguments(ex)
        found = _find_trace(arg)
        found === nothing || return found
    end
    return nothing
end

"""
    assert_no_traces(ex)

Throw an `ArgumentError` if `ex`, an expression or an `Equation`, still holds a trace (a
call with a [`BoundaryPiece`](@ref) argument) or a [`NormalDerivative`](@ref), and return
`ex` otherwise. Run it after replacing traces and before `expand_derivatives`, which
silently turns `Dx(u(t, Γ))` into `0`.
"""
function assert_no_traces(ex)
    found = _find_trace(ex)
    found === nothing && return ex
    throw(
        ArgumentError(
            "$(found) remains in $(ex) after trace substitution; traces and normal derivatives " *
                "must be replaced before the expression is expanded"
        )
    )
end
assert_no_traces(eq::Equation) = (assert_no_traces(eq.lhs); assert_no_traces(eq.rhs); eq)
