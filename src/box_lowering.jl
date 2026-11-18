### BOX LOWERING ###

# Rewrites trace conditions on the faces of a box into the endpoint form that box
# discretizations parse. Locations come from the classifier and every rewrite goes through
# `substitute_traces`; nothing here derivative-expands a boundary equation.

# `VariableMap` takes the time variable from a discretization; this carries only that.
struct _BoxTime
    time::Any
end
get_time(d::_BoxTime) = d.time

function _check_box_domains(pdesys)
    for pairing in get_domain(pdesys)
        d = pairing.domain
        kind = d isa MeshDomain ? "MeshDomain" :
            d isa DomainSets.Domain && !(d isa DomainSets.AbstractInterval || _is_interval_product(d)) ?
            string(nameof(typeof(d))) : nothing
        kind === nothing || throw(
            ArgumentError(
                "MethodOfLines needs an interval or box domain; this system uses a $(kind), use FEMDiscretization"
            )
        )
    end
    return nothing
end

_has_trace(bc::Union{Equation, Pair}) = any(side -> _find_trace(side) !== nothing || _has_facet_normal(side), _get_lhs_rhs(bc))
_has_trace(bcs::AbstractVector) = any(_has_trace, bcs)
_has_trace(_) = false

# Every call with a `BoundaryPiece` argument in the conditions `bc`.
function _trace_calls!(calls, ex)
    ex = safe_unwrap(ex)
    iscall(ex) || return calls
    any(_ispiece, arguments(ex)) && push!(calls, ex)
    foreach(a -> _trace_calls!(calls, a), arguments(ex))
    return calls
end
_trace_calls!(calls, bc::Union{Equation, Pair}) = foldl(_trace_calls!, _get_lhs_rhs(bc); init = calls)
_trace_calls!(calls, bcs::AbstractVector) = foldl(_trace_calls!, bcs; init = calls)

function _check_trace_types(calls)
    for call in calls
        T = SymbolicUtils.symtype(call)
        T <: Complex || T <: AbstractArray || continue
        throw(
            ArgumentError(
                "$(call) is a trace of a $(T <: Complex ? "complex" : "array") dependent variable; traces " *
                    "on array or complex dependent variables are not supported yet, write the endpoint form"
            )
        )
    end
    return nothing
end

_or_nothing(f) =
try
    f()
catch e
    e isa Union{ArgumentError, ErrorException} || rethrow()
    nothing
end

# The pieces of the traces and of the normal-derivative operators in `ex`.
function _pieces!(pieces, ex)
    ex = safe_unwrap(ex)
    iscall(ex) || return pieces
    operation(ex) isa Union{NormalDerivative, FacetNormal} && push!(pieces, operation(ex).piece)
    for a in arguments(ex)
        _ispiece(a) ? push!(pieces, unwrap_const(safe_unwrap(a))) : _pieces!(pieces, a)
    end
    return pieces
end

# The conditions holding a trace or a normal derivative, out of nested condition vectors.
_trace_conditions(bc::Union{Equation, Pair}) = _has_trace(bc) ? Any[bc] : Any[]
_trace_conditions(bcs::AbstractVector) = reduce(vcat, map(_trace_conditions, bcs); init = Any[])
_trace_conditions(_) = Any[]

# The time variable, when the caller gives none: the independent variable (or `nothing`)
# whose removal makes the spatial domain the parent of every piece in the conditions, as the
# classifier sees it. Two candidates can fit the domain when intervals coincide, as `t` and
# `x` on `Interval(0, 1)`; then the one under which every condition classifies wins, and after
# that the one declared by a pairing of its own.
function _infer_time(pdesys, conds)
    pieces = unique(foldl((ps, bc) -> foldl(_pieces!, _get_lhs_rhs(bc); init = ps), conds; init = Any[]))
    function fit(τ)
        v = _or_nothing(() -> VariableMap(pdesys, _BoxTime(τ)))
        Ω = v === nothing ? nothing : _or_nothing(() -> spatial_domain(v))
        Ω !== nothing && all(p -> isequal(parent(p), Ω), pieces) || return 0
        return _or_nothing(() -> (foreach(bc -> classify_condition(v, bc), conds); true)) === nothing ? 1 : 2
    end
    candidates = Any[nothing; map(safe_unwrap, get_ivs(pdesys))]
    scores = map(fit, candidates)
    best = maximum(scores)
    best == 0 && return nothing
    winners = candidates[scores .== best]
    length(winners) == 1 && return only(winners)
    # Still tied: prefer the one variable with a pairing of its own, as `t ∈ Interval(0, 1)`.
    alone = filter(τ -> any(p -> isequal(safe_unwrap(p.variables), τ), get_domain(pdesys)), winners)
    length(alone) == 1 && return only(alone)
    throw(
        ArgumentError(
            "the time variable cannot be inferred from the traces: each of $(winners) fits the " *
                "domain; pass it as to_box_bcs(pdesys; time = t)"
        )
    )
end

function _trace_face(v, piece)
    face = boxface(v, piece)
    face === nothing && throw(
        ArgumentError(
            "$(piece) is not a face of the box $(spatial_domain(v)); to_box_bcs lowers traces " *
                "on the faces from boundarypieces only"
        )
    )
    return face
end

# The trace `call` with its piece replaced by the declared spatial arguments, and the
# coordinate of `face` fixed at its endpoint.
function _endpoint_call(v, call, face)
    u_op = operation(call)
    any(a -> isequal(a, face.iv), v.args[u_op]) || throw(
        ArgumentError(
            "$(call) is a trace on the face where $(face.iv) is fixed, but $(_declared_call(v, u_op)) " *
                "does not depend on $(face.iv); write the condition on a face of its own coordinates"
        )
    )
    value = v.intervals[face.iv][face.side === :lo ? 1 : 2]
    kept = Iterators.Stateful(filter(a -> !_ispiece(a), arguments(call)))
    args = map(v.args[u_op]) do a
        isequal(a, v.time) ? popfirst!(kept) : isequal(a, face.iv) ? value : a
    end
    return u_op(args...)
end

# With `declared`, every call of a dependent variable becomes the call with its declared variables,
# for the derivative check only; the lowered output never uses it.
function _box_rules!(rules, v, ex; declared = false)
    ex = safe_unwrap(ex)
    iscall(ex) || return rules
    op = operation(ex)
    op isa Union{Differential, NormalDerivative} && _check_normal_not_differentiated(ex)
    if op isa FacetNormal
        iv, sign = normal(_trace_face(v, op.piece))
        rules[ex] = isequal(indvars(v)[_normal_component(ex)], iv) ? sign : 0
        return rules
    elseif op isa NormalDerivative
        iv, sign = normal(_trace_face(v, op.piece))
        operand = safe_unwrap(only(arguments(ex)))
        derivative = Differential(iv)(substitute_traces(operand, _box_rules!(Dict{Any, Any}(), v, operand; declared)))
        rules[ex] = sign < 0 ? -derivative : derivative
        return rules
    elseif declared && haskey(v.args, op)
        rules[ex] = op(v.args[op]...)
        return rules
    elseif haskey(v.args, op) && any(_ispiece, arguments(ex))
        piece = unwrap_const(safe_unwrap(only(filter(_ispiece, arguments(ex)))))
        rules[ex] = _endpoint_call(v, ex, _trace_face(v, piece))
        return rules
    end
    foreach(a -> _box_rules!(rules, v, a; declared), arguments(ex))
    return rules
end

# Does `ex` hold a call of a dependent variable of `v`?
_holds_field(v, ex) = (ex = safe_unwrap(ex); iscall(ex) && (haskey(v.args, operation(ex)) || any(a -> _holds_field(v, a), arguments(ex))))

# The derivative chains of dependent variables left in `ex` that have a spatial derivative, outermost
# first; a derivative of data is none.
function _spatial_derivatives!(nodes, ex, v)
    ex = safe_unwrap(ex)
    iscall(ex) || return nodes
    if operation(ex) isa Differential
        chain, inner = Any[], ex
        while iscall(inner) && operation(inner) isa Differential
            push!(chain, operation(inner).x)
            inner = safe_unwrap(only(arguments(inner)))
        end
        if any(x -> any(y -> isequal(x, y), indvars(v)), chain) && _holds_field(v, inner)
            any(n -> isequal(n, ex), nodes) || push!(nodes, ex)
            return nodes
        end
    end
    foreach(a -> _spatial_derivatives!(nodes, a, v), arguments(ex))
    return nodes
end

# The derivative chain `ex` with its variables in the order of the system's variables, as Symbolics does
# not identify the mixed partials `Dx(Dy(f))` and `Dy(Dx(f))`.
function _canonical_derivative(v, ex)
    vars, inner = Any[], safe_unwrap(ex)
    while iscall(inner) && operation(inner) isa Differential
        append!(vars, fill(operation(inner).x, operation(inner).order))
        inner = safe_unwrap(only(arguments(inner)))
    end
    order = vcat(Any[v.time], indvars(v))
    position(x) = something(findfirst(y -> isequal(x, y), order), length(order) + 1)
    sort!(vars; by = position)
    return foldr((x, e) -> Differential(x)(e), vars; init = inner)
end

# A condition with a facet normal must keep a derivative of a field: the normal is allowed only in a
# natural or Robin-type condition. The fields are put in their declared variables and the derivatives
# are really expanded (product and chain rules), for this check only, so a derivative that is zero
# (in a parameter, or of a field lacking the variable) does not count and one that survives a
# coefficient (`Dy(x * w)`) does. At the face the endpoint is substituted into the coefficients; a
# residual that is then zero but held a derivative is the zero natural condition of a degenerate face.
function _check_box_normal_derivative(v, lhs, rhs, bc)
    rules = _box_rules!(_box_rules!(Dict{Any, Any}(), v, lhs; declared = true), v, rhs; declared = true)
    residual = Symbolics.expand_derivatives(Num(substitute_traces(lhs, rules)) - Num(substitute_traces(rhs, rules)))
    nodes = _spatial_derivatives!(Any[], residual, v)
    canonical = [_canonical_derivative(v, n) for n in nodes]
    keys = unique(canonical)
    taken = Symbolics.get_variables(residual)
    function placeholder(j)
        name = "##derivative#$j"
        while any(var -> isequal(var, unwrap(Symbolics.variable(Symbol(name)))), taken)
            name *= "'"
        end
        return Symbolics.variable(Symbol(name))
    end
    slots = [placeholder(j) for j in eachindex(keys)]
    replacements = Dict(unwrap(n) => unwrap(slots[findfirst(k -> isequal(k, c), keys)]) for (n, c) in zip(nodes, canonical))
    residual = Symbolics.simplify(Symbolics.expand(substitute_traces(residual, replacements)))
    holds_derivative(r) = (vars = Set(unwrap.(Symbolics.get_variables(r))); any(g -> unwrap(g) in vars, slots))
    written = holds_derivative(residual)
    for piece in unique(_pieces!(_pieces!(Any[], lhs), rhs))
        face = boxface(v, piece)
        face === nothing || (residual = Symbolics.substitute(residual, Dict(face.iv => v.intervals[face.iv][face.side === :lo ? 1 : 2])))
    end
    _check_finite_constants(residual, bc)
    residual = Symbolics.simplify(Symbolics.expand(residual))
    return _check_normal_has_derivative(true, holds_derivative(residual), isequal(residual, 0) && written, bc)
end

function _lower_bc(v, bc::Union{Equation, Pair})
    _has_trace(bc) || return bc
    classify_condition(v, bc)
    lhs, rhs = _get_lhs_rhs(bc)
    rules = _box_rules!(_box_rules!(Dict{Any, Any}(), v, lhs), v, rhs)
    lowered = map(side -> assert_no_traces(substitute_traces(side, rules)), (lhs, rhs))
    if _has_facet_normal(lhs) || _has_facet_normal(rhs)
        for side in lowered
            _check_finite_constants(side, bc)
            for piece in unique(_pieces!(_pieces!(Any[], lhs), rhs))
                face = boxface(v, piece)
                face === nothing || _check_finite_constants(Symbolics.substitute(side, Dict(face.iv => v.intervals[face.iv][face.side === :lo ? 1 : 2])), bc)
            end
        end
        _check_box_normal_derivative(v, lhs, rhs, bc)
    end
    return bc isa Pair ? Pair(lowered...) : Equation(lowered...)
end
_lower_bc(v, bcs::AbstractVector) = map(bc -> _lower_bc(v, bc), bcs)
_lower_bc(v, bc) = bc

"""
    to_box_bcs(pdesys::PDESystem; time = <inferred>) -> PDESystem

Rewrite every trace condition of `pdesys` on the faces of an interval or box into the
endpoint form that box discretizations such as MethodOfLines parse, and leave everything
else untouched. For the face of axis `x` at its lower end `x_lo` (from
[`boundarypieces`](@ref), for example `:left` of a `Rectangle` or `:x1_lo` of an interval):

| Trace form | Lowered form |
|:--- |:--- |
| `u(t, Γ) ~ g` | `u(t, x_lo, y) ~ g` |
| `Dn_Γ(u(t, Γ)) ~ h` | `-Dx(u(t, x_lo, y)) ~ h`, and `+Dx(...)` on an upper face |
| `k * Dn_Γ(u(t, Γ)) + α * u(t, Γ) ~ β` | `-k * Dx(u(t, x_lo, y)) + α * u(t, x_lo, y) ~ β` |
| `u(t, Γ_lo) ~ u(t, Γ_hi)` | `u(t, x_lo, y) ~ u(t, x_hi, y)` |
| `Dx(u(t, Γ)) ~ h` | `Dx(u(t, x_lo, y)) ~ h` |
| `k * (n[1] * Dx(u(t, Γ)) + n[2] * Dy(u(t, Γ))) + α * u(t, Γ) ~ β`, `n = normal(Γ)` | `-k * Dx(u(t, x_lo, y)) + α * u(t, x_lo, y) ~ β` |

A component of [`normal`](@ref PDEBase.normal)`(Γ)` is replaced by the numeric `-1`, `+1` or `0` of the
face's outward normal, in the order of `indvars`.

The endpoint value is the interval's own endpoint, `DomainSets.infimum` or `supremum`, so
`0.0` on `Interval(0.0, 1.0)`; it is `isequal` to the `0` a user writes. Conditions that
hold no trace, including existing endpoint-form conditions, are passed through as they
are, so the two styles may be mixed. A system without traces is returned as the same
object. Every other field of the system is kept. Each lowered condition is checked with
[`classify_condition`](@ref) and [`assert_no_traces`](@ref) and is never
derivative-expanded.

`time` is the time variable. When it is not given, it is the independent variable (or
`nothing`, for a steady system) that leaves the parent domain of the traced pieces as the
spatial domain; when two fit, as `t` and `x` both on `Interval(0, 1)`, the one under which
every trace classifies is taken, then the one declared by a pairing of its own, and
otherwise `time` must be passed. A given `time` is used as is. `lower_boundaries` calls `to_box_bcs` with the discretization's time for
every `AbstractBoxDiscretization`. Conditions may be `Equation`s or `Pair`s.

Throws an `ArgumentError` when a domain of `pdesys` is a `MeshDomain` or another domain
that is not an interval or interval product, such as a ball ("MethodOfLines needs an
interval or box domain; this system uses a MeshDomain, use FEMDiscretization"), and with
the classifier's message for a trace missing its time argument, a periodic pair that is
not the two opposite faces of one axis, or a normal derivative applied to a trace on
another piece. Traces of array or complex dependent variables are not supported yet and are
an `ArgumentError` asking for the endpoint form, as is a trace on a face whose coordinate
the variable does not depend on.
"""
function to_box_bcs(pdesys::PDESystem; time = nothing)
    _check_box_domains(pdesys)
    bcs = get_bcs(pdesys)
    any(_has_trace, bcs) || return pdesys
    calls = _trace_calls!(Any[], bcs)
    _check_trace_types(calls)
    time === nothing && (time = _infer_time(pdesys, _trace_conditions(bcs)))
    v = VariableMap(pdesys, _BoxTime(time))
    return PDESystem(
        get_eqs(pdesys), map(bc -> _lower_bc(v, bc), bcs), get_domain(pdesys), get_ivs(pdesys),
        get_dvs(pdesys), get_ps(pdesys);
        initial_conditions = pdesys.initial_conditions, systems = get_systems(pdesys),
        connector_type = get_connector_type(pdesys), metadata = get_metadata(pdesys),
        analytic = getfield(pdesys, :analytic), analytic_func = getfield(pdesys, :analytic_func),
        gui_metadata = get_gui_metadata(pdesys), description = getfield(pdesys, :description),
        checks = false, name = getfield(pdesys, :name)
    )
end
