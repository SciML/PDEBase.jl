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

_has_trace(bc::Equation) = _find_trace(bc.lhs) !== nothing || _find_trace(bc.rhs) !== nothing
_has_trace(bcs::AbstractVector) = any(_has_trace, bcs)
_has_trace(_) = false

# Independent variables kept as plain arguments of trace calls; by the trace rule they are
# the non-spatial ones, that is time.
function _trace_ivs!(ivs, ex)
    ex = safe_unwrap(ex)
    iscall(ex) || return ivs
    args = arguments(ex)
    if any(_ispiece, args)
        for a in args
            a = unwrap_const(safe_unwrap(a))
            SymbolicUtils.issym(a) && push!(ivs, a)
        end
    end
    foreach(a -> _trace_ivs!(ivs, a), args)
    return ivs
end
_trace_ivs!(ivs, bc::Equation) = _trace_ivs!(_trace_ivs!(ivs, bc.lhs), bc.rhs)
_trace_ivs!(ivs, bcs::AbstractVector) = (foreach(bc -> _trace_ivs!(ivs, bc), bcs); ivs)

function _trace_time(bcs)
    ivs = unique(_trace_ivs!(Any[], bcs))
    length(ivs) <= 1 && return isempty(ivs) ? nothing : only(ivs)
    throw(
        ArgumentError(
            "the traces keep more than one independent variable, $(ivs), so the time variable " *
                "cannot be inferred; pass it as to_box_bcs(pdesys; time = t)"
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
    value = v.intervals[face.iv][face.side === :lo ? 1 : 2]
    kept = Iterators.Stateful(filter(a -> !_ispiece(a), arguments(call)))
    args = map(v.args[u_op]) do a
        isequal(a, v.time) ? popfirst!(kept) : isequal(a, face.iv) ? value : a
    end
    return u_op(args...)
end

function _box_rules!(rules, v, ex)
    ex = safe_unwrap(ex)
    iscall(ex) || return rules
    op = operation(ex)
    if op isa NormalDerivative
        face = _trace_face(v, op.piece)
        iv, sign = normal(face)
        derivative = Differential(iv)(_endpoint_call(v, safe_unwrap(only(arguments(ex))), face))
        rules[ex] = sign < 0 ? -derivative : derivative
    elseif haskey(v.args, op) && any(_ispiece, arguments(ex))
        piece = unwrap_const(safe_unwrap(only(filter(_ispiece, arguments(ex)))))
        rules[ex] = _endpoint_call(v, ex, _trace_face(v, piece))
        return rules
    end
    foreach(a -> _box_rules!(rules, v, a), arguments(ex))
    return rules
end

function _lower_bc(v, bc::Equation)
    _has_trace(bc) || return bc
    classify_condition(v, bc)
    rules = _box_rules!(_box_rules!(Dict{Any, Any}(), v, bc.lhs), v, bc.rhs)
    return assert_no_traces(substitute_traces(bc, rules))
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

The endpoint value is the interval's own endpoint, `DomainSets.infimum` or `supremum`, so
`0.0` on `Interval(0.0, 1.0)`; it is `isequal` to the `0` a user writes. Conditions that
hold no trace, including existing endpoint-form conditions, are passed through as they
are, so the two styles may be mixed. A system without traces is returned as the same
object. Every other field of the system is kept. Each lowered condition is checked with
[`classify_condition`](@ref) and [`assert_no_traces`](@ref) and is never
derivative-expanded.

`time` is the time variable; by default it is the independent variable the traces keep as
an argument, as `t` in `u(t, Γ)`, or `nothing` when they keep none. Pass it when the
traces do not show it. `lower_boundaries` calls `to_box_bcs` with the discretization's time
for every `AbstractBoxDiscretization`.

Throws an `ArgumentError` when a domain of `pdesys` is a `MeshDomain` or another domain
that is not an interval or interval product, such as a ball ("MethodOfLines needs an
interval or box domain; this system uses a MeshDomain, use FEMDiscretization"), and with
the classifier's message for a trace missing its time argument, a periodic pair that is
not the two opposite faces of one axis, or a normal derivative applied to a trace on
another piece.
"""
function to_box_bcs(pdesys::PDESystem; time = _trace_time(get_bcs(pdesys)))
    _check_box_domains(pdesys)
    bcs = get_bcs(pdesys)
    any(_has_trace, bcs) || return pdesys
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
