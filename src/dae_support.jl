# Helpers for building a `DAEProblem` from a discretized system without `mtkcompile`:
# classifying unknowns as differential or algebraic, checking whether
# `BrownFullBasicInit()` is safe for the system's initialization equations, and building
# the operating point.

"""
    BrownFullBasicInitUnsafeError(offenders)

Raised when building a `DAEProblem` from a discretized `PDESystem` whose system carries
initialization equations that `BrownFullBasicInit()` would silently discard, so no
default initialization algorithm can be chosen safely. `offenders` pairs each such
equation with the reason it is not honored.
"""
struct BrownFullBasicInitUnsafeError <: Exception
    offenders::Vector{Pair{Equation, String}}
end

function Base.showerror(io::IO, e::BrownFullBasicInitUnsafeError)
    print(
        io,
        """
        This discretization cannot build a `DAEProblem` for this system without `mtkcompile`.

        `BrownFullBasicInit()` is the only initialization algorithm that reproduces the
        `discretize` (`ODEProblem` + `mtkcompile`) result on the array-equation path. It
        takes the differential variables' initial values as given and solves for the
        algebraic variables and the derivatives, so it would silently ignore these
        initialization equations:
        """
    )
    for (eq, reason) in e.offenders
        print(io, "\n  ", eq, "\n    ", reason)
    end
    print(
        io,
        """


        Use `discretize(pdesys, disc)` instead. That path runs `mtkcompile`, which honors
        these constraints, at the cost of scalarizing the array equations: the slice form
        cannot currently be preserved for this system.
        """
    )
    return
end

_issym(ex) = ex isa SymbolicUtils.BasicSymbolic && iscall(ex)

function _has_differential(ex)
    ex = safe_unwrap(ex)
    _issym(ex) || return false
    operation(ex) isa Differential && return true
    return any(_has_differential, arguments(ex))
end

"""
    time_derivative_order(ex, t)

The order of the highest derivative with respect to `t` appearing in `ex`.
"""
function time_derivative_order(ex, t)
    ex = safe_unwrap(ex)
    _issym(ex) || return 0
    op = operation(ex)
    inner = maximum(a -> time_derivative_order(a, t), arguments(ex); init = 0)
    return op isa Differential && isequal(op.x, t) ? inner + op.order : inner
end


"""
    differential_unknowns(sys, unks)

The unknowns of `sys` that appear under a time derivative: a scalar unknown maps to
`true`, an array unknown to the mask of its differentiated elements. The rest are
algebraic. `unks` is the set of unwrapped unknowns of `sys`.
"""
function differential_unknowns(sys, unks)
    out = Dict{Any, Any}()
    for part in keys(time_differentiated(get_eqs(sys), safe_unwrap(ModelingToolkitBase.get_iv(sys))))
        if part in unks
            out[part] = _is_array_valued(part) ? trues(size(part)) : true
            continue
        end
        unknown, idxs = _parent_and_indices(part)
        (idxs !== nothing && unknown in unks) || continue
        mask = get!(() -> falses(size(unknown)), out, unknown)
        mask[_region(idxs)] .= true
    end
    return out
end

# Whether `x`, an unknown or an element or slice of one, is differential throughout, or
# `nothing` if it is none of these.
function _is_differential(x, diffvars, unks)
    if x in unks
        d = get(diffvars, x, false)
        return d isa Bool ? d : all(d)
    end
    unknown, idxs = _parent_and_indices(x)
    (idxs !== nothing && unknown in unks) || return nothing
    d = get(diffvars, unknown, false)
    return d isa Bool ? d : all(view(d, _region(idxs)))
end

# An expression may mention an unknown, or an element or slice of one.
function _mentions_unknown(ex, unks)
    for v in Symbolics.get_variables(ex)
        v = safe_unwrap(v)
        v in unks && return true
        unknown, idxs = _parent_and_indices(v)
        idxs !== nothing && unknown in unks && return true
    end
    return false
end

"""
    brown_init_offenders(sys)

Return the initialization equations of `sys` that `BrownFullBasicInit()` would not honor,
each paired with the reason, and an empty vector when the algorithm is safe for `sys`.

`BrownFullBasicInit()` holds the differential variables at their given values and solves
for everything else, so an initialization equation survives it only when it fixes
differential unknowns (an unknown, or an element or slice of one) to values involving no
other unknown. Everything else is reported, including equations this predicate cannot
classify.
"""
function brown_init_offenders(sys)
    ieqs = initialization_equations(sys)
    offenders = Pair{Equation, String}[]
    isempty(ieqs) && return offenders
    unks = Set{Any}(safe_unwrap.(get_unknowns(sys)))
    diffvars = differential_unknowns(sys, unks)
    for eq in ieqs
        reason = _brown_init_offence(eq, diffvars, unks)
        reason === nothing || push!(offenders, eq => reason)
    end
    return offenders
end

function _brown_init_offence(eq, diffvars, unks)
    if _has_differential(eq.lhs) || _has_differential(eq.rhs)
        return "constrains a time derivative, which BrownFullBasicInit solves for."
    end
    lhs = safe_unwrap(eq.lhs)
    differential = _is_differential(lhs, diffvars, unks)
    if differential === nothing
        return "does not fix a single unknown of the discretized system, so it cannot be checked against what BrownFullBasicInit holds fixed."
    end
    if !differential
        return "constrains $lhs, which is algebraic rather than differential; BrownFullBasicInit solves for algebraic variables instead of honoring a given value."
    end
    if _mentions_unknown(eq.rhs, unks)
        return "relates $lhs to other unknowns; BrownFullBasicInit fixes each differential unknown independently and does not enforce the relation."
    end
    return nothing
end

"""
    dae_operating_point(sys, u0, t)

The operating point for building a `DAEProblem` from `sys`. `varmap_to_vars` needs a value
for every unknown and every derivative of an unknown. The discretized initial conditions
`u0` cover the grid points of every variable given one; the rest are quantities the
initialization solves for, for which this is only a guess (zero).
"""
function dae_operating_point(sys, u0, t)
    D = Differential(t)
    op = Dict{Any, Any}()
    for x in get_unknowns(sys)
        x = safe_unwrap(x)
        op[x] = _zero_value(x)
        op[D(x)] = _zero_value(x)
    end
    for (k, v) in u0
        op[safe_unwrap(k)] = v
    end
    return op
end

_zero_value(x) = _is_array_valued(x) ? zeros(size(x)) : 0.0

# Each expression differentiated in time in `eqs`, with the highest order it is
# differentiated to.
function time_differentiated(eqs, t)
    out = Dict{Any, Int}()
    for eq in eqs
        _time_differentiated!(out, eq.lhs, t)
        _time_differentiated!(out, eq.rhs, t)
    end
    return out
end

function _time_differentiated!(out, ex, t)
    ex = safe_unwrap(ex)
    _issym(ex) || return out
    if _is_time_derivative(ex, t)
        part, order = _time_derivative_of(ex, t)
        out[part] = max(order, get(out, part, 0))
        return out
    end
    for arg in arguments(ex)
        _time_differentiated!(out, arg, t)
    end
    return out
end

_is_time_derivative(ex, t) = _issym(ex) && operation(ex) isa Differential &&
    isequal(operation(ex).x, t)

# `Dt^k(x)` as `(x, k)`.
function _time_derivative_of(ex, t)
    ex = safe_unwrap(ex)
    order = 0
    while _is_time_derivative(ex, t)
        order += operation(ex).order
        ex = safe_unwrap(only(arguments(ex)))
    end
    return ex, order
end

# The array an element or slice `u[i, 2:n-1]` indexes and its indices, or `x` and
# `nothing` for anything else.
function _parent_and_indices(x)
    iscall(x) && operation(x) === getindex || return x, nothing
    args = arguments(x)
    return safe_unwrap(args[1]), Tuple(unwrap_const(safe_unwrap(i)) for i in args[2:end])
end

_region(idxs) = CartesianIndices(map(i -> i isa Integer ? (i:i) : i, idxs))
