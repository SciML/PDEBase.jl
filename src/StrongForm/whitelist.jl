function is_state(ex, depvars)
    iscall(ex) || return false
    any(u -> isequal(operation(ex), operation(unwrap(u))), depvars) && return true
    return any(a -> is_state(a, depvars), arguments(ex))
end

has_differential(ex) = iscall(ex) &&
    (operation(ex) isa Differential || any(has_differential, arguments(ex)))

# Calls `op` on fresh symbols so that derivatives inside sibling arguments cannot be
# mistaken for a missing rule; Bool arguments need no derivative.
function has_derivative_rules(op, args, depvars)
    fresh = [Symbolics.variable(gensym(); T = symtype(a) <: Bool ? Bool : Real) for a in args]
    call = op(fresh...)
    return all(eachindex(args)) do i
        symtype(args[i]) <: Bool || !is_state(args[i], depvars) ||
            !has_differential(unwrap(Symbolics.derivative(call, fresh[i])))
    end
end

unsupported(op) = error("operator $op is not supported by the strong-form analysis")

"""
    check_whitelist(ex, ivs, depvars)

Stage 0 of the strong-form analysis: error unless every operation in `ex` (an expression
or an `Equation`) is supported. Supported are `Differential` in one of `ivs`,
`NormalDerivative` applied to a trace, calls of dependent variables and other symbolic
functions, `ifelse`, and any function Symbolics can differentiate in each of its non-`Bool`
arguments (arithmetic, comparisons and logical operators included). A registered function
without a derivative rule is accepted only while no state from `depvars` reaches it.
Returns `nothing`.
"""
function check_whitelist(ex, ivs, depvars)
    ex = unwrap(ex)
    iscall(ex) || return nothing
    op = operation(ex)
    args = arguments(ex)
    if op isa Differential
        any(iv -> isequal(op.x, iv), unwrap.(ivs)) || unsupported(op)
    elseif nameof(typeof(op)) === :NormalDerivative
        # NormalDerivative is not in Symbolics yet; once it is, also check that the
        # operator's piece is the trace's piece.
        iscall(only(args)) && is_state(only(args), depvars) || unsupported(op)
    elseif op isa Function
        op === ifelse || has_derivative_rules(op, args, depvars) || unsupported(op)
    elseif !(op isa SymbolicUtils.BasicSymbolic)
        unsupported(op)
    end
    foreach(a -> check_whitelist(a, ivs, depvars), args)
    return nothing
end

function check_whitelist(eq::Equation, ivs, depvars)
    check_whitelist(eq.lhs, ivs, depvars)
    return check_whitelist(eq.rhs, ivs, depvars)
end
