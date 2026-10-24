function is_state(ex, depvars)
    iscall(ex) || return false
    any(u -> isequal(operation(ex), operation(unwrap(u))), depvars) && return true
    return any(a -> is_state(a, depvars), arguments(ex))
end

differentiates_in(ex, s) = iscall(ex) &&
    (operation(ex) isa Differential && isequal(operation(ex).x, s) || any(a -> differentiates_in(a, s), arguments(ex)))

# Probes each state-carrying argument with its own fresh symbol and leaves the others (array
# indices, Bool conditions) in place, so the call keeps its argument shapes.
function has_derivative_rules(op, args, depvars)
    probed = findall(a -> !(symtype(a) <: Bool) && is_state(a, depvars), args)
    fresh = Any[i in probed ? unwrap(Symbolics.variable(gensym())) : args[i] for i in eachindex(args)]
    call = op(Symbolics.wrap.(fresh)...)
    return all(i -> !differentiates_in(unwrap(Symbolics.derivative(call, fresh[i])), fresh[i]), probed)
end

unsupported(op) = error("operator $op is not supported by the strong-form analysis")

"""
    check_whitelist(ex, ivs, depvars)

Stage 0 of the strong-form analysis: error unless every operation in `ex` (an expression
or an `Equation`) is supported. Supported are `Differential` in one of `ivs`,
`NormalDerivative` applied to a trace, calls of dependent variables and other symbolic
functions, `ifelse`, `getindex`, and any function Symbolics can differentiate in each of its non-`Bool`
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
        op === ifelse || op === getindex || has_derivative_rules(op, args, depvars) || unsupported(op)
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
