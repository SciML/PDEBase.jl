"""
    pde_substitute(expr, dict)

Custom substitute function for PDE discretization that properly handles Differential operators.
In Symbolics v7, `Differential(constant)` throws an error, and the standard SymbolicUtils
Substituter doesn't respect filter functions when processing Differential operators.

This implementation manually traverses the expression tree and:
- Substitutes variables according to the dict
- Does NOT substitute inside `Differential(x)` operators (the differential variable stays symbolic)
- DOES substitute inside the arguments of derivative calls like `Dx(u(t,x))`

This allows expressions like `Dx(u(t, boundary_value))` to work correctly.
"""
@inline function pde_substitute(expr, dict)
    return _pde_sub(safe_unwrap(expr), Dict(safe_unwrap(k) => safe_unwrap(v) for (k, v) in pairs(dict)))
end

const _PDE_SUB_DEBUG = Ref(false)

function _pde_sub(expr, dict::Dict)
    # If expr is a Differential operator itself, don't substitute it
    # This prevents Differential(x) from becoming Differential(constant)
    if expr isa Differential
        return expr
    end

    # Check if expr is directly in the dict
    # Important: check ALL keys for a match, including derivative expressions
    for (k, v) in pairs(dict)
        if isequal(expr, k)
            # Found exact match - return the substitution value directly
            if _PDE_SUB_DEBUG[] && iscall(expr) && operation(expr) isa Differential
                println("DEBUG _pde_sub: Found derivative match")
                println("  expr: ", expr)
                println("  returning: ", v)
            end
            return v
        end
    end

    # If it's a call expression, recurse into it
    if iscall(expr)
        op = operation(expr)
        args = arguments(expr)

        # Handle operation - if it's a Differential, don't substitute inside it
        if op isa Differential
            if _PDE_SUB_DEBUG[]
                println("DEBUG _pde_sub: Processing Differential call (no match found in dict)")
                println("  expr: ", expr)
                println("  dict keys: ", collect(keys(dict)))
            end
            # Substitute only in the arguments (the stuff inside Dx(...))
            # but not in the Differential operator itself
            new_args = [_pde_sub(arg, dict) for arg in args]
            return op(new_args...)
        else
            # For other operations, substitute in both operation and arguments
            new_op = _pde_sub(op, dict)
            new_args = [_pde_sub(arg, dict) for arg in args]
            # Reconstruct the term
            return maketerm(typeof(expr), new_op, new_args, metadata(expr))
        end
    end

    # Base case - return as-is (numbers, symbols not in dict, etc.)
    return expr
end

"""
    count_differentials(term, x)

Counts the number of `Differential` operators with respect to variable `x` in a term.
This is used to determine the order of a PDE.

# Arguments
- `term`: A symbolic expression to analyze
- `x`: The symbolic variable to count differentials for

# Returns
The total count of differential operators with respect to `x` in the term.
"""
function count_differentials(term, x)::Int
    if !iscall(term)
        return 0
    end
    op = operation(term)
    count_children = 0
    for arg in arguments(term)
        count_children += count_differentials(arg, x)
    end
    if op isa Differential && isequal(op.x, x)
        # In Symbolics v7, Differential has an `order` field for higher-order derivatives
        # e.g., Differential(x)^2 creates Differential(x, 2) with order=2
        return op.order + count_children
    end
    return count_children
end

"""
    differential_order(eq, x)

Returns a set of all differential orders present in the equation with respect to variable `x`.

# Arguments
- `eq`: A symbolic equation or expression to analyze
- `x`: The symbolic variable to compute differential orders for

# Returns
A set of integers representing all differential orders found in the equation (excluding zero).
"""
function differential_order(eq, x)
    orders = Set{Int}()
    _differential_order!(orders, eq, x)
    return filter(!iszero, orders)
end

# Internal helper to avoid allocations from closures
function _differential_order!(orders::Set{Int}, eq, x)
    if iscall(eq)
        op = operation(eq)
        if op isa Differential
            push!(orders, count_differentials(eq, x))
        else
            for ch in arguments(eq)
                _differential_order!(orders, ch, x)
            end
        end
    end
    return nothing
end

"""
    has_derivatives(term)

Determines whether a term contains any `Differential` operators.

# Arguments
- `term`: A symbolic expression to check

# Returns
`true` if the term contains any derivatives, `false` otherwise.
"""
function has_derivatives(term)::Bool
    if iscall(term)
        op = operation(term)
        if op isa Differential
            return true
        end
        for arg in arguments(term)
            if has_derivatives(arg)
                return true
            end
        end
        return false
    else
        return false
    end
end

"""
    find_derivative(term, depvar_op)

Finds the first `Differential` operator or dependent variable matching `depvar_op` within a term.

# Arguments
- `term`: A symbolic expression to search
- `depvar_op`: The operation (function symbol) of the dependent variable to find

# Returns
The first matching derivative or dependent variable found, or `nothing` if none exists.
"""
function find_derivative(term, depvar_op)
    if iscall(term)
        op = operation(term)
        if (op isa Differential) || isequal(op, depvar_op)
            return term
        end
        for arg in arguments(term)
            res = find_derivative(arg, depvar_op)
            if res !== nothing
                return res
            end
        end
    end
    return nothing
end

"""
    subs_alleqs!(eqs, bcs, rules)

Applies substitution rules to all equations and boundary conditions in place.

# Arguments
- `eqs`: Vector of equations to modify
- `bcs`: Vector of boundary conditions to modify
- `rules`: Substitution rules to apply
"""
function subs_alleqs!(eqs, bcs, rules)
    subs_alleqs!(eqs, rules)
    return subs_alleqs!(bcs, rules)
end

function subs_alleqs!(eqs, rules)
    return map!(eq -> substitute(eq.lhs, rules) ~ substitute(eq.rhs, rules), eqs, eqs)
end

"""
    get_depvars(eq, depvar_ops)

Finds all dependent variables matching the given operations in an expression.

# Arguments
- `eq`: A symbolic expression to search
- `depvar_ops`: Collection of dependent variable operations to match

# Returns
A set of all dependent variables found in the expression that match the given operations.
"""
function get_depvars(eq, depvar_ops)
    depvars = Set()
    _get_depvars!(depvars, safe_unwrap(eq), depvar_ops)
    return depvars
end

# Internal helper to avoid allocations from closures and intermediate sets
function _get_depvars!(depvars::Set, eq, depvar_ops)
    if iscall(eq)
        op = operation(eq)
        found = false
        for u in depvar_ops
            if isequal(op, u)
                found = true
                break
            end
        end
        if found
            push!(depvars, eq)
        else
            for arg in arguments(eq)
                _get_depvars!(depvars, arg, depvar_ops)
            end
        end
    end
    return nothing
end

function get_indvars(eq, v)
    ivs = Set()
    _get_indvars!(ivs, safe_unwrap(eq), v)
    return ivs
end

# Internal helper to avoid allocations from closures and intermediate sets
function _get_indvars!(ivs::Set, eq, v)
    if iscall(eq)
        for arg in arguments(eq)
            _get_indvars!(ivs, arg, v)
        end
    else
        for x in v.x̄
            if isequal(eq, x)
                push!(ivs, eq)
                break
            end
        end
    end
    return nothing
end

# Helper function to get left and right sides from Equation or Pair
# In MTK v11, initial conditions can be specified as Pairs (=>) instead of Equations (~)
@inline function _get_lhs_rhs(x::Equation)
    return (x.lhs, x.rhs)
end
@inline function _get_lhs_rhs(x::Pair)
    return (x.first, x.second)
end

@inline function get_all_depvars(pdeeqs, depvar_ops)
    return collect(
        mapreduce(x -> get_depvars(_get_lhs_rhs(x)[1], depvar_ops), union, pdeeqs) ∪
            mapreduce(x -> get_depvars(_get_lhs_rhs(x)[2], depvar_ops), union, pdeeqs)
    )
end

@inline function get_all_depvars(pdesys::PDESystem, depvar_ops)
    pdeeqs = get_eqs(pdesys)
    return collect(
        mapreduce(x -> get_depvars(x.lhs, depvar_ops), union, pdeeqs) ∪
            mapreduce(x -> get_depvars(x.rhs, depvar_ops), union, pdeeqs)
    )
end

get_ops(depvars) = map(u -> operation(safe_unwrap(u)), depvars)

function split_terms(eq::Equation)
    lhs = _split_terms(eq.lhs)
    rhs = _split_terms(eq.rhs)
    return vcat(lhs, rhs)
end

function split_terms(eq::Pair)
    # Handle initial conditions specified as Pairs (u(0,x) => value)
    lhs = _split_terms(eq.first)
    rhs = _split_terms(eq.second)
    return vcat(lhs, rhs)
end

function _split_terms(term)
    # TODO: Update this to be exclusive of derivatives and depvars rather than inclusive of +-/*
    if iscall(term)
        op = operation(term)
        if op === (+) || op === (-) || op === (*) || op === (/)
            result = Any[]
            for arg in arguments(term)
                append!(result, _split_terms(arg))
            end
            return result
        end
    end
    return [term]
end

# Additional handling to get around limitations in rules
# Splits out derivatives from containing math expressions for ingestion by the rules
function _split_terms(term, x̄)
    # TODO: Update this to handle more ops e.g. exp sin tanh etc.
    # TODO: Handle cases where two nonlinear laplacians are multiplied together
    if !iscall(term)
        return [term]
    end

    op = operation(term)
    args = collect(arguments(term))  # Collect to mutable array for v4 compatibility

    # Additional handling for upwinding
    if op === (*)
        for (i, arg) in enumerate(args)
            # Incase of upwinding, we need to keep the original term
            if iscall(arg) && operation(arg) isa Differential
                # Flatten the arguments of the differential to make nonlinear laplacian work in more cases
                try
                    args[i] = operation(arg)(flatten_division.(arguments(arg))...)
                catch e
                    println("Argument to derivative in $term is not a dependant variable, is trivially differentiable or is otherwise not differentiable.")
                    throw(e)
                end
                return [*(flatten_division.(args)...)]
            end
        end
        result = Any[]
        for arg in args
            append!(result, _split_terms(arg, x̄))
        end
        return result
    elseif op === (/)
        # Incase of upwinding or spherical, we need to keep the original term
        if iscall(args[1])
            if args[1] isa Differential
                try
                    args[1] = operation(args[1])(flatten_division.(arguments(args[1]))...)
                catch e
                    println("Argument to derivative in $term is not a dependant variable, is trivially differentiable or is otherwise not differentiable.")
                    throw(e)
                end
                return [/(flatten_division.(args)...)]
                # Handle with care so that spherical still works
            elseif operation(args[1]) === (*)
                subargs = collect(arguments(args[1]))  # Collect to mutable array
                # look for a differential in the arguments
                for (i, arg) in enumerate(subargs)
                    if iscall(arg) && operation(arg) isa Differential
                        # Flatten the arguments of the differential to make nonlinear laplacian/spherical work in more cases
                        try
                            subargs[i] = operation(arg)(flatten_division.(arguments(arg))...)
                            args[1] = operation(args[1])(flatten_division.(subargs)...)
                        catch e
                            println("Argument to derivative in $term is not a dependant variable, is trivially differentiable or is otherwise not differentiable.")
                            throw(e)
                        end
                        return [/(flatten_division.(args)...)]
                    end
                end
            end
        end
        # Basecase for division
        return vcat(_split_terms(args[1], x̄), _split_terms(args[2], x̄))
    elseif op === (+) || op === (-)
        result = Any[]
        for arg in args
            append!(result, _split_terms(arg, x̄))
        end
        return result
    elseif op isa Differential
        return [op(flatten_division.(arguments(term))...)]
    else
        return [term]
    end
end

function split_terms(eq::Equation, x̄)
    lhs = _split_terms(eq.lhs, x̄)
    rhs = _split_terms(eq.rhs, x̄)
    return filter(term -> !isequal(term, Num(0)), flatten_division.(vcat(lhs, rhs)))
end

function split_additive_terms(eq)
    # Calling the methods from symbolicutils matches the expressions
    rhs_arg = iscall(eq.rhs) && (SymbolicUtils.operation(eq.rhs) == +) ?
        SymbolicUtils.arguments(eq.rhs) : [eq.rhs]
    lhs_arg = iscall(eq.lhs) && (SymbolicUtils.operation(eq.lhs) == +) ?
        SymbolicUtils.arguments(eq.lhs) : [eq.lhs]

    return vcat(lhs_arg, rhs_arg)
end

# Filthy hack to get around limitations in rules and avoid simplification to a dividing expression
@inline function flatten_division(term)
    #=rules = [@rule(/(~a, ~b) => *(~a, b^(-1.0))),
             @rule(/(*(~~a), ~b) => *(~a..., b^(-1.0))),
             @rule(/(~a, *(~~b)) => *(~a, *(~b...)^(-1.0))),
             @rule(/(*(~~a), *(~~b)) => *(~a..., *(~b...)^(-1.0)))]
    for r in rules
        if r(term) !== nothing
            return r(term)
        end
    end=#
    return term
end

subsmatch(eq::Equation, rule) = subsmatch(eq.lhs, rule) | subsmatch(eq.rhs, rule)

function subsmatch(expr, rule)
    if isequal(expr, rule.first)
        return true
    end
    if iscall(expr)
        return any(ex -> subsmatch(ex, rule), arguments(expr))
    end
    return false
end
#substitute(eq::Equation, rules) = substitute(eq.lhs, rules) ~ substitute(eq.rhs, rules)

"""
    ex2term(x::Term) -> Symbolic
    ex2term(x) -> x

Convert a Term to a variable `Term`. Note that it only takes a `Term`
not a `Num`.
```
"""
function ex2term(term, v)
    iscall(term) || return term
    termdvs = collect(get_depvars(term, v.depvar_ops))
    symdvs = filter(u -> all(x -> !(unwrap_const(safe_unwrap(x)) isa Number), arguments(u)), termdvs)
    exdv = last(sort(symdvs, by = u -> length(arguments(u))))
    name = Symbol("⟦" * string(term) * "⟧")
    return setname(maketerm(typeof(exdv), rename(operation(exdv), name), arguments(exdv), metadata(exdv)), name)
end

safe_unwrap(x) = x isa Num ? unwrap(x) : x

function recursive_unwrap(ex)
    if !iscall(ex)
        return safe_unwrap(ex)
    end

    op = operation(ex)
    args = arguments(ex)
    return safe_unwrap(op(recursive_unwrap.(args)...))
end

hascomplex(eq::Equation) = hascomplex(eq.lhs) || hascomplex(eq.rhs)
hascomplex(eq::Pair) = hascomplex(eq.first) || hascomplex(eq.second)
# Complex{T} types are obviously complex
hascomplex(term::Complex) = true
hascomplex(term) = !isequal(term, real(term))

split_complex(eq::Vector) = eq
split_complex(eq::Equation) = split_complex(eq.lhs) .~ split_complex(eq.rhs)
function split_complex(term)
    if hascomplex(term)
        return [real(term), imag(term)]
    else
        return [term, term]
    end
end
