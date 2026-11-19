### VECTOR OPERATORS ###

# Symbolics 7.49 adds the lazy vector-calculus operators `Grad`, `Divergence`, `Curl` and
# `Laplacian`, which are not `SymbolicUtils.Operator`s and which the rest of PDEBase does not
# recognise. `expand_vector_operators` rewrites them into `Differential` terms before any
# analysis. Older Symbolics lacks them, and the function then returns its argument.

"""
    expand_vector_operators(pdesys::PDESystem) -> PDESystem

Rewrite the vector-calculus operators `Grad`, `Divergence`, `Curl` and `Laplacian` of
Symbolics (7.49 and later) in the equations and boundary conditions of `pdesys` into
`Differential` terms, so that the rest of the pipeline sees only the syntax it already
handles. `parsed_boundary_map` and `symbolic_discretize` call it first.

The operators are applied eagerly, so `∇ ⋅ (k * ∇u)` becomes `Dx(k * Dx(u)) + Dy(k * Dy(u))`
and keeps its divergence form. This also holds for a trace such as `u(t, Γ)`, whose
derivatives `expand_derivatives` would replace by zero. Every variable of an operator must
be an independent variable of `pdesys`, otherwise an `ArgumentError` is thrown. A system
that uses none of the operators is returned as is, and so is every system when the installed
Symbolics lacks them.

The operand of a `Divergence` or `Curl` needs a declared length, since lowering it takes the components
named by the coordinates and must not drop further ones, so an unsized vector, or a composite flux that
contains one (for example a sized matrix times an unsized vector, whose two-element result Symbolics infers
but cannot write out), throws an `ArgumentError`. Other operations on operator results that contain an
unsized vector, such as `sum(g + p)` or `dot(g, p)`, are kept as calls on the lowered components and are
evaluated once the vector is bound. The exception is a derivative that would act on the unknown elements: a
`Differential` applied to an array expression that contains an unsized vector as a whole, or to an index into
such a composite such as `(g + p)[1]`, throws an `ArgumentError`. A derivative of a scalar reduction such as
`dot(g, p)` or `sum(g + p)` is deferred like any other operation, and so is a bare parameter indexed by an
integer-typed expression, as in `p[1]`, `p[i]` or `p[i + 1]` with an `Int` symbol `i`, which counts as a scalar.
"""
function expand_vector_operators end

@static if isdefined(Symbolics, :Grad)
    const _VectorOperator = Union{Symbolics.Grad, Symbolics.Divergence, Symbolics.Curl}

    _eager(op::Symbolics.Grad) = Symbolics.Grad{true}(op.vars; check = false)
    _eager(op::Symbolics.Divergence) = Symbolics.Divergence{true}(op.vars; check = false)
    _eager(op::Symbolics.Curl{E, N}) where {E, N} = Symbolics.Curl{true, N}(op.vars; check = false)

    _has_vector_operator(ex) = false
    function _has_vector_operator(ex::Union{Num, SymbolicUtils.BasicSymbolic})
        ex = unwrap(ex)
        iscall(ex) || return false
        operation(ex) isa _VectorOperator && return true
        return any(_has_vector_operator, arguments(ex))
    end
    _has_vector_operator(z::Complex) = _has_vector_operator(real(z)) || _has_vector_operator(imag(z))
    _has_vector_operator(eq::Equation) = _has_vector_operator(eq.lhs) || _has_vector_operator(eq.rhs)
    _has_vector_operator(eq::Pair) = _has_vector_operator(first(eq)) || _has_vector_operator(last(eq))
    _has_vector_operator(eqs::Union{AbstractVector, Tuple}) = any(_has_vector_operator, eqs)

    function _check_operator_vars(op, ivs)
        for var in op.vars
            any(iv -> isequal(var, iv), ivs) && continue
            throw(
                ArgumentError(
                    "$(var) in the vector operator $(nameof(typeof(op))) is not an independent variable of the " *
                        "system, whose independent variables are $(join(ivs, ", "))"
                )
            )
        end
        return nothing
    end

    # Applies the eager form of every operator in `ex` to its already-rewritten argument, so the
    # result holds `Differential`s that are not simplified: `∇ ⋅ (k * ∇u)` becomes
    # `Dx(k * Dx(u)) + Dy(k * Dy(u))`. `expand_derivatives` would also do this but it zeroes the
    # derivatives of traces, and expands products that should stay in divergence form.
    function _expand_vector_operators(ex, ivs)
        ex = unwrap(ex)
        iscall(ex) || return ex
        op = operation(ex)
        args = arguments(ex)
        if op isa _VectorOperator
            _check_operator_vars(op, ivs)
            arg = only(args)
            if SymbolicUtils.shape(arg) isa SymbolicUtils.ShapeVecT && length(SymbolicUtils.shape(arg)) == 1
                scalarized = try
                    Symbolics.scalarize(Symbolics.wrap(arg))
                catch err
                    # A sum with an unsized vector has an inferred shape but cannot be written out per component.
                    err isa TypeError || rethrow()
                    throw(
                        ArgumentError(
                            "the operand of the vector operator $(nameof(typeof(op))) contains a vector without a " *
                                "declared length, declare the length of the vector contracted with the operator"
                        )
                    )
                end
                components = [_expand_vector_operators(c, ivs) for c in scalarized]
                return _concrete(_eager(op)(components))
            end
            if op isa Union{Symbolics.Divergence, Symbolics.Curl} && !(SymbolicUtils.shape(arg) isa SymbolicUtils.ShapeVecT)
                # Taking the components named by the coordinates would drop any further component of the bound
                # vector without the dimension error the operator raises for a sized one.
                throw(
                    ArgumentError(
                        "the operand of the vector operator $(nameof(typeof(op))) has no declared length, " *
                            "declare the length of the vector contracted with the operator"
                    )
                )
            end
            return _concrete(_eager(op)(Symbolics.wrap(_expand_vector_operators(arg, ivs))))
        elseif op === getindex && iscall(first(args)) && operation(first(args)) isa _VectorOperator
            components = _expand_vector_operators(first(args), ivs)
            indices = map(unwrap_const, args[2:end])
            all(i -> i isa Union{Integer, Colon}, indices) || throw(
                ArgumentError(
                    "indexing the result of the vector operator $(nameof(typeof(operation(first(args))))) " *
                        "with the symbolic index $(join(indices, ", ")) is not supported, use integer indices"
                )
            )
            return unwrap(components[indices...])
        end
        any(_has_vector_operator, args) || return ex
        newargs = map(a -> _expand_vector_operators(a, ivs), args)
        # An operand of unknown extent anywhere beneath defers every enclosing operation: it stays a call on the
        # lowered components, so nothing above it is scalarized or evaluated before the extent is known and
        # Julia's own length rules apply once the vector is bound.
        if any(_has_unsized, newargs)
            if _is_linalg_dot(op)
                lengths = collect(skipmissing(map(_extent, newargs)))
                allequal(lengths) || throw(ArgumentError("the vectors contracted with `dot` have different lengths $(lengths)"))
            end
            if op isa Differential && (!_scalar_shaped_value(only(newargs)) || _composite_lookup(only(newargs)))
                throw(
                    ArgumentError(
                        "a derivative of an expression that contains a vector without a declared length " *
                            "is not supported, declare the length of the vector"
                    )
                )
            end
            return maketerm(typeof(ex), op, map(unwrap, newargs), metadata(ex))
        end
        # An operator result is a Julia array of expressions, and applying `op` to it (`2g`, `dot(g, g)`,
        # `sum(g)`, `(2g)[1]`) evaluates per component, which a rebuilt symbolic array term would not.
        # A derivative of an operator result is taken per component.
        if op isa Differential && only(newargs) isa AbstractArray
            return map(c -> unwrap(op(c)), only(newargs))
        end
        if _is_linalg_dot(op)
            lengths = collect(skipmissing(map(_extent, newargs)))
            allequal(lengths) || throw(ArgumentError("the vectors contracted with `dot` have different lengths $(lengths)"))
        end
        if any(a -> a isa AbstractArray, newargs)
            # `broadcast` carries its function (`^`, `sin`) as a constant argument, which must be unwrapped to be called.
            callable = try
                map(newargs) do a
                    a isa SymbolicUtils.BasicSymbolic || return a
                    a = unwrap_const(a)
                    a isa SymbolicUtils.BasicSymbolic || return a
                    _scalar_shaped(a) && return Symbolics.wrap(a)
                    known = SymbolicUtils.shape(a) isa SymbolicUtils.ShapeVecT
                    return known ? Symbolics.scalarize(Symbolics.wrap(a)) : Symbolics.wrap(a)
                end
            catch err
                # A sized operand built from an unsized parameter cannot be written out before binding.
                err isa TypeError || rethrow()
                return maketerm(typeof(ex), op, map(unwrap, newargs), metadata(ex))
            end
            result = op(callable...)
            return _unwrap_result(result isa Symbolics.Arr ? Symbolics.scalarize(result) : result)
        end
        rebuilt = maketerm(typeof(ex), op, newargs, metadata(ex))
        # An eager operator can also return a symbolic array, whose operations are only evaluated by `scalarize`.
        arrays = filter(a -> a isa SymbolicUtils.BasicSymbolic && !_scalar_shaped(a), newargs)
        # An array of unknown shape cannot be scalarized and stays as an indexed term.
        (isempty(arrays) || any(a -> !(SymbolicUtils.shape(a) isa SymbolicUtils.ShapeVecT), arrays)) && return rebuilt
        return _unwrap_result(Symbolics.scalarize(Symbolics.wrap(rebuilt)))
    end
    # Whether an unsized symbolic array occurs in `a` or anywhere beneath it.
    _has_unsized(a::AbstractArray) = any(_has_unsized, a)
    _has_unsized(a::Num) = _has_unsized(unwrap(a))
    _has_unsized(a::Symbolics.Arr) = _has_unsized(unwrap(a))
    function _has_unsized(a::SymbolicUtils.BasicSymbolic)
        SymbolicUtils.shape(a) isa SymbolicUtils.ShapeVecT || return true
        # An integer-valued index (a literal or an `Int` symbol) of a bare unsized array is a scalar, so the array
        # does not flow on as a whole. An index into a composite expression is not exempt.
        iscall(a) && operation(a) === getindex && _scalar_shaped(a) && !iscall(first(arguments(a))) &&
            all(_integer_index, arguments(a)[2:end]) && return false
        return iscall(a) && any(_has_unsized, arguments(a))
    end
    _has_unsized(a) = false
    function _integer_index(i)
        i = unwrap_const(i)
        # An integer-typed index expression, which must not itself contain an unsized array.
        return i isa Integer || (i isa SymbolicUtils.BasicSymbolic && SymbolicUtils.symtype(i) <: Integer && !_has_unsized(i))
    end
    # Whether `a` contains a lookup into a composite expression that holds an unsized array as a whole.
    _composite_lookup(a::AbstractArray) = any(_composite_lookup, a)
    _composite_lookup(a::Union{Num, Symbolics.Arr}) = _composite_lookup(unwrap(a))
    function _composite_lookup(a::SymbolicUtils.BasicSymbolic)
        iscall(a) || return false
        if operation(a) === getindex
            indexed = first(arguments(a))
            iscall(indexed) && _has_unsized(indexed) && return true
        end
        return any(_composite_lookup, arguments(a))
    end
    _composite_lookup(a) = false
    _scalar_shaped_value(a) = a isa SymbolicUtils.BasicSymbolic ? _scalar_shaped(a) : !(a isa AbstractArray)
    # The length of a one-dimensional operand when it is known, `missing` otherwise.
    _extent(a::AbstractVector) = a isa Symbolics.Arr ? _extent(unwrap(a)) : length(a)
    function _extent(a::SymbolicUtils.BasicSymbolic)
        shape = SymbolicUtils.shape(a)
        return shape isa SymbolicUtils.ShapeVecT && length(shape) == 1 ? length(only(shape)) : missing
    end
    _extent(a) = missing
    # LinearAlgebra's `dot` by identity, so a user's own function called `dot` is an ordinary function.
    _is_linalg_dot(op) = op === LinearAlgebra.dot
    _scalar_shaped(a) = SymbolicUtils.shape(a) isa SymbolicUtils.ShapeVecT && isempty(SymbolicUtils.shape(a))
    # A vector-valued operator result becomes a Julia array of expressions, so that products with constant
    # matrices, adjoints and contractions (`g' * A * g`, `dot(g, A * g)`) evaluate per component.
    function _concrete(r::Symbolics.Arr)
        SymbolicUtils.shape(unwrap(r)) isa SymbolicUtils.ShapeVecT || return unwrap(r)
        return Symbolics.scalarize(r)
    end
    _concrete(r) = _unwrap_result(r)
    _unwrap_result(r::AbstractArray) = r
    _unwrap_result(r) = unwrap(r)

    # Returns unwrapped sides: wrapping a complex side makes `~` split it into two real equations,
    # which is `handle_complex`'s job.
    _expand_side(side::Complex, ivs) = complex(
        Symbolics.wrap(_expand_side(real(side), ivs)), Symbolics.wrap(_expand_side(imag(side), ivs))
    )
    function _expand_side(side, ivs)
        side = unwrap(side)
        side isa SymbolicUtils.BasicSymbolic && _scalar_shaped(side) &&
            return _expand_vector_operators(side, ivs)
        scalarized = Symbolics.scalarize(Symbolics.wrap(side))
        scalarized isa AbstractArray && return map(c -> _expand_side(c, ivs), scalarized)
        return _expand_vector_operators(scalarized, ivs)
    end
    # An array-valued equation becomes one scalar equation per component, which is all the rest of the
    # pipeline handles.
    function _expand_condition(eq::Equation, ivs)
        lhs, rhs = _expand_side(eq.lhs, ivs), _expand_side(eq.rhs, ivs)
        (lhs isa AbstractArray || rhs isa AbstractArray) || return Equation(lhs, rhs)
        if lhs isa AbstractArray && rhs isa AbstractArray && size(lhs) != size(rhs)
            throw(
                ArgumentError(
                    "the sides of the equation $(eq) have different sizes $(size(lhs)) and $(size(rhs)) " *
                        "once the vector operators are expanded"
                )
            )
        end
        n = length(lhs isa AbstractArray ? lhs : rhs)
        components(side) = side isa AbstractArray ? vec(side) : fill(side, n)
        return Equation[Equation(l, r) for (l, r) in zip(components(lhs), components(rhs))]
    end
    _expand_condition(eq::Pair, ivs) = _expand_side(first(eq), ivs) => _expand_side(last(eq), ivs)

    # Nested groups of equations keep their nesting. The rewrite does not flatten them, and the later complex
    # splitting does not handle nested groups.
    function _expand_conditions(conditions, ivs)
        expanded = Any[]
        for condition in conditions
            if !_has_vector_operator(condition)
                push!(expanded, condition)
            elseif condition isa Union{AbstractVector, Tuple}
                push!(expanded, _expand_conditions(condition, ivs))
            else
                result = _expand_condition(condition, ivs)
                result isa AbstractVector ? append!(expanded, result) : push!(expanded, result)
            end
        end
        return conditions isa Tuple ? Tuple(expanded) : map(identity, expanded)
    end

    function expand_vector_operators(pdesys::PDESystem)
        eqs, bcs = get_eqs(pdesys), get_bcs(pdesys)
        (_has_vector_operator(eqs) || _has_vector_operator(bcs)) || return pdesys
        ivs = unwrap.(get_ivs(pdesys))
        return PDESystem(
            _expand_conditions(eqs, ivs), _expand_conditions(bcs, ivs), get_domain(pdesys), get_ivs(pdesys),
            get_dvs(pdesys), get_ps(pdesys), initial_conditions = pdesys.initial_conditions,
            systems = get_systems(pdesys), connector_type = get_connector_type(pdesys),
            metadata = get_metadata(pdesys), analytic = getfield(pdesys, :analytic),
            analytic_func = getfield(pdesys, :analytic_func),
            gui_metadata = get_gui_metadata(pdesys), description = getfield(pdesys, :description),
            name = getfield(pdesys, :name)
        )
    end
else
    expand_vector_operators(pdesys::PDESystem) = pdesys
end
