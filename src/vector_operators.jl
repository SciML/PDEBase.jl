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
    _has_vector_operator(eq::Equation) = _has_vector_operator(eq.lhs) || _has_vector_operator(eq.rhs)
    _has_vector_operator(eq::Pair) = _has_vector_operator(first(eq)) || _has_vector_operator(last(eq))
    _has_vector_operator(eqs::AbstractVector) = any(_has_vector_operator, eqs)

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
                components = [_expand_vector_operators(c, ivs) for c in Symbolics.scalarize(Symbolics.wrap(arg))]
                return _unwrap_result(_eager(op)(components))
            end
            return _unwrap_result(_eager(op)(Num(_expand_vector_operators(arg, ivs))))
        elseif op === getindex && iscall(first(args)) && operation(first(args)) isa _VectorOperator
            components = _expand_vector_operators(first(args), ivs)
            return unwrap(components[(Int(unwrap_const(i)) for i in args[2:end])...])
        end
        any(_has_vector_operator, args) || return ex
        return maketerm(typeof(ex), op, map(a -> _expand_vector_operators(a, ivs), args), metadata(ex))
    end
    _unwrap_result(r::AbstractArray) = r
    _unwrap_result(r) = unwrap(r)

    function _expand_side(side, ivs)
        scalarized = Symbolics.scalarize(Symbolics.wrap(side))
        scalarized isa AbstractArray && return map(c -> _expand_side(c, ivs), scalarized)
        return Symbolics.wrap(_expand_vector_operators(scalarized, ivs))
    end
    _expand_condition(eq::Equation, ivs) = _expand_side(eq.lhs, ivs) ~ _expand_side(eq.rhs, ivs)
    _expand_condition(eq::Pair, ivs) = _expand_side(first(eq), ivs) => _expand_side(last(eq), ivs)

    function expand_vector_operators(pdesys::PDESystem)
        eqs, bcs = get_eqs(pdesys), get_bcs(pdesys)
        (_has_vector_operator(eqs) || _has_vector_operator(bcs)) || return pdesys
        ivs = unwrap.(get_ivs(pdesys))
        expand(eq) = _has_vector_operator(eq) ? _expand_condition(eq, ivs) : eq
        return PDESystem(
            map(expand, eqs), map(expand, bcs), get_domain(pdesys), get_ivs(pdesys),
            get_dvs(pdesys), get_ps(pdesys), initial_conditions = pdesys.initial_conditions,
            systems = get_systems(pdesys), connector_type = get_connector_type(pdesys),
            metadata = get_metadata(pdesys), analytic = getfield(pdesys, :analytic),
            analytic_func = getfield(pdesys, :analytic_func),
            gui_metadata = get_gui_metadata(pdesys), name = getfield(pdesys, :name)
        )
    end
else
    expand_vector_operators(pdesys::PDESystem) = pdesys
end
