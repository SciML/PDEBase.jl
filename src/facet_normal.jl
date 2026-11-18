### FACET NORMAL ###

"""
    FacetNormal(Γ::BoundaryPiece)

The outward unit normal of the boundary piece `Γ`, as an operator that gives its components.
`FacetNormal(Γ)(d)` is the symbolic term for component `d`, printed `n_<name>(d)`, where `d`
is a position in `indvars(v)`, the spatial independent variables of the system in their
declared order, as for the gradient slots of the strong-form analysis. Use
[`normal`](@ref PDEBase.normal)`(Γ)` to get all components as a vector.

A component is an opaque real scalar for Symbolics: `isequal` and `hash` compare the pieces,
it carries no derivative and counts as order 0. Do not differentiate it: Symbolics evaluates
`Differential(x)(n_Γ(1))` to `0`, so PDEBase rejects it in boundary conditions.

In a boundary condition on `Γ` each component becomes the normal slot `N[d]` of the
strong-form analysis, or the numeric `±1` or `0` when `Γ` is a face of a box. A Cartesian
derivative of a trace on a piece, such as `Dx(E(t, Γ))`, is admitted only in a condition that
also contains a facet normal of that piece.

# Fields
- `piece`: The boundary piece whose outward normal this is.
"""
struct FacetNormal <: SymbolicUtils.Operator
    piece::BoundaryPiece
end

function (n::FacetNormal)(d::Integer)
    dim = DomainSets.dimension(parent(n.piece))
    1 <= d <= dim || throw(
        ArgumentError(
            "component $(d) of the normal of $(n.piece), which is a piece of a $(dim)-dimensional domain; " *
                "valid components are 1:$(dim)"
        )
    )
    T = SymbolicUtils.SymReal
    return Num(
        SymbolicUtils.BSImpl.Term{T}(
            n, SymbolicUtils.ArgsT{T}((SymbolicUtils.Const{T}(Int(d)),));
            type = Real, shape = SymbolicUtils.ShapeVecT()
        )
    )
end
SymbolicUtils.isbinop(::FacetNormal) = false

Base.:(==)(a::FacetNormal, b::FacetNormal) = isequal(a, b)
Base.isequal(a::FacetNormal, b::FacetNormal) = isequal(a.piece, b.piece)
Base.hash(a::FacetNormal, h::UInt) = hash(a.piece, hash(FacetNormal, h))
Base.show(io::IO, n::FacetNormal) = print(io, "n_", n.piece.name)

"The component index `d` of the facet-normal term `ex`."
_normal_component(ex) = Int(unwrap_const(only(arguments(unwrap(ex)))))

_is_facet_normal(ex) = iscall(ex) && operation(ex) isa FacetNormal

"Does `ex` hold a facet normal anywhere? With `piece`, one of that piece only."
function _has_facet_normal(ex, test = _ -> true)
    ex = safe_unwrap(ex)
    iscall(ex) || return false
    _is_facet_normal(ex) && test(operation(ex).piece) && return true
    operation(ex) isa Differential && _has_facet_normal(operation(ex).x, test) && return true
    return any(a -> _has_facet_normal(a, test), arguments(ex))
end

"Reject a `Differential` node `ex` whose operand holds a facet normal: the normal is not a field to differentiate."
function _check_normal_not_differentiated(ex)
    _has_facet_normal(ex) && error(
        "derivative of an expression containing the facet normal is not supported in a boundary condition: $(Symbolics.wrap(ex))"
    )
    return nothing
end

"""
Throw if a condition that held a facet normal has no derivative left in its residual, unless the
residual is identically zero there (a zero natural condition, as `x * n[1] * Dx(u) ~ 0` at `x = 0`).
"""
function _check_normal_has_derivative(had_normal, has_derivative, residual_zero, label)
    had_normal && !has_derivative && !residual_zero && _normal_without_derivative_error(label)
    return nothing
end

"The error for a facet normal in a gradient-free (essential) condition `label`."
_normal_without_derivative_error(label) = error(
    "the condition $label contains the facet normal but no derivative of a trace; the normal can appear only in a natural or Robin-type condition"
)

"""
Throw an `ArgumentError` if `ex`, a condition `label` after the numeric facet normal replaced the symbolic one,
holds a constant that is not finite: it divided by a component of the normal that is zero on this face.
"""
function _check_finite_constants(ex, label)
    ex = safe_unwrap(ex)
    c = unwrap_const(ex)
    c isa Number && !isfinite(c) && throw(
        ArgumentError(
            "the condition $label divides by a component of the facet normal that is zero on this face, giving $c"
        )
    )
    iscall(ex) && foreach(a -> _check_finite_constants(a, label), arguments(ex))
    return nothing
end

"The facet-normal terms in `ex`."
function _normal_terms!(terms, ex)
    ex = safe_unwrap(ex)
    iscall(ex) || return terms
    _is_facet_normal(ex) && push!(terms, ex)
    operation(ex) isa Differential && _normal_terms!(terms, operation(ex).x)
    foreach(a -> _normal_terms!(terms, a), arguments(ex))
    return terms
end

# The number `ex` takes when the symbols in `rules` get their values, or `nothing` if it also holds other
# symbols. Evaluated bottom-up on the written tree, so no simplification can cancel a pole first.
function _numeric_value(ex, rules)
    ex = safe_unwrap(ex)
    haskey(rules, ex) && return rules[ex]
    c = unwrap_const(ex)
    c isa Number && return c
    iscall(ex) || return nothing
    values = map(a -> _numeric_value(a, rules), arguments(ex))
    any(isnothing, values) && return nothing
    return try
        operation(ex)(values...)
    catch e
        e isa Union{DomainError, DivideError, ArgumentError} || rethrow()
        NaN
    end
end

"""
Throw an `ArgumentError` if a written side of the condition `label` divides by zero, or is not finite, once
`rules` (the numeric facet normal and the face coordinate) give values to its symbols. Each side is checked on
its own and as written, so a pole cannot cancel against the other side or against itself (`x / (x + n[2])`).
"""
function _check_sides_finite(sides, rules, label)
    undefined(why) = throw(ArgumentError("the condition $label is undefined on this face: $why"))
    function walk(ex)
        ex = safe_unwrap(ex)
        iscall(ex) || return nothing
        op, args = operation(ex), arguments(ex)
        den = op === (/) ? args[2] :
            op === (^) && (e = unwrap_const(safe_unwrap(args[2]))) isa Real && e < 0 ? args[1] : nothing
        if den !== nothing
            value = _numeric_value(den, rules)
            value isa Number && iszero(value) && undefined("it divides by $(Symbolics.wrap(den)), which is zero there")
        end
        value = _numeric_value(ex, rules)
        value isa Number && !isfinite(value) && undefined("$(Symbolics.wrap(ex)) is $value there")
        return foreach(walk, args)
    end
    foreach(walk, map(unwrap, sides))
    return nothing
end
