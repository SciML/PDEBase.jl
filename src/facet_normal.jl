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
    return any(a -> _has_facet_normal(a, test), arguments(ex))
end
