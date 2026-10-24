"""
    meshdim(mesh)

Return the spatial dimension of the coordinates of `mesh`. For a surface mesh embedded
in three dimensions this is the embedding dimension, 3.

Part of the mesh-domain interface: PDEBase defines no methods, mesh backends add them.
A `PDESystem` pairing `(x, y) ∈ MeshDomain(mesh)` must name exactly `meshdim(mesh)`
coordinates.
"""
function meshdim end

"""
    facetset_names(mesh)

Return the `Symbol`s naming the boundary facet sets of `mesh`, the named pieces of the
mesh boundary.

Part of the mesh-domain interface: PDEBase defines no methods, mesh backends add them.
"""
function facetset_names end

"""
    cellset_names(mesh)

Return the `Symbol`s naming the cell sets of `mesh`, the named regions of the mesh.

Part of the mesh-domain interface: PDEBase defines no methods, mesh backends add them.
"""
function cellset_names end

"""
    contains_point(mesh, x)

Return whether the point `x` lies in `mesh`. Used by `x ∈ MeshDomain(mesh)`.

Part of the mesh-domain interface. Implementing it is optional; the fallback throws an
error naming the mesh type.
"""
function contains_point(mesh, x)
    return error(
        "PDEBase.contains_point is not implemented for meshes of type $(typeof(mesh)); " *
            "the mesh backend must add a method to test point membership"
    )
end

"""
    backend_mesh(mesh)

Return the native mesh object of the backend that `mesh` wraps. The fallback returns
`mesh` itself.

Part of the mesh-domain interface.
"""
backend_mesh(mesh) = mesh

"""
    MeshDomain(mesh)

A `DomainSets.Domain` whose points are the points of the backend mesh `mesh`, with
element type `SVector{meshdim(mesh), Float64}`.

A `MeshDomain` carries no coordinate names. They come only from the `PDESystem` domain
pairing `(x, y) ∈ Ω`, where tuple position `i` names coordinate component `i` of the mesh;
`VariableMap` checks that the tuple has `meshdim(mesh)` entries.

The mesh type must implement the mesh-domain interface: [`PDEBase.meshdim`](@ref),
[`PDEBase.facetset_names`](@ref), [`PDEBase.cellset_names`](@ref), and optionally
[`PDEBase.contains_point`](@ref) and [`PDEBase.backend_mesh`](@ref).

# Examples
```julia
Ω = MeshDomain(grid)
domains = [(x, y) ∈ Ω, t ∈ Interval(0.0, 1.0)]
```
"""
struct MeshDomain{M, T <: AbstractVector} <: DomainSets.Domain{T}
    mesh::M
end

# The element type of a static unit cube is `SVector{N, Float64}`; reading it from there
# avoids depending on StaticArrays, whose `SVector` DomainSets does not make public.
MeshDomain(mesh) = MeshDomain{typeof(mesh), eltype(DomainSets.UnitCube(Val(meshdim(mesh))))}(mesh)

DomainSets.dimension(d::MeshDomain) = meshdim(d.mesh)
DomainSets.indomain(x, d::MeshDomain) = contains_point(d.mesh, x)

"""
    Region(d::MeshDomain, name::Symbol)

The cell set `name` of the mesh domain `d`, as a `DomainSets.Domain` with the element
type of `d`. Fields: `parent` (the mesh domain) and `name`.

Throws an `ArgumentError` unless `name ∈ cellset_names(d.mesh)`.
"""
struct Region{T, P} <: DomainSets.Domain{T}
    parent::P
    name::Symbol
end

function Region(d::MeshDomain{M, T}, name::Symbol) where {M, T}
    names = cellset_names(d.mesh)
    name in names ||
        throw(ArgumentError("mesh of type $(typeof(d.mesh)) has no cell set $name; its cell sets are $(collect(names))"))
    return Region{T, typeof(d)}(d, name)
end

"""
    regions(d::MeshDomain)

Return a [`Region`](@ref) for every cell set of the mesh of `d`, in the order of
`cellset_names(d.mesh)`.
"""
regions(d::MeshDomain) = [Region(d, name) for name in cellset_names(d.mesh)]

# Qualified so that these stay methods of the upstream functions once DomainSets defines them.
function PDEBase.boundarypieces(d::MeshDomain)
    return [BoundaryPiece(d, name, nothing) for name in facetset_names(d.mesh)]
end

function PDEBase.boundarypiece(d::MeshDomain, name::Symbol)
    names = facetset_names(d.mesh)
    name in names || throw(
        ArgumentError(
            "mesh of type $(typeof(d.mesh)) has no facet set $(repr(name)); " *
                "valid names are $(join(map(repr, collect(names)), ", "))"
        )
    )
    return BoundaryPiece(d, name, nothing)
end

function DomainSets.boundary(d::MeshDomain)
    pieces = boundarypieces(d)
    # DomainSets rejects a `UnionDomain` of no components.
    isempty(pieces) && return DomainSets.EmptySpace{eltype(d)}()
    return DomainSets.UnionDomain(pieces)
end

"""
    BCKind

Supertype of the syntactic boundary-condition kinds returned by [`PDEBase.bckind`](@ref):
[`PDEBase.Essential`](@ref), [`PDEBase.Natural`](@ref) and [`PDEBase.Robin`](@ref).
"""
abstract type BCKind end

"""
    Essential()

The [`PDEBase.BCKind`](@ref) of a condition on trace values only, such as `u(t, Γ) ~ g`.
"""
struct Essential <: BCKind end

"""
    Natural()

The [`PDEBase.BCKind`](@ref) of a condition on derivatives of traces only, such as
`Dn_Γ(u(t, Γ)) ~ h`.
"""
struct Natural <: BCKind end

"""
    Robin()

The [`PDEBase.BCKind`](@ref) of a condition on both trace values and derivatives of
traces, such as `k * Dn_Γ(u(t, Γ)) + α * u(t, Γ) ~ β`.
"""
struct Robin <: BCKind end

"""
    bckind(eq::Equation)

Classify the boundary condition `eq` by its syntax. A trace is a call with a
[`BoundaryPiece`](@ref) argument, such as `u(t, Γ)`. A condition whose traces all sit
inside a `Differential` or [`NormalDerivative`](@ref) is [`PDEBase.Natural`](@ref), one
with no trace inside a derivative is [`PDEBase.Essential`](@ref), and one with both is
[`PDEBase.Robin`](@ref).

The classification is advisory: it does not check coefficients or which pieces the
operators act on. Throws an `ArgumentError` when `eq` contains no trace.
"""
function bckind(eq::Equation)
    under_derivative = Bool[]
    _collect_traces!(under_derivative, eq.lhs, false)
    _collect_traces!(under_derivative, eq.rhs, false)
    isempty(under_derivative) &&
        throw(ArgumentError("$(eq) contains no trace u(..., Γ) on a boundary piece"))
    all(under_derivative) && return Natural()
    any(under_derivative) && return Robin()
    return Essential()
end

function _collect_traces!(under_derivative, ex, in_derivative)
    ex = unwrap_const(safe_unwrap(ex))
    iscall(ex) || return under_derivative
    args = arguments(ex)
    any(a -> unwrap_const(a) isa BoundaryPiece, args) &&
        return push!(under_derivative, in_derivative)
    op = operation(ex)
    in_derivative |= op isa Differential || op isa NormalDerivative
    foreach(a -> _collect_traces!(under_derivative, a, in_derivative), args)
    return under_derivative
end
