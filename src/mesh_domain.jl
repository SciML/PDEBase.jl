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
