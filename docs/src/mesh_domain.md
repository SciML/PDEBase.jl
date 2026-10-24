# [Mesh Domains](@id mesh_domain)

A `MeshDomain` lets a `PDESystem` be posed on an unstructured mesh. PDEBase defines
the type and a five-function interface; mesh backends (for example a Ferrite adapter)
add the methods for their mesh type.

The mesh carries no coordinate names. The domain pairing `(x, y) ∈ Ω` binds them:
tuple position `i` names coordinate component `i` of the mesh, and `VariableMap`
checks that the tuple has `meshdim(mesh)` entries. `spatial_domain(v)` returns the
mesh domain (or, for analytic problems, the domain of the spatial coordinates).

```julia
Ω = MeshDomain(grid)
domains = [(x, y) ∈ Ω, t ∈ Interval(0.0, 1.0)]
```

A tuple pairing on a `Rectangle` or another product of intervals also fills
`v.intervals` per axis, so box discretizers accept `(x, y) ∈ Rectangle(...)`.

## Types

```@docs
PDEBase.MeshDomain
PDEBase.Region
PDEBase.regions
PDEBase.spatial_domain
```

## Backend interface

```@docs
PDEBase.meshdim
PDEBase.facetset_names
PDEBase.cellset_names
PDEBase.contains_point
PDEBase.backend_mesh
```

## Boundary pieces and condition kinds

`boundarypiece(Ω, name)` returns the facet set `name` of the mesh as a `BoundaryPiece`
with no geometry, `boundarypieces(Ω)` returns one per facet set, and
`DomainSets.boundary(Ω)` is their `UnionDomain`. `PDEBase.bckind` classifies a
condition on traces by its syntax.

```@docs
PDEBase.BCKind
PDEBase.Essential
PDEBase.Natural
PDEBase.Robin
PDEBase.bckind
```
