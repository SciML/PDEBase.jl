# [Domains](@id domains)

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

## Design notes

The mesh-domain types build on the proposal in
[PDEBase.jl#61](https://github.com/SciML/PDEBase.jl/pull/61), whose domain wrapper,
named boundary regions and discretization hierarchy shaped this design; its two commits
are the base of this branch. The design here differs in three ways.

- Named facet sets and cell sets are read from the mesh through the interface
  ([`PDEBase.facetset_names`](@ref), [`PDEBase.cellset_names`](@ref)), not stored as
  markers on the domain. A mesh is the single source of truth for its sets, so a domain
  cannot disagree with its mesh, and a backend that already builds sets from predicates
  (for example Ferrite's `addfacetset!`) needs nothing more.
- A named piece of the boundary is a [`BoundaryPiece`](@ref), a `DomainSets.Domain` that
  can appear as an argument of a trace `u(t, Γ)`, so boundary conditions are written in
  the same syntax on a box and on a mesh.
- PDEBase makes no assumption about the discretization method. It does not assume a
  weak form, a set of fields, a quadrature rule for each part of a weak form, an ansatz
  for each field, or that the problem is linear. Those choices belong to the discretizer,
  which reads the `PDESystem` and the mesh through the interface above.
