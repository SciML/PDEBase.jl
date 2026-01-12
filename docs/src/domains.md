# Domain Types for Unstructured Meshes

PDEBase.jl provides abstract types and interfaces for working with discretized domains (meshes) from various FEM/FVM/DG frameworks.

## Overview

While `DomainSets.jl` provides analytical domain descriptions (intervals, products, set operations), many PDE solvers work with numerical meshes that may not have analytical descriptions. PDEBase provides types to bridge this gap.

## Type Hierarchy

```
DomainSets.Domain{T}
    AbstractDiscretizedDomain{T}
        MeshDomain{T,M}
```

## Abstract Types

### AbstractDiscretizedDomain

```@docs
PDEBase.AbstractDiscretizedDomain
```

## Concrete Types

### MeshDomain

```@docs
PDEBase.MeshDomain
```

The `MeshDomain` type wraps mesh objects from external libraries and provides:
- Storage for the underlying mesh
- Named boundary region markers
- A common interface for PDEBase operations

### Example Usage

```julia
using PDEBase

# Assume `grid` is a mesh object from an external library (Ferrite, Gridap, etc.)
# Create a MeshDomain wrapper
domain = MeshDomain(grid)

# Mark boundary regions
mark_boundary!(domain, "inlet", inlet_facets)
mark_boundary!(domain, "outlet", outlet_facets)
mark_boundary!(domain, "wall", wall_facets)

# Access boundary markers
markers = get_boundary_markers(domain)
# markers["inlet"] => inlet_facets

# Access the underlying mesh
mesh = get_mesh(domain)
```

## Boundary Regions

PDEBase provides types for referencing named boundary regions on domains.

### BoundaryRegion

```@docs
PDEBase.BoundaryRegion
```

Use `BoundaryRegion` to reference named boundaries when specifying boundary conditions:

```julia
using PDEBase

domain = MeshDomain(grid)
mark_boundary!(domain, "left", left_faces)
mark_boundary!(domain, "right", right_faces)

# Create boundary region references
left_boundary = BoundaryRegion(domain, "left")
right_boundary = BoundaryRegion(domain, "right")

# These can be used in boundary condition specifications
# (syntax depends on downstream discretization package)
```

### ConditionalBoundary

```@docs
PDEBase.ConditionalBoundary
```

Use `ConditionalBoundary` when boundaries are defined by conditions rather than explicit markers:

```julia
using PDEBase, Symbolics
@parameters x y

# For analytical domains, select boundaries by condition
domain = (x in Interval(0, 1)) * (y in Interval(0, 1))
left_bc = ConditionalBoundary(domain, [x ~ 0])
right_bc = ConditionalBoundary(domain, [x ~ 1])
```

## Interface Functions

### Mesh Access

```@docs
PDEBase.get_mesh
```

### Boundary Markers

```@docs
PDEBase.get_boundary_markers
PDEBase.mark_boundary!
```

### Boundary Region Access

```@docs
PDEBase.get_parent_domain
PDEBase.get_boundary_name
PDEBase.get_condition
```

## Integration with External Mesh Libraries

### Ferrite.jl

```julia
using Ferrite, PDEBase

# Create a Ferrite grid
grid = generate_grid(Quadrilateral, (10, 10))

# Wrap in MeshDomain
domain = MeshDomain(grid)

# Use Ferrite's face sets as boundary markers
mark_boundary!(domain, "left", getfacetset(grid, "left"))
mark_boundary!(domain, "right", getfacetset(grid, "right"))
```

### Gridap.jl

```julia
using Gridap, PDEBase

# Create a Gridap model
model = CartesianDiscreteModel((0,1,0,1), (10,10))

# Wrap in MeshDomain
domain = MeshDomain(model)

# Use boundary tags
mark_boundary!(domain, "boundary", "boundary")
```

### GMSH Import

```julia
using Gmsh, PDEBase

# Load mesh from GMSH
gmsh.initialize()
gmsh.open("domain.msh")
# ... extract mesh data ...

# Wrap mesh with boundary markers from GMSH physical groups
domain = MeshDomain(mesh_data)
for (name, tag) in gmsh_physical_groups
    mark_boundary!(domain, name, tag)
end
```

## Design Notes

The domain types in PDEBase are intentionally minimal and abstract. The philosophy is:

1. **Framework-agnostic**: PDEBase doesn't depend on any specific mesh library
2. **Type wrapper**: `MeshDomain` stores mesh objects opaquely; framework-specific packages interpret them
3. **Named boundaries**: Support for named boundary regions enables clean BC specification
4. **Extensible**: Downstream packages can subtype `AbstractDiscretizedDomain` for specialized needs

For actual mesh operations (element iteration, DOF mapping, assembly), users should use the native APIs of their mesh library, accessed via `get_mesh()`.
