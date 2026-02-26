############################################################################################
# Domain types for unstructured meshes and general domains
############################################################################################

"""
    AbstractDiscretizedDomain{T} <: DomainSets.Domain{T}

Abstract supertype for domains that have been discretized into meshes or grids.
This provides an interface layer between analytical domain descriptions (from DomainSets.jl)
and numerical mesh representations from various FEM/FVM/DG frameworks.

Concrete implementations should wrap mesh objects from packages like Ferrite.jl, Gridap.jl,
Trixi.jl, or other meshing libraries.

# Type Parameter
- `T`: The element type of the domain coordinates (e.g., `Float64`, `SVector{3,Float64}`)

# Interface Functions
Implementations should define:
- `DomainSets.dimension(d::AbstractDiscretizedDomain)` - Spatial dimension of the domain
- `Base.in(point, d::AbstractDiscretizedDomain)` - Point membership test
- `get_boundary_markers(d::AbstractDiscretizedDomain)` - Named boundary regions

# See Also
- [`MeshDomain`](@ref) for a generic mesh wrapper
- [`BoundaryRegion`](@ref) for named boundary regions
"""
abstract type AbstractDiscretizedDomain{T} <: DomainSets.Domain{T} end

"""
    MeshDomain{T, M} <: AbstractDiscretizedDomain{T}

A generic wrapper for mesh objects from external meshing libraries.

This type allows PDEBase to work with meshes from any FEM/FVM/DG framework by
wrapping the native mesh object and providing a common interface.

# Type Parameters
- `T`: The element type of domain coordinates
- `M`: The type of the wrapped mesh object

# Fields
- `mesh::M`: The underlying mesh object from an external library
- `boundary_markers::Dict{String, Any}`: Named boundary regions and their markers

# Example
```julia
using PDEBase
# Assuming `grid` is a mesh from some external library
mesh_domain = MeshDomain{Float64}(grid)
mark_boundary!(mesh_domain, "inlet", inlet_facets)
mark_boundary!(mesh_domain, "outlet", outlet_facets)
```

# Note
The mesh object is stored opaquely. Framework-specific discretization packages
should access it via `get_mesh(domain)` and handle it appropriately.
"""
struct MeshDomain{T, M} <: AbstractDiscretizedDomain{T}
    mesh::M
    boundary_markers::Dict{String, Any}
    function MeshDomain{T}(mesh::M) where {T, M}
        new{T, M}(mesh, Dict{String, Any}())
    end
    function MeshDomain{T, M}(mesh::M, markers::Dict{String, Any}) where {T, M}
        new{T, M}(mesh, markers)
    end
end

"""
    MeshDomain(mesh; eltype=Float64)

Construct a MeshDomain from a mesh object with default coordinate type.

# Arguments
- `mesh`: The mesh object from an external library
- `eltype`: The coordinate element type (default: `Float64`)
"""
MeshDomain(mesh; eltype::Type{T} = Float64) where {T} = MeshDomain{T}(mesh)

"""
    get_mesh(domain::MeshDomain)

Retrieve the underlying mesh object from a MeshDomain.
"""
get_mesh(domain::MeshDomain) = domain.mesh

"""
    get_boundary_markers(domain::AbstractDiscretizedDomain)

Get the dictionary of named boundary regions for a discretized domain.
Returns a `Dict{String, Any}` mapping boundary names to their markers/identifiers.

Default implementation returns an empty dictionary.
"""
get_boundary_markers(::AbstractDiscretizedDomain) = Dict{String, Any}()
get_boundary_markers(domain::MeshDomain) = domain.boundary_markers

"""
    mark_boundary!(domain::MeshDomain, name::String, marker)

Add a named boundary region to a MeshDomain.

# Arguments
- `domain`: The MeshDomain to modify
- `name`: A string identifier for the boundary (e.g., "inlet", "wall", "left")
- `marker`: The boundary marker/identifier (framework-specific)

# Returns
The modified domain (for chaining).

# Example
```julia
domain = MeshDomain(grid)
mark_boundary!(domain, "inlet", inlet_faces)
mark_boundary!(domain, "wall", wall_faces)
```
"""
function mark_boundary!(domain::MeshDomain, name::String, marker)
    domain.boundary_markers[name] = marker
    return domain
end

############################################################################################
# Boundary region types for named subdomains
############################################################################################

"""
    BoundaryRegion{D}

Represents a named boundary region on a domain.

This is distinct from boundary *conditions* (like `LowerBoundary`, `UpperBoundary`).
A `BoundaryRegion` identifies a portion of the domain boundary by name, which can
then be used when specifying boundary conditions.

# Type Parameters
- `D`: The type of the parent domain

# Fields
- `domain::D`: The parent domain containing this boundary
- `name::String`: Identifier for this boundary region

# Example
```julia
using PDEBase

@parameters x
Omega = MeshDomain(grid)
mark_boundary!(Omega, "left", left_facets)

# Create a boundary region reference
partial_Omega_left = BoundaryRegion(Omega, "left")

# Use in boundary conditions (syntax supported by downstream packages)
# bcs = [u(partial_Omega_left) ~ 0.0]
```

# Note
For analytical domains (from DomainSets), you can also create boundary regions
using selector conditions. See [`ConditionalBoundary`](@ref).
"""
struct BoundaryRegion{D}
    domain::D
    name::String
end

"""
    BoundaryRegion(domain, name::Symbol)

Convenience constructor accepting a Symbol for the name.
"""
BoundaryRegion(domain::D, name::Symbol) where {D} = BoundaryRegion{D}(domain, String(name))

"""
    get_parent_domain(br::BoundaryRegion)

Get the parent domain of a boundary region.
"""
get_parent_domain(br::BoundaryRegion) = br.domain

"""
    get_boundary_name(br::BoundaryRegion)

Get the name/identifier of a boundary region.
"""
get_boundary_name(br::BoundaryRegion) = br.name

"""
    ConditionalBoundary{D, C}

Represents a boundary region selected by a condition on coordinates.

This allows specifying boundaries without explicit markers, useful for
analytical domains or simple geometries.

# Type Parameters
- `D`: The type of the parent domain
- `C`: The type of the selection condition

# Fields
- `domain::D`: The parent domain
- `condition::C`: A condition identifying the boundary (e.g., equations like `x ~ 0`)

# Example
```julia
using PDEBase, Symbolics
@parameters x y

Omega = (x in Interval(0, 1)) × (y in Interval(0, 1))

# Select boundary where x = 0
left_boundary = ConditionalBoundary(Omega, [x ~ 0])
```
"""
struct ConditionalBoundary{D, C}
    domain::D
    condition::C
end

"""
    get_parent_domain(cb::ConditionalBoundary)

Get the parent domain of a conditional boundary.
"""
get_parent_domain(cb::ConditionalBoundary) = cb.domain

"""
    get_condition(cb::ConditionalBoundary)

Get the selection condition for a conditional boundary.
"""
get_condition(cb::ConditionalBoundary) = cb.condition
