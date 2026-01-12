############################################################################################
# Abstract types for variational/weak-form PDE discretization methods
############################################################################################

"""
    AbstractVariationalMethod

Abstract supertype for variational formulation methods used in PDE discretization.

Variational methods reformulate PDEs in weak form, seeking solutions that satisfy
integral equations rather than pointwise equations. This includes finite element
methods (FEM), spectral methods, and other Galerkin-type approaches.

# Subtypes
- [`AbstractGalerkinMethod`](@ref): Methods where solution is projected onto a finite basis

# See Also
- [`AbstractRitzGalerkin`](@ref): Standard Galerkin with same test and trial spaces
- [`AbstractPetrovGalerkin`](@ref): Galerkin with different test and trial spaces
"""
abstract type AbstractVariationalMethod end

"""
    AbstractGalerkinMethod <: AbstractVariationalMethod

Abstract type for Galerkin-type methods that approximate solutions by projection
onto finite-dimensional subspaces.

The general Galerkin approach:
1. Choose finite-dimensional trial space Vₕ for the solution uₕ
2. Choose finite-dimensional test space Wₕ for the test functions vₕ
3. Find uₕ ∈ Vₕ such that B(uₕ, vₕ) = L(vₕ) for all vₕ ∈ Wₕ

where B(·,·) is a bilinear (or semilinear) form and L(·) is a linear form.

# Subtypes
- [`AbstractRitzGalerkin`](@ref): Vₕ = Wₕ (standard Galerkin)
- [`AbstractPetrovGalerkin`](@ref): Vₕ ≠ Wₕ (e.g., SUPG, DG)
"""
abstract type AbstractGalerkinMethod <: AbstractVariationalMethod end

"""
    AbstractRitzGalerkin <: AbstractGalerkinMethod

Abstract type for Ritz-Galerkin methods where the test and trial spaces are identical.

Standard finite element methods use Ritz-Galerkin formulations. The weak form is:
Find uₕ ∈ Vₕ such that a(uₕ, vₕ) = b(vₕ) for all vₕ ∈ Vₕ

where:
- a(·,·) is the bilinear form from the PDE
- b(·) is the linear form (source terms, boundary terms)
- Vₕ is the finite element space

# Example Implementation
```julia
struct ContinuousGalerkin <: AbstractRitzGalerkin
    polynomial_order::Int
    quadrature_order::Int
end
```
"""
abstract type AbstractRitzGalerkin <: AbstractGalerkinMethod end

"""
    AbstractPetrovGalerkin <: AbstractGalerkinMethod

Abstract type for Petrov-Galerkin methods where test and trial spaces differ.

Petrov-Galerkin methods are useful for:
- Stabilization (SUPG for advection-dominated problems)
- Discontinuous Galerkin (DG) methods
- Streamline-upwind methods

The weak form is:
Find uₕ ∈ Vₕ such that a(uₕ, wₕ) = b(wₕ) for all wₕ ∈ Wₕ

where Vₕ (trial space) ≠ Wₕ (test space).
"""
abstract type AbstractPetrovGalerkin <: AbstractGalerkinMethod end

"""
    RitzGalerkin

A basic marker type for standard Ritz-Galerkin finite element discretization.

This type can be used as a tag to indicate that a discretization should use
standard Galerkin formulation. Framework-specific discretizers can dispatch
on this type to apply appropriate transformations.

# Usage
```julia
# In a downstream package
struct MyFEMDiscretization{M} <: AbstractEquationSystemDiscretization
    method::M
    # ... other fields
end

# Create a standard Galerkin discretization
disc = MyFEMDiscretization(RitzGalerkin())
```
"""
struct RitzGalerkin <: AbstractRitzGalerkin end

############################################################################################
# Abstract types for unstructured mesh discretizations
############################################################################################

"""
    AbstractUnstructuredDiscretization <: AbstractEquationSystemDiscretization

Abstract type for discretizations on unstructured meshes.

Unstructured discretizations work with meshes that don't have a regular Cartesian
structure. This includes triangular/tetrahedral meshes, quadrilateral/hexahedral
meshes, and general polyhedral meshes.

Discretizations subtyping this should implement the standard PDEBase interface
functions, with the discrete space being an [`AbstractUnstructuredDiscreteSpace`](@ref).

# Example
```julia
# In a downstream package (e.g., for Ferrite.jl or Gridap.jl integration)
struct FerriteDiscretization{M, Q} <: AbstractUnstructuredDiscretization
    method::M       # Variational method (e.g., RitzGalerkin)
    quadrature::Q   # Quadrature specification
    # ... other fields
end
```
"""
abstract type AbstractUnstructuredDiscretization <: AbstractEquationSystemDiscretization end

"""
    AbstractUnstructuredDiscreteSpace <: AbstractDiscreteSpace

Abstract type for discrete spaces on unstructured meshes.

This represents the discretized spatial domain and associated finite element
spaces for unstructured mesh discretizations. Unlike Cartesian discrete spaces
which have regular grid structure, unstructured spaces may have:
- Variable element shapes (triangles, quads, tetrahedra, hexahedra, etc.)
- Non-uniform element sizes
- Arbitrary connectivity

# Interface
Implementations should provide:
- Storage for discrete variables (DOF values)
- Mapping between local (element) and global DOF numbering
- Access to the underlying mesh geometry

# See Also
- `AbstractCartesianDiscreteSpace`: For structured Cartesian grids
"""
abstract type AbstractUnstructuredDiscreteSpace <: AbstractDiscreteSpace end

############################################################################################
# Interface function stubs for variational methods
############################################################################################

"""
    get_variational_method(disc::AbstractUnstructuredDiscretization)

Get the variational method (e.g., `RitzGalerkin`) used by a discretization.

Default implementation returns `nothing`. Discretizations using variational
methods should override this.
"""
get_variational_method(::AbstractUnstructuredDiscretization) = nothing

"""
    is_galerkin(method::AbstractVariationalMethod) -> Bool

Check if a variational method is a Galerkin-type method.
"""
is_galerkin(::AbstractGalerkinMethod) = true
is_galerkin(::AbstractVariationalMethod) = false

"""
    is_symmetric_galerkin(method::AbstractVariationalMethod) -> Bool

Check if a Galerkin method uses symmetric test/trial spaces (Ritz-Galerkin).
"""
is_symmetric_galerkin(::AbstractRitzGalerkin) = true
is_symmetric_galerkin(::AbstractVariationalMethod) = false

############################################################################################
# Traits for FEM discretization capabilities
############################################################################################

"""
    abstract type AbstractFEMTrait end

Base type for traits describing FEM discretization capabilities.
"""
abstract type AbstractFEMTrait end

"""
    SupportsWeakForm <: AbstractFEMTrait

Trait indicating a discretization can handle PDEs in weak form.
"""
struct SupportsWeakForm <: AbstractFEMTrait end

"""
    SupportsStrongForm <: AbstractFEMTrait

Trait indicating a discretization can handle PDEs in strong form.
"""
struct SupportsStrongForm <: AbstractFEMTrait end

"""
    SupportsBothForms <: AbstractFEMTrait

Trait indicating a discretization can handle both weak and strong form PDEs.
"""
struct SupportsBothForms <: AbstractFEMTrait end

"""
    get_supported_form(disc::AbstractDiscretization) -> AbstractFEMTrait

Query what form (weak/strong/both) a discretization supports.

Default returns `SupportsStrongForm()`. FEM discretizations should override
to return `SupportsWeakForm()` or `SupportsBothForms()`.
"""
get_supported_form(::AbstractDiscretization) = SupportsStrongForm()
