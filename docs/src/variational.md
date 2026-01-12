# Variational Methods and FEM Interface

PDEBase.jl provides abstract types and interfaces for finite element methods (FEM) and other variational discretization approaches.

## Overview

Variational methods reformulate PDEs from their strong (pointwise) form into a weak (integral) form. This is the foundation of:

- Finite Element Methods (FEM)
- Spectral Galerkin Methods
- Discontinuous Galerkin (DG) Methods
- Finite Volume Methods (FVM) with variational interpretation

PDEBase provides abstract types that discretization packages can subtype to indicate their variational formulation.

## Type Hierarchy

```
AbstractVariationalMethod
    AbstractGalerkinMethod
        AbstractRitzGalerkin
            RitzGalerkin
        AbstractPetrovGalerkin
```

## Variational Method Types

### AbstractVariationalMethod

```@docs
PDEBase.AbstractVariationalMethod
```

### AbstractGalerkinMethod

```@docs
PDEBase.AbstractGalerkinMethod
```

### AbstractRitzGalerkin

```@docs
PDEBase.AbstractRitzGalerkin
```

### AbstractPetrovGalerkin

```@docs
PDEBase.AbstractPetrovGalerkin
```

### RitzGalerkin

```@docs
PDEBase.RitzGalerkin
```

## Unstructured Discretization Types

For discretizations on unstructured meshes (triangular, tetrahedral, general polyhedral):

### AbstractUnstructuredDiscretization

```@docs
PDEBase.AbstractUnstructuredDiscretization
```

### AbstractUnstructuredDiscreteSpace

```@docs
PDEBase.AbstractUnstructuredDiscreteSpace
```

## Query Functions

### get_variational_method

```@docs
PDEBase.get_variational_method
```

### is_galerkin

```@docs
PDEBase.is_galerkin
```

### is_symmetric_galerkin

```@docs
PDEBase.is_symmetric_galerkin
```

## FEM Traits

PDEBase provides traits to query discretization capabilities:

```@docs
PDEBase.AbstractFEMTrait
PDEBase.SupportsWeakForm
PDEBase.SupportsStrongForm
PDEBase.SupportsBothForms
PDEBase.get_supported_form
```

## Example: Implementing a FEM Discretizer

Here's how a downstream package might use these types:

```julia
using PDEBase, ModelingToolkit, SciMLBase

# Define a FEM discretization using Ritz-Galerkin
struct MyFEMDiscretization{M<:AbstractVariationalMethod, Q} <: AbstractUnstructuredDiscretization
    method::M           # Variational method (e.g., RitzGalerkin())
    polynomial_order::Int
    quadrature::Q
    # ... other fields
end

# Default constructor with standard Galerkin
function MyFEMDiscretization(; order=1, quadrature=GaussQuadrature(order+1))
    MyFEMDiscretization(RitzGalerkin(), order, quadrature)
end

# Implement the method query
PDEBase.get_variational_method(d::MyFEMDiscretization) = d.method

# Indicate we support weak form
PDEBase.get_supported_form(::MyFEMDiscretization) = SupportsWeakForm()

# Define discrete space for unstructured meshes
struct MyFEMSpace{M} <: AbstractUnstructuredDiscreteSpace
    mesh::M
    dof_handler::Any
    # ... DOF mapping, element info, etc.
end

# Implement PDEBase interface functions
function PDEBase.construct_discrete_space(v::VariableMap, disc::MyFEMDiscretization)
    # Build finite element space from variable map
    # ...
    return MyFEMSpace(mesh, dof_handler)
end

function PDEBase.discretize_equation!(
    disc_state,
    pde::Equation,
    vareqmap,
    eqvar,
    bcmap,
    depvars,
    s::MyFEMSpace,
    derivweights,
    indexmap,
    disc::MyFEMDiscretization
)
    # Apply weak form transformation
    # Assemble element matrices/vectors
    # Handle boundary conditions
    # ...
end
```

## Strong Form to Weak Form Transformation

For discretizations that need to transform PDEs from strong to weak form, implement `should_transform` and `transform_pde_system!`:

```julia
using PDEBase

# Indicate transformation is needed for Galerkin methods
function PDEBase.should_transform(pdesys::PDESystem, disc::MyFEMDiscretization, boundarymap)
    return is_galerkin(disc.method)
end

# Perform the transformation
function PDEBase.transform_pde_system!(
    v::VariableMap,
    boundarymap,
    pdesys::PDESystem,
    disc::MyFEMDiscretization
)
    # 1. Introduce test functions v for each trial function u
    # 2. Multiply equations by test functions
    # 3. Integrate over domain
    # 4. Apply integration by parts where needed
    # 5. Return transformed PDESystem
    # ...
end
```

## Weak Form Representation

When working with weak forms, the PDE is expressed as:

**Find u in V such that:** `a(u, v) = L(v)` for all v in V

Where:
- `a(u, v)` is the bilinear form (depends on both u and v)
- `L(v)` is the linear form (depends only on v)
- `V` is the function space (trial = test for Ritz-Galerkin)

Example for the heat equation `-div(grad(u)) = f`:

**Strong form:**
```julia
eqs = [-div(grad(u(x,y))) ~ f]
```

**Weak form (after integration by parts):**
```julia
# grad(u) dot grad(v) integrated over domain = f*v integrated over domain
eqs = [Integral(domain)(grad(u(x,y)) dot grad(v(x,y))) ~ Integral(domain)(f*v(x,y))]
```

## Design Philosophy

The variational types in PDEBase are designed to be:

1. **Method descriptors**: Types like `RitzGalerkin` describe the mathematical method, not the implementation
2. **Dispatch hooks**: Enable method dispatch for different formulations
3. **Framework-agnostic**: Work with any FEM library (Ferrite, Gridap, etc.)
4. **Composable**: Can be combined with domain types for complete problem specification

The actual finite element assembly, basis function evaluation, and numerical integration are left to downstream packages that have access to specific mesh and element implementations.
