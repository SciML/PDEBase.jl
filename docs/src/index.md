# PDEBase.jl

PDEBase.jl provides common types and an interface for building discretizers of [ModelingToolkit.jl](https://docs.sciml.ai/ModelingToolkit/stable/) PDESystems. It serves as a foundation for PDE discretization packages in the SciML ecosystem.

## Installation

To install PDEBase.jl, use the Julia package manager:

```julia
using Pkg
Pkg.add("PDEBase")
```

## Overview

PDEBase.jl defines the core abstractions and utilities needed for discretizing partial differential equations. The main components include:

### Abstract Types

- `AbstractEquationSystemDiscretization`: Base type for discretizations that produce equation systems (ODEs, DAEs)
- `AbstractOptimizationSystemDiscretization`: Base type for discretizations targeting optimization problems
- `AbstractDiscreteSpace`: Base type for discrete spatial representations
- `AbstractCartesianDiscreteSpace`: Specialized type for Cartesian grid discretizations
- `AbstractVarEqMapping`: Type for variable-to-equation mappings
- `AbstractDifferentialDiscretizer`: Type for differential operator discretizers
- `AbstractDiscretizationState`: Type for storing discretization state during processing

### Boundary Condition Types

- `LowerBoundary`: Represents a boundary condition at the lower boundary
- `UpperBoundary`: Represents a boundary condition at the upper boundary
- `InterfaceBoundary`: Represents an interface boundary condition between regions
- `HigherOrderInterfaceBoundary`: Represents higher-order interface conditions

### Utility Types

- `VariableMap`: Maps dependent and independent variables with their properties
- `PeriodicMap`: Tracks periodic boundary conditions for variables

## Interface Functions

Packages building on PDEBase.jl should implement these interface functions:

```@docs
PDEBase.interface_errors
PDEBase.check_boundarymap
PDEBase.should_transform
PDEBase.transform_pde_system!
PDEBase.construct_disc_state
PDEBase.construct_discrete_space
PDEBase.construct_var_equation_mapping
PDEBase.construct_differential_discretizer
PDEBase.discretize_equation!
PDEBase.generate_ic_defaults
PDEBase.generate_metadata
PDEBase.generate_system
PDEBase.get_time
PDEBase.get_discvars
PDEBase.get_eqvar
PDEBase.add_metadata!
```

## Related Packages

- [MethodOfLines.jl](https://github.com/SciML/MethodOfLines.jl): A finite difference discretization package built on PDEBase.jl
- [ModelingToolkit.jl](https://docs.sciml.ai/ModelingToolkit/stable/): The symbolic modeling framework that defines PDESystem

## Contributing

If you are using this package or want to contribute, please open an issue on [GitHub](https://github.com/SciML/PDEBase.jl).
