# PDEBase.jl

PDEBase.jl provides common types and an interface for building discretizers of [ModelingToolkit.jl](https://docs.sciml.ai/ModelingToolkit/stable/) PDESystems. It serves as a foundation for PDE discretization packages in the SciML ecosystem.

## Installation

To install PDEBase.jl, use the Julia package manager:

```julia
using Pkg
Pkg.add("PDEBase")
```

## Purpose

PDEBase.jl is designed to be extended by discretization packages, not used directly for solving PDEs. It provides:

- **Abstract types** for discretizations, discrete spaces, and differential operators
- **Interface functions** that discretization packages implement
- **Utility functions** for symbolic manipulation of PDEs
- **Boundary condition handling** with a structured type hierarchy

If you want to solve PDEs, use a package built on PDEBase:
- [MethodOfLines.jl](https://github.com/SciML/MethodOfLines.jl) - Finite difference discretization

## Key Concepts

### The PDESystem Interface

PDEBase defines how a `ModelingToolkit.PDESystem` is converted into a discretized system (typically an `ODESystem`). Key assumptions about the input PDESystem:

- Well-defined domains with finite intervals for all variables
- Consistent dependent variable signatures throughout
- Boundary conditions specified at domain boundaries
- Initial conditions for time-dependent problems

See [PDESystem Interface](@ref) for complete details.

### The Discretization Workflow

The `symbolic_discretize` function orchestrates the conversion:

1. **Preprocessing** - Handle complex numbers and array variables
2. **Parsing** - Create a `VariableMap` with all variable information
3. **Boundary parsing** - Structure boundary conditions by variable and dimension
4. **Space construction** - Build the discrete grid and variables
5. **Equation processing** - Convert each PDE to discrete equations
6. **System generation** - Combine into the output system

See [Discretization Workflow](@ref workflow) for the complete pipeline.

### Type Hierarchy

```
AbstractDiscretization (from SciMLBase)
├── AbstractEquationSystemDiscretization
└── AbstractOptimizationSystemDiscretization

AbstractDiscreteSpace
└── AbstractCartesianDiscreteSpace

AbstractBoundary
└── AbstractTruncatingBoundary
    ├── LowerBoundary
    ├── UpperBoundary
    └── AbstractInterfaceBoundary
        ├── InterfaceBoundary
        └── HigherOrderInterfaceBoundary
```

## Documentation

```@contents
Pages = [
    "getting_started.md",
    "interface.md",
    "workflow.md",
    "variablemap.md",
    "boundaries.md",
    "api.md",
]
Depth = 2
```

## Quick Reference

### Abstract Types to Extend

| Type | Purpose |
|------|---------|
| `AbstractEquationSystemDiscretization` | For ODE/DAE-producing discretizations |
| `AbstractOptimizationSystemDiscretization` | For optimization problems |
| `AbstractDiscreteSpace` | Discrete spatial representation |
| `AbstractDifferentialDiscretizer` | Differential operator handling |

### Key Interface Functions

| Function | Purpose |
|----------|---------|
| `construct_discrete_space` | Build the discrete grid and variables |
| `construct_differential_discretizer` | Prepare derivative approximations |
| `discretize_equation!` | Convert a PDE to discrete equations |
| `generate_system` | Create the final symbolic system |

### Exported Utility Functions

| Category | Functions |
|----------|-----------|
| Differential analysis | `count_differentials`, `differential_order`, `has_derivatives`, `find_derivative`, `d_orders` |
| Expression manipulation | `get_depvars`, `split_terms`, `split_additive_terms`, `subs_alleqs!` |
| Variable access | `depvars`, `indvars`, `all_ivs`, `ivs`, `depvar`, `x2i` |
| Boundary utilities | `isupper`, `isinterface`, `isperiodic`, `filter_interfaces`, `haslowerupper`, `has_interfaces` |

## Related Packages

- [ModelingToolkit.jl](https://docs.sciml.ai/ModelingToolkit/stable/) - Defines the `PDESystem` type
- [MethodOfLines.jl](https://github.com/SciML/MethodOfLines.jl) - Finite difference discretization built on PDEBase
- [DifferentialEquations.jl](https://docs.sciml.ai/DiffEqDocs/stable/) - Solves the resulting ODE/DAE systems

## Contributing

Contributions are welcome! Please open an issue or pull request on [GitHub](https://github.com/SciML/PDEBase.jl).
