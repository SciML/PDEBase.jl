# Getting Started

PDEBase.jl is a foundational package that provides common types and interfaces for building PDE discretizers in the SciML ecosystem. It is not meant to be used directly for solving PDEs, but rather serves as a base for discretization packages like [MethodOfLines.jl](https://github.com/SciML/MethodOfLines.jl).

## Installation

```julia
using Pkg
Pkg.add("PDEBase")
```

## When to Use PDEBase.jl

You should use PDEBase.jl if you are:

1. **Building a new PDE discretization method** - PDEBase provides the framework for converting `ModelingToolkit.PDESystem` to discretized systems
2. **Extending an existing discretizer** - Understanding PDEBase helps you extend packages like MethodOfLines.jl
3. **Working with symbolic PDE transformations** - The utility functions are helpful for symbolic manipulation of PDEs

If you just want to **solve PDEs**, use a higher-level package:
- [MethodOfLines.jl](https://github.com/SciML/MethodOfLines.jl) - Finite difference methods
- [NeuralPDE.jl](https://github.com/SciML/NeuralPDE.jl) - Physics-informed neural networks

## Package Structure

PDEBase.jl provides:

### Abstract Types

Base types for discretization packages to extend:
- `AbstractEquationSystemDiscretization` - For ODE/DAE-producing discretizations
- `AbstractOptimizationSystemDiscretization` - For optimization discretizations
- `AbstractDiscreteSpace` - For discrete spatial representations
- `AbstractDifferentialDiscretizer` - For differential operator handling

### Concrete Types

Data structures used during discretization:
- `VariableMap` - Organizes variable information from a PDESystem
- `PeriodicMap` - Tracks periodic boundary conditions
- Boundary types (`LowerBoundary`, `UpperBoundary`, etc.)

### Interface Functions

Functions that discretization packages must implement:
- `construct_discrete_space` - Build the discrete grid
- `discretize_equation!` - Convert PDEs to discrete equations
- `generate_system` - Create the final symbolic system

### Utility Functions

Helpers for symbolic manipulation:
- `count_differentials`, `differential_order` - Analyze derivatives
- `get_depvars`, `split_terms` - Parse expressions
- `safe_unwrap`, `recursive_unwrap` - Handle Symbolics types

## Example: Inspecting a PDESystem

While you typically do not use PDEBase directly, you can use it to analyze PDESystems:

```julia
using ModelingToolkit, PDEBase

# Define a simple heat equation
@parameters t x
@variables u(..)
Dt = Differential(t)
Dxx = Differential(x)^2

eq = Dt(u(t,x)) ~ Dxx(u(t,x))
bcs = [u(0,x) ~ sin(pi*x),
       u(t,0) ~ 0,
       u(t,1) ~ 0]
domains = [t ∈ (0, 1), x ∈ (0, 1)]

@named pdesys = PDESystem([eq], bcs, domains, [t, x], [u(t,x)])

# Create a VariableMap to inspect the system
v = PDEBase.VariableMap(pdesys)

# Inspect variables
println("Dependent variables: ", PDEBase.depvars(v))
println("Independent variables: ", PDEBase.indvars(v))
println("Domain intervals: ", v.intervals)

# Analyze the equation
println("Derivative orders in x: ", PDEBase.differential_order(eq, x))
println("Has derivatives: ", PDEBase.has_derivatives(eq.lhs))
```

## Creating a Simple Discretization

Here's a conceptual outline of implementing a discretization:

```julia
using PDEBase, ModelingToolkit, SciMLBase

# 1. Define your discretization type
struct SimpleFiniteDifference <: PDEBase.AbstractEquationSystemDiscretization
    dx::Float64
    t::Any  # time variable
end

# 2. Define your discrete space
struct SimpleSpace <: PDEBase.AbstractCartesianDiscreteSpace
    grid::Vector{Float64}
    u_disc::Vector  # discretized variable
end

# 3. Implement get_time (required)
PDEBase.get_time(d::SimpleFiniteDifference) = d.t

# 4. Implement construct_discrete_space
function PDEBase.construct_discrete_space(v::PDEBase.VariableMap, disc::SimpleFiniteDifference)
    x = first(PDEBase.indvars(v))
    lo, hi = v.intervals[x]
    grid = collect(lo:disc.dx:hi)
    # Create discrete variables...
    u_disc = []  # placeholder for discretized variables
    return SimpleSpace(grid, u_disc)
end

# 5. Implement discretize_equation!
function PDEBase.discretize_equation!(
    disc_state,
    pde,
    vareqmap,
    eqvar,
    bcmap,
    depvars,
    s::SimpleSpace,
    derivweights,
    indexmap,
    disc::SimpleFiniteDifference
)
    # Convert PDE to finite difference equations
    # Add to disc_state
end

# 6. Implement generate_system
function PDEBase.generate_system(
    disc_state,
    s::SimpleSpace,
    u0,
    tspan,
    metadata,
    disc::SimpleFiniteDifference
)
    # Build and return an ODESystem
end
```

## Key Concepts

### The Discretization Pipeline

1. **PDESystem** → Input from ModelingToolkit
2. **VariableMap** → Parsed variable information
3. **Boundary parsing** → Structured boundary conditions
4. **Discrete space** → Grid and discrete variables
5. **Equation discretization** → Convert each PDE
6. **System generation** → Output ODESystem/DAESystem

### Variable Handling

- **Dependent variables** like `u(t,x)` become arrays of unknowns
- **Independent variables** define the grid dimensions
- **Time** is special - it becomes the independent variable of the output ODE

### Boundary Conditions

- Parsed into `LowerBoundary`, `UpperBoundary`, etc.
- Organized by variable and dimension
- Applied during equation discretization

## Next Steps

- Read the [PDESystem Interface](@ref) for detailed interface documentation
- Study [VariableMap](@ref variablemap) to understand variable handling
- Review [Boundary Conditions](@ref boundaries) for BC types
- See [Discretization Workflow](@ref workflow) for the complete pipeline

## Related Packages

- [ModelingToolkit.jl](https://docs.sciml.ai/ModelingToolkit/stable/) - Defines `PDESystem`
- [MethodOfLines.jl](https://github.com/SciML/MethodOfLines.jl) - Production discretizer using PDEBase
- [DifferentialEquations.jl](https://docs.sciml.ai/DiffEqDocs/stable/) - Solves the resulting ODEs
