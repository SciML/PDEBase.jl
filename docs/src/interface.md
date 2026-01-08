# PDESystem Interface

PDEBase.jl provides a set of interface functions that discretization packages must implement to convert a `ModelingToolkit.PDESystem` into a discretized system (typically an `ODESystem` or `OptimizationSystem`).

## Key Assumptions about PDESystem

When working with PDEBase.jl, the following assumptions are made about the input `PDESystem`:

### Domain Requirements
- The PDESystem must have a well-defined domain with intervals for all independent variables
- Domain boundaries must be finite (using `DomainSets.jl` interval notation)
- Each independent variable must have exactly one domain interval defined

### Variable Requirements
- **Dependent variables** (e.g., `u(t,x,y)`) are functions of independent variables
- **Independent variables** include spatial variables (e.g., `x`, `y`, `z`) and optionally a time variable
- All dependent variables appearing in equations must be registered in the PDESystem
- Variable arguments must be consistent across all uses (e.g., `u(t,x)` cannot appear as `u(x,t)` elsewhere)

### Equation Requirements
- The system must include PDEs (partial differential equations) as the main equations
- Boundary conditions must be specified for all spatial boundaries
- Initial conditions are required if the system is time-dependent
- All boundary conditions must be defined at the domain boundaries (not in the interior)

### Boundary Condition Requirements
- Each dependent variable must have boundary conditions for each spatial dimension
- Boundary conditions must evaluate on the lower or upper bound of each spatial domain
- Interface conditions (connecting multiple regions) must be properly structured
- Periodic boundary conditions must appear in pairs (lower and upper)

## Abstract Types

PDEBase.jl defines several abstract types that discretization packages should subtype:

### Discretization Types

```julia
abstract type AbstractEquationSystemDiscretization <: AbstractDiscretization end
```
Use this as a supertype for discretizations that produce systems of ODEs or DAEs.

```julia
abstract type AbstractOptimizationSystemDiscretization <: AbstractDiscretization end
```
Use this for discretizations targeting optimization problems (e.g., optimal control).

### Space Representation Types

```julia
abstract type AbstractDiscreteSpace end
abstract type AbstractCartesianDiscreteSpace <: AbstractDiscreteSpace end
```
These represent the discretized spatial domain. `AbstractCartesianDiscreteSpace` is for regular Cartesian grids.

### Auxiliary Types

```julia
abstract type AbstractVarEqMapping end
```
Maps variables to the equations they are solved from.

```julia
abstract type AbstractDifferentialDiscretizer end
```
Stores information needed to discretize differential operators.

```julia
abstract type AbstractDiscretizationState end
```
Holds mutable state during the discretization process.

## Interface Functions to Implement

To create a new discretization package, you must implement the following interface functions. Default implementations are provided that return `nothing`, but most must be overridden for a functional discretizer.

### Validation Functions

```@docs
PDEBase.interface_errors
PDEBase.check_boundarymap
```

### Transformation Functions

```@docs
PDEBase.should_transform
PDEBase.transform_pde_system!
```

### Construction Functions

```@docs
PDEBase.construct_disc_state
PDEBase.construct_discrete_space
PDEBase.construct_var_equation_mapping
PDEBase.construct_differential_discretizer
```

### Discretization Functions

```@docs
PDEBase.discretize_equation!
PDEBase.generate_ic_defaults
```

### Finalization Functions

```@docs
PDEBase.generate_metadata
PDEBase.generate_system
```

### Utility Functions

```@docs
PDEBase.get_time
PDEBase.get_discvars
PDEBase.get_eqvar
PDEBase.add_metadata!
```

## Implementation Example

Here is a minimal skeleton for implementing a new discretization:

```julia
using PDEBase, ModelingToolkit, SciMLBase

# Define your discretization type
struct MyDiscretization <: PDEBase.AbstractEquationSystemDiscretization
    dx::Float64  # grid spacing
    # ... other parameters
end

# Define your discrete space
struct MyDiscreteSpace <: PDEBase.AbstractCartesianDiscreteSpace
    grid::Vector{Float64}
    vars::Dict
    # ... other fields
end

# Implement required interface functions
function PDEBase.construct_discrete_space(v::PDEBase.VariableMap, disc::MyDiscretization)
    # Build the discrete grid and variables
    # ...
    return MyDiscreteSpace(grid, vars)
end

function PDEBase.discretize_equation!(
    disc_state,
    pde::Equation,
    vareqmap,
    eqvar,
    bcmap,
    depvars,
    s::MyDiscreteSpace,
    derivweights,
    indexmap,
    disc::MyDiscretization
)
    # Convert the PDE to discretized equations
    # Add equations to disc_state
    # ...
end

function PDEBase.generate_system(
    disc_state,
    s::MyDiscreteSpace,
    u0,
    tspan,
    metadata,
    disc::MyDiscretization
)
    # Combine all equations into an ODESystem or similar
    # Extract time variable and build system
    t = PDEBase.get_time(disc)
    # eqs = ... (build equations from disc_state)
    return ODESystem(eqs, t; name=:discretized)
end
```

## See Also

- [Discretization Workflow](@ref workflow) for the complete processing pipeline
- [VariableMap](@ref variablemap) for understanding variable handling
- [Boundary Conditions](@ref boundaries) for boundary condition types
