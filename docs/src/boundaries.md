# [Boundary Conditions](@id boundaries)

PDEBase.jl provides a type hierarchy for representing boundary conditions that arise during the parsing of a PDESystem. Understanding these types is essential for implementing discretization methods that handle boundaries correctly.

## Type Hierarchy

```
AbstractBoundary
└── AbstractTruncatingBoundary
    ├── LowerBoundary
    ├── UpperBoundary
    └── AbstractInterfaceBoundary
        ├── InterfaceBoundary
        └── HigherOrderInterfaceBoundary
```

## Boundary Types

### `AbstractBoundary`

The root abstract type for all boundary conditions.

### `AbstractTruncatingBoundary`

Abstract type for boundaries that "truncate" the computational domain - i.e., boundaries at the edges of the spatial domain.

### `LowerBoundary`

Represents a boundary condition at the lower edge of a spatial domain (e.g., `x = 0`).

```julia
struct LowerBoundary <: AbstractTruncatingBoundary
    u        # Dependent variable this BC applies to
    x        # Independent variable (spatial dimension)
    depvars  # All dependent variables appearing in the BC equation
    indvars  # Independent variables in the BC
    eq       # The boundary condition equation
    order    # Order of the highest derivative in the BC
end
```

**Example**: For `u(t, 0) = 0`, this would create a `LowerBoundary` for variable `u` at `x = 0`.

### `UpperBoundary`

Represents a boundary condition at the upper edge of a spatial domain (e.g., `x = 1`).

```julia
struct UpperBoundary <: AbstractTruncatingBoundary
    u        # Dependent variable this BC applies to
    x        # Independent variable (spatial dimension)
    depvars  # All dependent variables appearing in the BC equation
    indvars  # Independent variables in the BC
    eq       # The boundary condition equation
    order    # Order of the highest derivative in the BC
end
```

**Example**: For `Dx(u(t, 1)) = 0` (Neumann BC), this would be an `UpperBoundary` with `order = 1`.

### `InterfaceBoundary`

Represents an interface condition connecting two regions or variables. These are used for multi-domain problems where solutions must match at internal boundaries.

```julia
struct InterfaceBoundary{IsUpper_u, IsUpper_u2} <: AbstractInterfaceBoundary
    u    # First variable
    u2   # Second variable
    x    # First interface coordinate
    x2   # Second interface coordinate
    eq   # Interface condition equation
end
```

The type parameters `IsUpper_u` and `IsUpper_u2` indicate whether each variable is at its upper boundary.

**Key Assumptions**:
- Variables in an interface BC must have the same argument signature
- They differ in only one variable (the interface variable)
- Interfaces are typically at the lower boundary of one domain and upper of another

### `HigherOrderInterfaceBoundary`

For interface conditions involving derivatives (e.g., flux continuity).

```julia
struct HigherOrderInterfaceBoundary <: AbstractInterfaceBoundary
    u        # First variable
    u2       # Second variable
    x        # First interface coordinate
    x2       # Second interface coordinate
    depvars  # Dependent variables in the condition
    indvars  # Independent variables
    eq       # The interface equation
    order    # Derivative order
end
```

## Boundary Map Structure

During discretization, boundary conditions are organized into a nested dictionary structure:

```julia
boundarymap = Dict(
    operation(u) => Dict(
        x => [LowerBoundary(...), UpperBoundary(...)],
        y => [LowerBoundary(...), UpperBoundary(...)],
        t => [initial_condition, ...]
    ),
    operation(v) => Dict(...)
)
```

This structure maps each dependent variable to its boundaries for each independent variable.

## Utility Functions

### `isupper(b::AbstractBoundary)`

Returns `true` if the boundary is an upper boundary.

```julia
isupper(LowerBoundary(...))  # false
isupper(UpperBoundary(...))  # true
```

### `isinterface(b)`

Returns `Val(true)` if the boundary is an interface boundary, `Val(false)` otherwise.

```julia
isinterface(InterfaceBoundary(...))  # Val(true)
isinterface(LowerBoundary(...))      # Val(false)
```

### `isperiodic(b1::InterfaceBoundary, b2::InterfaceBoundary)`

Checks if two interface boundaries form a periodic pair.

```julia
isperiodic(bc_left, bc_right)  # true if periodic
```

### `filter_interfaces(bs)`

Filters a collection of boundaries to return only interface boundaries.

```julia
interfaces = filter_interfaces(all_boundaries)
```

### `haslowerupper(bs, x)`

Checks if a collection of boundaries has both lower and upper boundaries for variable `x`.

```julia
has_lower, has_upper = haslowerupper(boundaries, x)
```

### `has_interfaces(bmps)`

Returns `true` if the boundary map contains any interface boundaries.

```julia
if has_interfaces(boundarymap)
    # Handle multi-domain case
end
```

### `getvars(b::AbstractBoundary)`

Returns the dependent variable and independent variable for a boundary.

```julia
u, x = getvars(boundary)
```

### `flatten_vardict(bmps)`

Flattens a nested boundary map dictionary into a single vector of boundaries.

```julia
all_boundaries = flatten_vardict(boundarymap)
```

## Parsing Boundary Conditions

The `parse_bcs` function (internal) processes the raw boundary conditions from the PDESystem:

1. Generates matching rules for lower and upper boundaries
2. Iterates through each boundary condition equation
3. Determines which boundary (lower/upper/interface) each condition applies to
4. Detects periodic boundaries (conditions that reference both lower and upper bounds)
5. Constructs the appropriate boundary type
6. Builds the boundary map

### How Boundaries are Matched

A boundary condition is matched based on:
1. Whether dependent variables are evaluated at domain boundaries
2. The presence of derivatives and their orders
3. Whether the condition spans multiple boundaries (interface/periodic)

**Example of boundary matching**:
```julia
# This is detected as a LowerBoundary at x = 0:
u(t, 0) ~ 0

# This is detected as an UpperBoundary at x = 1:
Dx(u(t, 1)) ~ 0

# This is detected as an InterfaceBoundary (periodic):
u(t, 0) ~ u(t, 1)
```

## Implementation Notes

When implementing a discretization:

1. **Check boundary validity** in `check_boundarymap`:
```julia
function PDEBase.check_boundarymap(bmap, v, disc::MyDiscretization)
    for u in keys(bmap)
        for x in keys(bmap[u])
            if isempty(bmap[u][x])
                error("Missing boundary condition for $u at $x")
            end
        end
    end
end
```

2. **Handle different boundary types** in `discretize_equation!`:
```julia
for bc in bcmap[operation(eqvar)][x]
    if bc isa LowerBoundary
        # Apply at first grid point
    elseif bc isa UpperBoundary
        # Apply at last grid point
    elseif bc isa InterfaceBoundary
        # Connect multiple regions
    end
end
```

3. **Separate initial conditions from spatial BCs**:
```julia
# Initial conditions are boundaries on the time variable
ics = boundarymap[u][t]

# Spatial boundaries are on spatial variables
spatial_bcs = Dict(x => boundarymap[u][x] for x in indvars(v))
```

## PeriodicMap

The `PeriodicMap` type tracks which variables have periodic boundary conditions:

```julia
struct PeriodicMap
    map  # Dict tracking periodic dimensions for each variable
end
```

This is used to determine when to wrap grid indices in finite difference stencils.

## See Also

- [PDESystem Interface](@ref) for how boundaries integrate with discretization
- [Discretization Workflow](@ref workflow) for the complete processing pipeline
- [VariableMap](@ref variablemap) for variable handling
