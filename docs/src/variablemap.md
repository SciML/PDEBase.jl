# [VariableMap](@id variablemap)

The `VariableMap` is a central data structure in PDEBase.jl that stores information about the variables in a PDESystem and provides utilities for accessing and manipulating them.

## Overview

When a PDESystem is passed to `symbolic_discretize`, a `VariableMap` is constructed to organize all the variable information needed during discretization.

## Structure

```julia
struct VariableMap
    ū           # Dependent variables (e.g., [u(t,x), v(t,x)])
    x̄           # Independent spatial variables (e.g., [x, y])
    ps          # Parameters
    time        # Time variable (or nothing if time-independent)
    intervals   # Domain intervals: Dict(x => (lower, upper))
    args        # Argument signatures: Dict(operation(u) => [t, x, y])
    depvar_ops  # Dependent variable operations (the function symbols)
    x2i         # Spatial variable to index mapping
    i2x         # Index to spatial variable mapping
    replaced_vars  # Variables replaced during preprocessing
end
```

## Constructor

A `VariableMap` is typically constructed automatically during discretization:

```julia
v = VariableMap(pdesys, discretization)
```

It can also be constructed without a discretization for inspection:

```julia
v = VariableMap(pdesys)
```

## Field Descriptions

### `ū` - Dependent Variables

A vector of all dependent variables appearing in the system, filtered to exclude boundary evaluations (where arguments are numbers).

**Example**: For a 2D heat equation with `u(t,x,y)`, this would contain the symbolic expression `u(t,x,y)`.

### `x̄` - Independent Spatial Variables

A vector of spatial (non-time) independent variables. These are the dimensions of the domain.

**Example**: For `u(t,x,y)`, this would be `[x, y]` (excluding `t`).

### `ps` - Parameters

Parameters of the PDESystem (from `get_ps(pdesys)`).

### `time` - Time Variable

The time variable from the discretization (via `get_time(disc)`). This is `nothing` for steady-state problems.

### `intervals` - Domain Intervals

A dictionary mapping each independent variable to its domain interval as a tuple `(lower_bound, upper_bound)`.

**Example**:
```julia
v.intervals[x]  # Returns (0.0, 1.0) for x ∈ [0, 1]
v.intervals[t]  # Returns (0.0, 10.0) for t ∈ [0, 10]
```

### `args` - Argument Signatures

A dictionary mapping each dependent variable's operation (function symbol) to its argument list.

**Example**:
```julia
v.args[operation(u)]  # Returns [t, x, y] for u(t,x,y)
```

### `depvar_ops` - Dependent Variable Operations

The symbolic operations (function symbols) for all dependent variables. Used to identify dependent variables in expressions.

### `x2i` and `i2x` - Index Mappings

Bidirectional mappings between spatial variables and their indices:
- `x2i`: `Dict(x => 1, y => 2)` - variable to index
- `i2x`: `Dict(1 => x, 2 => y)` - index to variable

## Accessor Functions

PDEBase exports several functions for working with `VariableMap`:

### `depvars(v::VariableMap)`

Returns all dependent variables.

```julia
depvars(v)  # Returns [u(t,x,y), v(t,x,y), ...]
```

### `indvars(v::VariableMap)`

Returns spatial independent variables (excludes time).

```julia
indvars(v)  # Returns [x, y, ...]
```

### `all_ivs(v::VariableMap)`

Returns all independent variables including time.

```julia
all_ivs(v)  # Returns [x, y, t] or [x, y] if no time
```

### `ivs(u, v::VariableMap)`

Returns the independent variables for a specific dependent variable, excluding time.

```julia
ivs(u, v)  # Returns [x, y] for u(t,x,y)
```

### `all_ivs(u, v::VariableMap)`

Returns all independent variables for a specific dependent variable.

```julia
all_ivs(u, v)  # Returns [t, x, y] for u(t,x,y)
```

### `depvar(u, v::VariableMap)`

Returns the full dependent variable expression with all its arguments.

```julia
depvar(u, v)  # Returns u(t,x,y) from any partial reference
```

### `x2i(v::VariableMap, u, x)`

Returns the index of spatial variable `x` in the argument list of variable `u`.

```julia
x2i(v, u, x)  # Returns 1 if x is the first spatial argument
```

### `Base.ndims(u, v::VariableMap)`

Returns the number of spatial dimensions for a dependent variable.

```julia
ndims(u, v)  # Returns 2 for u(t,x,y)
```

## Usage Example

```julia
using ModelingToolkit, PDEBase

@parameters t x y
@variables u(..) v(..)
Dt = Differential(t)
Dx = Differential(x)
Dy = Differential(y)

# Define a simple system
eqs = [Dt(u(t,x,y)) ~ Dx(Dx(u(t,x,y))) + Dy(Dy(u(t,x,y)))]
bcs = [u(0,x,y) ~ sin(pi*x)*sin(pi*y),
       u(t,0,y) ~ 0, u(t,1,y) ~ 0,
       u(t,x,0) ~ 0, u(t,x,1) ~ 0]
domains = [t ∈ (0, 1), x ∈ (0, 1), y ∈ (0, 1)]

@named pdesys = PDESystem(eqs, bcs, domains, [t, x, y], [u(t,x,y)])

# Create VariableMap (without discretization for inspection)
v = VariableMap(pdesys)

# Access variable information
depvars(v)         # [u(t,x,y)]
indvars(v)         # [t, x, y] (all variables when created without discretization)
all_ivs(v)         # [t, x, y]
v.intervals[x]     # (0, 1)
v.intervals[y]     # (0, 1)
```

## Internal Functions

### `update_varmap!(v, newdv)`

Adds a new dependent variable to an existing VariableMap. Used internally when the system is transformed.

### `axiesvals(v, u, x, I)`

Returns axis values at grid positions. Used internally during discretization.

## See Also

- [PDESystem Interface](@ref) for how VariableMap is used in discretization
- [Boundary Conditions](@ref boundaries) for how boundaries reference variables
