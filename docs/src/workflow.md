# [Discretization Workflow](@id workflow)

This page describes the complete workflow of `symbolic_discretize`, the main entry point for converting a `PDESystem` into a discretized system.

## Overview

The `symbolic_discretize` function orchestrates the entire discretization process through several distinct phases:

1. **Preprocessing** - Handle complex numbers and array variables
2. **Parsing** - Extract and organize variable information
3. **Validation** - Check system compatibility with the discretization
4. **Transformation** (optional) - Modify the system for compatibility
5. **Discretization Setup** - Build discrete space and operators
6. **Equation Processing** - Convert PDEs to discrete equations
7. **Finalization** - Generate the output system

## Detailed Workflow

### Phase 1: Preprocessing

```julia
pdesys, complexmap = handle_complex(pdesys)
cardinalize_eqs!(pdesys)
pdesys, replaced_vars = make_pdesys_compatible(pdesys)
```

**`handle_complex`**: Splits complex-valued equations into real and imaginary parts.

**`cardinalize_eqs!`**: Normalizes equations to the form `lhs - rhs ~ 0`.

**`make_pdesys_compatible`**: Handles array-valued variables, expanding them into scalar equations.

### Phase 2: System Parsing

```julia
v = VariableMap(pdesys, discretization, replaced_vars = replaced_vars)
```

Creates a `VariableMap` that stores:
- All dependent and independent variables
- Domain intervals
- Variable argument signatures
- Index mappings

See [VariableMap](@ref variablemap) for details.

### Phase 3: Validation

```julia
interface_errors(pdesys, v, discretization)
```

Checks for fundamental incompatibilities between the PDESystem and the discretization method. Discretization packages should override this to catch unsupported features early.

### Phase 4: Boundary Parsing

```julia
tspan = t !== nothing ? v.intervals[t] : nothing
bcorders = Dict(map(x -> x => d_orders(x, get_bcs(pdesys)), all_ivs(v)))
boundarymap = parse_bcs(get_bcs(pdesys), v, bcorders)
check_boundarymap(boundarymap, v, discretization)
```

1. **Extract time span** (if time-dependent)
2. **Find derivative orders** in boundary conditions
3. **Parse boundary conditions** into structured `boundarymap`
4. **Validate boundaries** for the specific discretization

The boundary parsing identifies:
- Which variable each BC applies to
- Which boundary (lower/upper) it's on
- Whether it's a periodic or interface condition
- The derivative order involved

### Phase 5: System Transformation (Optional)

```julia
if should_transform(pdesys, discretization, boundarymap)
    pdesys = transform_pde_system!(v, boundarymap, pdesys, discretization)
end
```

Some discretizations may need to modify the PDESystem. For example:
- Adding auxiliary equations
- Transforming coordinate systems
- Splitting higher-order equations

### Phase 6: Discretization Setup

```julia
disc_state = construct_disc_state(discretization)
s = construct_discrete_space(v, discretization)
vareqmap = construct_var_equation_mapping(pdeeqs, boundarymap, s, discretization)
orders = Dict(map(x -> x => collect(union(pdeorders[x], bcorders[x])), indvars(v)))
derivweights = construct_differential_discretizer(pdesys, s, discretization, orders)
```

1. **`construct_disc_state`**: Creates a mutable container to accumulate discretized equations
2. **`construct_discrete_space`**: Builds the discrete grid/mesh and discretized variables
3. **`construct_var_equation_mapping`**: Determines which equation solves for which variable
4. **`construct_differential_discretizer`**: Pre-computes derivative approximations (e.g., finite difference stencils)

### Phase 7: Initial Condition Processing

```julia
ics = t === nothing ? [] : mapreduce(u -> boundarymap[u][t], vcat, operation.(depvars(v)))

bcmap = Dict(map(collect(keys(boundarymap))) do u
    u => Dict(map(indvars(v)) do x
        x => boundarymap[u][x]
    end)
end)
```

Separates initial conditions (on time) from spatial boundary conditions.

### Phase 8: Equation Discretization Loop

```julia
for pde in pdeeqs
    # Extract dependent variables from equation
    depvars_lhs = get_depvars(pde.lhs, v.depvar_ops)
    depvars_rhs = get_depvars(pde.rhs, v.depvar_ops)
    depvars = collect(depvars_lhs ∪ depvars_rhs)
    depvars = filter(u -> !any(map(x -> x isa Number, arguments(u))), depvars)

    # Get the variable this equation solves for
    eqvar = get_eqvar(vareqmap, pde)

    # Build index map for this equation
    args = ivs(eqvar, v)
    indexmap = Dict([args[i] => i for i in 1:length(args)])

    # Discretize the equation
    discretize_equation!(disc_state, pde, vareqmap, eqvar, bcmap,
        depvars, s, derivweights, indexmap, discretization)
end
```

For each PDE:
1. Find all dependent variables in the equation
2. Determine which variable the equation solves for
3. Build an index map for spatial dimensions
4. Call `discretize_equation!` to add discrete equations to `disc_state`

The `discretize_equation!` function is responsible for:
- Converting the PDE to discrete equations (e.g., at interior grid points)
- Handling boundary conditions for this variable
- Adding all generated equations to `disc_state`

### Phase 9: Finalization

```julia
u0 = generate_ic_defaults(ics, s, discretization)
metadata = generate_metadata(s, discretization, pdesys, boundarymap, complexmap)
return generate_system(disc_state, s, u0, tspan, metadata, discretization)
```

1. **`generate_ic_defaults`**: Creates initial condition values for the discrete variables
2. **`generate_metadata`**: Builds metadata for solution reconstruction (e.g., grid information)
3. **`generate_system`**: Combines everything into the final `ODESystem`, `DAESystem`, or other system type

## Data Flow Diagram

```
PDESystem
    │
    ▼
┌─────────────────────────────────────┐
│  Preprocessing                      │
│  - handle_complex                   │
│  - cardinalize_eqs!                 │
│  - make_pdesys_compatible           │
└─────────────────────────────────────┘
    │
    ▼
┌─────────────────────────────────────┐
│  Parsing                            │
│  - VariableMap creation             │
│  - Boundary condition parsing       │
└─────────────────────────────────────┘
    │
    ▼
┌─────────────────────────────────────┐
│  Setup                              │
│  - Discrete space construction      │
│  - Differential discretizer         │
│  - Variable-equation mapping        │
└─────────────────────────────────────┘
    │
    ▼
┌─────────────────────────────────────┐
│  Equation Loop                      │
│  For each PDE:                      │
│    - discretize_equation!           │
│    - Handle BCs                     │
└─────────────────────────────────────┘
    │
    ▼
┌─────────────────────────────────────┐
│  Finalization                       │
│  - generate_ic_defaults             │
│  - generate_metadata                │
│  - generate_system                  │
└─────────────────────────────────────┘
    │
    ▼
ODESystem / DAESystem / OptimizationSystem
```

## Key Extension Points

When implementing a new discretization, the main functions to implement are:

| Function | Purpose | Typical Implementation |
|----------|---------|------------------------|
| `construct_discrete_space` | Build grid and discrete variables | Create mesh, discretize variables |
| `construct_differential_discretizer` | Prepare derivative approximations | Compute finite difference weights |
| `discretize_equation!` | Convert PDE to discrete equations | Apply stencils, generate equations |
| `generate_system` | Create final symbolic system | Combine equations into `ODESystem` |

## Symbolic Utilities

PDEBase provides several utility functions for working with symbolic expressions during discretization:

### Differential Analysis

- `count_differentials(term, x)` - Count derivatives w.r.t. variable `x`
- `differential_order(eq, x)` - Get all derivative orders for `x` in an equation
- `has_derivatives(term)` - Check if term contains any derivatives
- `find_derivative(term, depvar_op)` - Find derivative or dependent variable in expression
- `d_orders(x, eqs)` - Get all derivative orders for `x` across equations

### Expression Manipulation

- `split_terms(eq)` - Split equation into individual terms
- `split_additive_terms(eq)` - Split only by addition
- `get_depvars(eq, depvar_ops)` - Find all dependent variables in expression
- `ex2term(term, v)` - Convert expression to a new term variable
- `subs_alleqs!(eqs, rules)` - Substitute rules in all equations

### Utility Functions

- `safe_unwrap(x)` - Safely unwrap `Num` types
- `recursive_unwrap(ex)` - Recursively unwrap nested expressions
- `subsmatch(expr, rule)` - Check if a substitution rule matches

## Error Handling

The workflow includes several validation points:

1. **`interface_errors`** - Catch unsupported PDESystem features
2. **`check_boundarymap`** - Validate boundary conditions
3. **Assertions during parsing** - Ensure boundary conditions are on domain boundaries

Implement these in your discretization to provide helpful error messages.

## See Also

- [PDESystem Interface](@ref) for interface function details
- [VariableMap](@ref variablemap) for variable handling
- [Boundary Conditions](@ref boundaries) for boundary types
