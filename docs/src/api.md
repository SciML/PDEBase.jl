# API Reference

This page provides a comprehensive reference for all exported types and functions in PDEBase.jl.

## Abstract Types

### Discretization Types

**`AbstractEquationSystemDiscretization <: AbstractDiscretization`**

Base type for discretizations that produce systems of ODEs or DAEs.

**`AbstractOptimizationSystemDiscretization <: AbstractDiscretization`**

Base type for discretizations targeting optimization problems.

### Space Types

**`AbstractDiscreteSpace`**

Base type for discrete spatial representations.

**`AbstractCartesianDiscreteSpace <: AbstractDiscreteSpace`**

Specialized type for Cartesian grid discretizations.

### Auxiliary Types

**`AbstractVarEqMapping`**

Abstract type for mappings from variables to the equations they are solved from.

**`AbstractDifferentialDiscretizer`**

Abstract type for objects that discretize differential operators.

**`AbstractDiscretizationState`**

Abstract type for mutable state containers used during discretization.

## Boundary Condition Types

**`AbstractBoundary`**

Root abstract type for all boundary conditions.

**`AbstractTruncatingBoundary <: AbstractBoundary`**

Abstract type for boundaries at domain edges.

**`AbstractInterfaceBoundary <: AbstractTruncatingBoundary`**

Abstract type for interface boundaries between regions.

**`LowerBoundary <: AbstractTruncatingBoundary`**

Boundary condition at the lower edge of a spatial domain. Fields:
- `u` - Dependent variable
- `x` - Independent variable (spatial dimension)
- `depvars` - All dependent variables in the BC
- `indvars` - Independent variables in the BC
- `eq` - The boundary condition equation
- `order` - Derivative order

**`UpperBoundary <: AbstractTruncatingBoundary`**

Boundary condition at the upper edge of a spatial domain. Same fields as `LowerBoundary`.

**`InterfaceBoundary{IsUpper_u, IsUpper_u2} <: AbstractInterfaceBoundary`**

Interface condition connecting two regions at zero-th order. Fields:
- `u`, `u2` - The two variables
- `x`, `x2` - Interface coordinates
- `eq` - Interface equation

**`HigherOrderInterfaceBoundary <: AbstractInterfaceBoundary`**

Interface condition with derivatives (e.g., flux continuity). Additional fields include `depvars`, `indvars`, and `order`.

## Concrete Types

**`VariableMap`**

Central data structure storing information about PDESystem variables. See [VariableMap](@ref variablemap) for details.

Fields:
- `ū` - Dependent variables
- `x̄` - Independent spatial variables
- `ps` - Parameters
- `time` - Time variable
- `intervals` - Domain intervals
- `args` - Argument signatures
- `depvar_ops` - Dependent variable operations
- `x2i`, `i2x` - Index mappings
- `replaced_vars` - Variables replaced during preprocessing

**`PeriodicMap`**

Tracks which variables have periodic boundary conditions.

## Interface Functions

For detailed documentation of interface functions, see [PDESystem Interface](@ref).

### Summary

| Function | Purpose |
|----------|---------|
| `interface_errors` | Check PDESystem compatibility |
| `check_boundarymap` | Validate boundary conditions |
| `should_transform` | Determine if system needs transformation |
| `transform_pde_system!` | Transform PDESystem |
| `construct_disc_state` | Create discretization state container |
| `construct_discrete_space` | Build discrete grid and variables |
| `construct_var_equation_mapping` | Map variables to equations |
| `construct_differential_discretizer` | Prepare derivative approximations |
| `discretize_equation!` | Convert PDE to discrete equations |
| `generate_ic_defaults` | Create initial conditions |
| `generate_metadata` | Build discretization metadata |
| `generate_system` | Create final symbolic system |
| `get_time` | Get time variable |
| `get_discvars` | Get discrete variables |
| `get_eqvar` | Get equation variable |
| `add_metadata!` | Add to metadata |

## VariableMap Accessor Functions

| Function | Purpose |
|----------|---------|
| `depvars(v)` | Get all dependent variables |
| `indvars(v)` | Get spatial independent variables |
| `all_ivs(v)` | Get all independent variables |
| `ivs(u, v)` | Get independent variables for a dependent variable |
| `depvar(u, v)` | Get full dependent variable expression |
| `x2i(v, u, x)` | Get index of spatial variable |

## Boundary Utility Functions

| Function | Purpose |
|----------|---------|
| `isupper(b)` | Check if upper boundary |
| `isinterface(b)` | Check if interface boundary |
| `isperiodic(b1, b2)` | Check if periodic pair |
| `filter_interfaces(bs)` | Get interface boundaries |
| `haslowerupper(bs, x)` | Check for lower/upper boundaries |
| `has_interfaces(bmps)` | Check for interfaces in boundary map |
| `getvars(b)` | Get boundary variables |
| `flatten_vardict(bmps)` | Flatten boundary map |

## Symbolic Utility Functions

### Differential Analysis

| Function | Purpose |
|----------|---------|
| `count_differentials(term, x)` | Count derivatives w.r.t. variable |
| `differential_order(eq, x)` | Get derivative orders in equation |
| `has_derivatives(term)` | Check for derivatives |
| `find_derivative(term, op)` | Find derivative in expression |
| `d_orders(x, eqs)` | Get derivative orders across equations |

### Expression Manipulation

| Function | Purpose |
|----------|---------|
| `get_depvars(eq, ops)` | Find dependent variables |
| `get_all_depvars(eqs, ops)` | Find all dependent variables |
| `split_terms(eq)` | Split equation into terms |
| `split_additive_terms(eq)` | Split by addition |
| `subs_alleqs!(eqs, rules)` | Substitute in all equations |
| `subsmatch(expr, rule)` | Check if rule matches |
| `ex2term(term, v)` | Convert to term variable |
| `safe_unwrap(x)` | Unwrap Num types |
| `recursive_unwrap(ex)` | Recursively unwrap |
