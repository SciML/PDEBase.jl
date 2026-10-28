module PDEBase
import DomainSets
using IntervalSets: var".."
import ModelingToolkitBase
import PrecompileTools
import SciMLBase
import Symbolics
import SymbolicUtils
using ModelingToolkitBase: PDESystem, ProblemTypeCtx, System, get_bcs, get_connector_type,
    get_domain, get_dvs, get_eqs, get_gui_metadata, get_iv, get_ivs, get_metadata, get_ps,
    get_systems, get_unknowns, initialization_equations, mtkcompile
using PrecompileTools: @compile_workload, @setup_workload
using SciMLBase: AbstractDiscretization, AbstractDiscretizationMetadata, NonlinearProblem,
    ODEFunction, ODEProblem
import SciMLPublic: @public
using Symbolics: @variables, Differential, Equation, Num, unwrap
using SymbolicUtils: @rule, arguments, getmetadata, iscall, operation, substitute,
    unwrap_const
using SymbolicUtils.Rewriters: Chain, Prewalk
using SymbolicIndexingInterface: is_time_dependent
using TermInterface: maketerm, metadata

"""
    AbstractEquationSystemDiscretization <: AbstractDiscretization

Supertype for discretizations that lower a `PDESystem` into an equation-based
SciML system, typically an `ODESystem` or `DAESystem`.

Discretizer packages subtype this type and implement the developer extension
hooks documented in the [Developer Extension API](@ref developer_extension_api).
"""
abstract type AbstractEquationSystemDiscretization <: AbstractDiscretization end

"""
    AbstractBoxDiscretization <: AbstractEquationSystemDiscretization

Supertype for equation-system discretizations whose spatial domain is an
interval or a product of intervals, such as finite difference and
pseudospectral methods.

Boundary conditions reach these discretizations in endpoint-substituted form,
for example `u(t, 0, y)` on the lower `x` face; `lower_boundaries` is the hook
that rewrites other boundary spellings into that form.
"""
abstract type AbstractBoxDiscretization <: AbstractEquationSystemDiscretization end

"""
    AbstractMeshDiscretization <: AbstractEquationSystemDiscretization

Supertype for equation-system discretizations whose spatial domain is a
general mesh, such as finite element methods.

Boundary conditions for these discretizations are keyed by named boundary
pieces of the mesh rather than by interval endpoints.
"""
abstract type AbstractMeshDiscretization <: AbstractEquationSystemDiscretization end

"""
    AbstractOptimizationSystemDiscretization <: AbstractDiscretization

Supertype for discretizations that lower a `PDESystem` into an optimization
system.

Use this type when the generated system is consumed by an optimization problem
rather than by an ODE or DAE solver.
"""
abstract type AbstractOptimizationSystemDiscretization <: AbstractDiscretization end

"""
    AbstractDiscreteSpace

Supertype for discretized spatial domain representations produced from a
`PDESystem` and a discretization.

Concrete spaces are returned by `construct_discrete_space` and are passed to
the remaining discretization hooks.
"""
abstract type AbstractDiscreteSpace end

"""
    AbstractCartesianDiscreteSpace <: AbstractDiscreteSpace

Supertype for discrete space representations whose coordinates form a
Cartesian product grid.

Subtype this type when `x2i` and Cartesian indexing describe the spatial
coordinates without a general mesh connectivity object.
"""
abstract type AbstractCartesianDiscreteSpace <: AbstractDiscreteSpace end

"""
    AbstractVarEqMapping

Supertype for mappings that assign dependent variables and boundary conditions
to the equations used to discretize them.

Concrete mappings are returned by `construct_var_equation_mapping`; their
`get_eqvar` method must identify the discrete variable for every PDE.
"""
abstract type AbstractVarEqMapping end

"""
    AbstractDifferentialDiscretizer

Supertype for objects that hold the data needed to discretize differential
operators on a discrete space.

Concrete derivative data is returned by `construct_differential_discretizer`
and consumed by `discretize_equation!`.
"""
abstract type AbstractDifferentialDiscretizer end

"""
    AbstractDiscretizationState

Supertype for mutable state accumulated while PDE equations and boundary
conditions are discretized.

Concrete state is returned by `construct_disc_state` and is updated in place by
`discretize_equation!` before `generate_system` consumes it.
"""
abstract type AbstractDiscretizationState end

include("fallbacks.jl")
include("symbolic_utils.jl")
include("utils.jl")
include("mesh_domain.jl")
include("variable_map.jl")
include("classify.jl")
include("box_lowering.jl")
include("interface_defaults.jl")
include("discretization_state.jl")
include("parse_boundaries.jl")
include("periodic_map.jl")
include("make_pdesys_compatible.jl")
include("symbolic_discretize.jl")
include("StrongForm/StrongForm.jl")
include("dae_support.jl")
include("precompilation.jl")

export AbstractDiscreteSpace, AbstractCartesianDiscreteSpace, AbstractVarEqMapping,
    AbstractDifferentialDiscretizer, AbstractDiscretizationState,
    AbstractEquationSystemDiscretization, AbstractOptimizationSystemDiscretization,
    AbstractBoxDiscretization, AbstractMeshDiscretization

export AbstractBoundary, AbstractTruncatingBoundary, AbstractInterfaceBoundary,
    LowerBoundary, UpperBoundary, InterfaceBoundary,
    HigherOrderInterfaceBoundary, PieceBoundary

export BoundaryPiece, boundarypiece, boundarypieces, NormalDerivative

export get_time
export count_differentials, differential_order, has_derivatives, find_derivative, d_orders,
    subs_alleqs!, get_depvars, getvars, pde_substitute,
    get_all_depvars, split_terms, split_additive_terms, subsmatch, ex2term, safe_unwrap,
    recursive_unwrap, flatten_vardict, filter_interfaces, isperiodic, isinterface,
    haslowerupper, has_interfaces, isupper
export VariableMap, spatial_domain
export MeshDomain, Region, regions
export ivs, all_ivs, depvar, depvars, indvars, replaced_vars, x2i
export PeriodicMap
export classify_argument, classify_call, classify_condition, substitute_traces, assert_no_traces
export to_box_bcs

# Argument classes, locations and condition classes returned by the classifier. Their
# names are generic, so they stay qualified.
@public ArgumentClass, Free, BoxFace, Piece, Interior, ConditionClass, location, boxface, normal

# Mesh-domain interface, implemented by mesh backends.
@public meshdim, facetset_names, cellset_names, contains_point, backend_mesh

# Syntactic boundary-condition classification.
@public BCKind, Essential, Natural, Robin, bckind

# Versioned extension points for PDE discretization packages. These remain
# qualified so that ordinary applications do not accidentally depend on them.
@public lower_boundaries, interface_errors, check_boundarymap, should_transform,
    transform_pde_system!, split_ics_bcs,
    construct_disc_state, construct_discrete_space, construct_var_equation_mapping,
    construct_differential_discretizer, discretize_equation!, generate_ic_defaults,
    generate_metadata, generate_system, get_discvars, get_system_unknowns, get_eqvar,
    add_metadata!, make_pdesys_compatible, StrongForm, supported_conditions, SupportedConditions,
    BrownFullBasicInitUnsafeError, brown_init_offenders, differential_unknowns,
    dae_operating_point, time_derivative_order

end # module PDEBase
