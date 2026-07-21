module PDEBase
using SciMLBase
using SciMLBase: AbstractDiscretization, AbstractDiscretizationMetadata

using ModelingToolkit

using ModelingToolkit: get_unknowns, ProblemTypeCtx, get_ps, get_bcs, get_dvs,
    get_eqs, get_iv,
    get_domain, get_ivs, get_systems, get_connector_type,
    get_metadata, get_gui_metadata

using Symbolics, SymbolicUtils
using Symbolics: unwrap, setname, rename
using SymbolicUtils: operation, arguments, Chain, Prewalk, maketerm, metadata,
    iscall, getmetadata, unwrap_const, substitute
using SymbolicIndexingInterface: is_time_dependent
using DomainSets

"""
    AbstractEquationSystemDiscretization <: AbstractDiscretization

Supertype for discretizations that lower a `PDESystem` into an equation-based
SciML system, typically an `ODESystem` or `DAESystem`.
"""
abstract type AbstractEquationSystemDiscretization <: AbstractDiscretization end

"""
    AbstractOptimizationSystemDiscretization <: AbstractDiscretization

Supertype for discretizations that lower a `PDESystem` into an optimization
system.
"""
abstract type AbstractOptimizationSystemDiscretization <: AbstractDiscretization end

"""
    AbstractDiscreteSpace

Supertype for discretized spatial domain representations produced from a
`PDESystem` and a discretization.
"""
abstract type AbstractDiscreteSpace end

"""
    AbstractCartesianDiscreteSpace <: AbstractDiscreteSpace

Supertype for discrete space representations whose coordinates form a
Cartesian product grid.
"""
abstract type AbstractCartesianDiscreteSpace <: AbstractDiscreteSpace end

"""
    AbstractVarEqMapping

Supertype for mappings that assign dependent variables and boundary conditions
to the equations used to discretize them.
"""
abstract type AbstractVarEqMapping end

"""
    AbstractDifferentialDiscretizer

Supertype for objects that hold the data needed to discretize differential
operators on a discrete space.
"""
abstract type AbstractDifferentialDiscretizer end

"""
    AbstractDiscretizationState

Supertype for mutable state accumulated while PDE equations and boundary
conditions are discretized.
"""
abstract type AbstractDiscretizationState end

include("symbolic_utils.jl")
include("utils.jl")
include("variable_map.jl")
include("interface_defaults.jl")
include("discretization_state.jl")
include("parse_boundaries.jl")
include("periodic_map.jl")
include("make_pdesys_compatible.jl")
include("symbolic_discretize.jl")
include("domains.jl")
include("variational.jl")
include("precompilation.jl")

export AbstractDiscreteSpace, AbstractCartesianDiscreteSpace, AbstractVarEqMapping,
    AbstractDifferentialDiscretizer, AbstractDiscretizationState,
    AbstractEquationSystemDiscretization, AbstractOptimizationSystemDiscretization

export AbstractBoundary, AbstractTruncatingBoundary, AbstractInterfaceBoundary,
    LowerBoundary, UpperBoundary, InterfaceBoundary,
    HigherOrderInterfaceBoundary

export get_time
export count_differentials, differential_order, has_derivatives, find_derivative, d_orders,
    subs_alleqs!, get_depvars, getvars, pde_substitute,
    get_all_depvars, split_terms, split_additive_terms, subsmatch, ex2term, safe_unwrap,
    recursive_unwrap, flatten_vardict, filter_interfaces, isperiodic, isinterface,
    haslowerupper, has_interfaces, isupper
export VariableMap
export ivs, all_ivs, depvar, depvars, indvars, x2i
export PeriodicMap

# Domain types for unstructured meshes
export AbstractDiscretizedDomain, MeshDomain
export BoundaryRegion, ConditionalBoundary
export get_mesh, get_boundary_markers, get_boundary_name, get_parent_domain, get_condition
export mark_boundary!

# Variational method types for FEM/FVM/DG discretizations
export AbstractVariationalMethod, AbstractGalerkinMethod
export AbstractRitzGalerkin, AbstractPetrovGalerkin
export RitzGalerkin
export AbstractUnstructuredDiscretization, AbstractUnstructuredDiscreteSpace
export get_variational_method, is_galerkin, is_symmetric_galerkin
export AbstractFEMTrait, SupportsWeakForm, SupportsStrongForm, SupportsBothForms
export get_supported_form

end # module PDEBase
