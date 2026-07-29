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
    get_systems, get_unknowns, mtkcompile
using PrecompileTools: @compile_workload, @setup_workload
using SciMLBase: AbstractDiscretization, AbstractDiscretizationMetadata, NonlinearProblem,
    ODEFunction, ODEProblem
import SciMLPublic: @public
using Symbolics: @variables, Differential, Equation, Num, unwrap
using SymbolicUtils: arguments, getmetadata, iscall, operation, substitute, unwrap_const
using SymbolicIndexingInterface: is_time_dependent
using TermInterface: maketerm, metadata

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

# Versioned extension points for PDE discretization packages. These remain
# qualified so that ordinary applications do not accidentally depend on them.
@public interface_errors, check_boundarymap, should_transform, transform_pde_system!,
    construct_disc_state, construct_discrete_space, construct_var_equation_mapping,
    construct_differential_discretizer, discretize_equation!, generate_ic_defaults,
    generate_metadata, generate_system, get_discvars, get_eqvar, add_metadata!

end # module PDEBase
