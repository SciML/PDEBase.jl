module PDEBase
using SciMLBase
using SciMLBase: AbstractDiscretization, AbstractDiscretizationMetadata

using ModelingToolkit

using ModelingToolkit: get_unknowns, ProblemTypeCtx,
                       parameters, defaults, varmap_to_vars, get_eqs, get_iv

using Symbolics, SymbolicUtils
using Symbolics: unwrap, solve_for, expand_derivatives, diff2term, setname, rename, variable
using SymbolicUtils: operation, arguments, Chain, Prewalk, Postwalk, maketerm, metadata,
      symtype, operation, iscall, arguments, getmetadata
using DomainSets

abstract type AbstractEquationSystemDiscretization <: AbstractDiscretization end
abstract type AbstractOptimizationSystemDiscretization <: AbstractDiscretization end

abstract type AbstractDiscreteSpace end
abstract type AbstractCartesianDiscreteSpace <: AbstractDiscreteSpace end

abstract type AbstractVarEqMapping end
abstract type AbstractDifferentialDiscretizer end
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

export AbstractDiscreteSpace, AbstractCartesianDiscreteSpace, AbstractVarEqMapping,
       AbstractDifferentialDiscretizer, AbstractDiscretizationState,
       AbstractEquationSystemDiscretization, AbstractOptimizationSystemDiscretization

export AbstractBoundary, AbstractTruncatingBoundary, AbstractInterfaceBoundary,
       LowerBoundary, UpperBoundary, InterfaceBoundary,
       HigherOrderInterfaceBoundary

export get_time
export count_differentials, differential_order, has_derivatives, find_derivative, d_orders,
       subs_alleqs!, get_depvars, getvars,
       get_all_depvars, split_terms, split_additive_terms, subsmatch, ex2term, safe_unwrap,
       recursive_unwrap, flatten_vardict, filter_interfaces, isperiodic, isinterface, haslowerupper, has_interfaces, isupper
export VariableMap
export ivs, all_ivs, depvar, depvars, indvars, x2i
export PeriodicMap

end # module PDEBase
