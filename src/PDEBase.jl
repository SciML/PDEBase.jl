module PDEBase
using SciMLBase
using SciMLBase: AbstractDiscretization, AbstractDiscretizationMetadata
using DomainSets

abstract type AbstractEquationSystemDiscretization <: AbstractDiscretization end
abstract type AbstractOptimizationSystemDiscretization <: AbstractDiscretization end

abstract type AbstractDiscreteSpace end
abstract type AbstractCartesianDiscreteSpace <: AbstractDiscreteSpace end
abstract type AbstractVarEqMapping end
abstract type AbstractDifferentialDiscretizer end
abstract type AbstractDiscretizationState end

include("symbolic_utils.jl")
include("variable_map.jl")
include("interface_defaults.jl")
include("discretization_state.jl")
include("parse_boundaries.jl")
include("symbolic_discretize.jl")

export AbstractDiscreteSpace, AbstractCartesianDiscreteSpace, AbstractVarEqMapping,
    AbstractDifferentialDiscretizer, AbstractDiscretizationState,
    AbstractEquattionSystemDiscretization, AbstractOptimizationSystemDiscretization

export get_time
export count_differentials, differential_order, has_derivatives, find_derivative, subs_alleqs!, get_depvars,
    get_all_depvars, split_terms, split_additive_terms, subsmatch, ex2term, safe_unwrap, recursive_unwrap

end # module PDEBase
