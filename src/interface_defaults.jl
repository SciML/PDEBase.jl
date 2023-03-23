############################################################################################
# Default interface functions for `symbolic_discretize`
############################################################################################

"""
    interface_errors(sys::PDESystem, v::VariableMap, disc::AbstractDiscretization)

Use this function to check that your discretization type can discretize the given PDESystem.
"""
interface_errors(sys::PDESystem, v::VariableMap, disc::AbstractDiscretization) = nothing

"""
    check_boundarymap(bmap, v::VariableMap, disc::AbstractDiscretization)

    Once the boundaries have been parsed, check that they are valid for the given discretization.
"""
check_boundarymap(bmap, v::VariableMap, disc::AbstractDiscretization) = nothing

"""
    should_transform(pdesys::PDESystem, disc::AbstractDiscretization)

Boolean function that determines whether a PDESystem should be transformed to
make it compatible with the given discretization.
"""
should_transform(pdesys::PDESystem, disc::AbstractDiscretization) = false

"""
    transform_pde_system!(v::VariableMap, boundarymap, pdesys::PDESystem, disc::AbstractDiscretization)

Transforms the PDESystem to make it compatible with the given discretization.
"""
function transform_pde_system!(v::VariableMap, boundarymap,
                               pdesys::PDESystem, disc::AbstractDiscretization)
    nothing
end

"""
    construct_disc_state(disc::AbstractDiscretization)

Constructs the container which will hold all information needed to construct
the discretized PDESystem. Must be mutable, or contain mutables.
"""
construct_disc_state(::AbstractDiscretization) = []

construct_disc_state(::AbstractEquationSystemDiscretization) = EquationState()

"""
    construct_discrete_space(v::VariableMap, disc::AbstractDiscretization)

"""
construct_discrete_space(v::VariableMap, disc::AbstractDiscretization) = nothing

"""
    construct_var_equation_mapping(pdeeqs, bmap, s::AbstractDiscreteSpace, disc::AbstractDiscretization)

Constructs the mapping from the variables in the PDESystem to which equations they are solved from,
given the equations, boundary conditions, and discretization. This is a good time to calculate any
extra information needed for the boundary condition handling.
"""
function construct_var_equation_mapping(pdeeqs, bmap,
                                        s::AbstractDiscreteSpace,
                                        disc::AbstractDiscretization)
    nothing
end

"""
    construct_differential_discretizer(pdesys::PDESystem, s::AbstractDiscreteSpace,
                                       disc::AbstractDiscretization, orders)

Constructs the differential discretizer for the given discretization, given the orders of
derivative for each independent variable that are present in the PDESystem.

This is a type useful for doing any precalculations that are needed for the discretization.
"""
construct_differential_discretizer(pdesys, s, discretization, orders) = nothing

"""
    discretize_equation!(disc_state::AbstractDiscretizationState, pde::Equaation, vareqmap::AbstractVarEqMap,
                         eqvar, bcmap, depvars, s::AbstractDiscreteSpace,
                         derivweights::AbstractDifferentialDiscretizer, indexmap, discretization::AbstractDiscretization)
Add the information from the given pde equation to the discretization state. Also, discretize the bcs associated with the eqvar.

eqvar: the variable that this equation is for, contains all free ivs in the equation.
bcmap: boundarymap with any conditions on time removed.
depvars: the dependent variables in the equation.
indexmap: dict mapping each iv in this equation to its index in the eqvar.
"""
function discretize_equation!(disc_state::AbstractDiscretizationState, pde::Equation,
                              vareqmap::AbstractVarEqMapping,
                              eqvar, bcmap, depvars, s::AbstractDiscreteSpace,
                              derivweights::AbstractDifferentialDiscretizer, indexmap,
                              discretization::AbstractDiscretization)
    nothing
end

"""
    generate_ic_defaults(ics, s::AbstractDiscreteSpace disc::AbstractDiscretization)

Generate the default initial conditions for the given discretization. This is a good time to check
that the initial conditions are valid for the given discretization.
ics is a list of all conditions on the time variable.
"""
generate_ic_defaults(ics, s::AbstractDiscreteSpace, disc::AbstractDiscretization) = nothing

"""
    generate_metadata(s, discretization, pdesys, boundarymap)

Generate the metadata for the discretization. This can be used to store any extra information that
may be needed to convert the system in to a problem, or to reshape the solution once it is solved.
It is a good idea for this to return a subtype of SciMLBase.AbstractDiscretizationMetadata.
"""
generate_metadata(s, discretization, pdesys, boundarymap) = nothing

"""
    generate_system(disc_state::AbstractDiscretizationState, s::AbstractDiscreteSpace, u0, tspan, metadata::AbstractDiscretizationMetadata, discretization::AbstractDiscretization)

Generate the discretized system from the discretization state, the discretization, and the metadata.
You likely want this to return one of the symbolic system types from ModelingToolkit.jl.
u0: the return value of generate_ic_defaults
tspan: the time span of the problem, if relevant. (else nothing)
"""
function generate_system(disc_state::AbstractDiscretizationState, s::AbstractDiscreteSpace,
                         u0, tspan, metadata::AbstractDiscretizationMetadata,
                         discretization::AbstractDiscretization)
    nothing
end

############################################################################################
# Default interface functions for `AbstractDiscretization`
############################################################################################
"""
    get_time(discretization::AbstractDiscretization)

Get the time variable for the given discretization.
"""
get_time(discretization::AbstractDiscretization) = nothing

############################################################################################
# Default interface functions for `AbstractDiscreteSpace`
############################################################################################

"""
    get_discvars(s::AbstractDiscreteSpace)

Get the discrete variable map, where relevant
"""
get_discvars(s::AbstractDiscreteSpace) = []

############################################################################################
# Default interface functions for `AbstractVarEqMapping`
############################################################################################
"""
    get_eqvar(vareqmap::AbstractVarEqMapping)

Get the variable that the given equation will be solved for.
"""
get_eqvar(vareqmap::AbstractVarEqMapping, eq) = nothing

############################################################################################
# Default interface functions for `AbstractDiscretizationMetadata`
############################################################################################
"""
    add_metadata!(metadata::AbstractDiscretizationMetadata, value)
"""
function add_metadata!(metadata::AbstractDiscretizationMetadata, value)
    metadata.metadata[] = value
end

"""
    get_metadata(metadata::AbstractDiscretizationMetadata)
"""
get_metadata(metadata::AbstractDiscretizationMetadata) = metadata.metadata[]
