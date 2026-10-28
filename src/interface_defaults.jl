############################################################################################
# Default interface functions for `symbolic_discretize`
############################################################################################

"""
    lower_boundaries(pdesys::PDESystem, disc::AbstractDiscretization)

Rewrite the boundary conditions of `pdesys` into the form `disc` consumes.

This is the first hook `symbolic_discretize` calls, before complex splitting,
equation normalization and `VariableMap` construction, so every later hook sees
the lowered system. The default returns `pdesys` unchanged. Extend it for a
discretization family whose downstream code expects one spelling of boundary
conditions, such as endpoint-substituted calls on a box.

# Arguments
- `pdesys::PDESystem`: The system supplied to `symbolic_discretize`.
- `disc::AbstractDiscretization`: The discretization being applied.

# Returns
- `PDESystem`: The system used by the rest of the pipeline; `pdesys` itself by default.

# Examples
```julia
function PDEBase.lower_boundaries(pdesys::PDESystem, disc::MyDiscretization)
    return rewrite_boundaries_for_my_backend(pdesys)
end
```
"""
lower_boundaries(pdesys::PDESystem, disc::AbstractDiscretization) = pdesys

"""
    interface_errors(sys::PDESystem, v::VariableMap, disc::AbstractDiscretization)

Check whether `disc` supports the structure of `sys` represented by `v`.

Discretization packages should extend this function with checks for features that
cannot be handled by their implementation. Return `nothing` when the system is
supported; otherwise throw an informative exception before any boundary parsing
or discretization work begins.

# Arguments
- `sys::PDESystem`: The PDE system being discretized.
- `v::VariableMap`: The variable map constructed from `sys`.
- `disc::AbstractDiscretization`: The discretization being applied.

# Returns
- `nothing`: The system is supported.

# Examples
```julia
function PDEBase.interface_errors(sys, v, disc::MyDiscretization)
    isempty(v.x̄) && throw(ArgumentError("MyDiscretization requires a spatial variable"))
    return nothing
end
```
"""
interface_errors(sys::PDESystem, v::VariableMap, disc::AbstractDiscretization) = nothing

"""
    check_boundarymap(bmap, v::VariableMap, disc::AbstractDiscretization)

Validate the parsed boundary map before discretizing equations.

Extend this function when a discretizer imposes requirements beyond the
`PDESystem` boundary syntax, such as requiring both sides of every spatial
dimension. The map is the nested representation produced by PDEBase, with
dependent-variable operations at the outer level and independent variables at
the inner level.

# Arguments
- `bmap`: The parsed boundary map.
- `v::VariableMap`: The variable map for the PDE system.
- `disc::AbstractDiscretization`: The discretization being applied.

# Returns
- `nothing`: The boundary map is valid.

# Examples
```julia
function PDEBase.check_boundarymap(bmap, v, disc::MyDiscretization)
    isempty(PDEBase.flatten_vardict(bmap)) &&
        throw(ArgumentError("MyDiscretization requires boundary conditions"))
    return nothing
end
```
"""
check_boundarymap(bmap, v::VariableMap, disc::AbstractDiscretization) = nothing

"""
    should_transform(pdesys::PDESystem, disc::AbstractDiscretization, boundarymap)

Report whether `pdesys` must be transformed before discretization.

This hook is called after boundary validation. A `true` result must be paired
with a `transform_pde_system!` method for the same discretizer; a `false` result
leaves the original system unchanged.

# Arguments
- `pdesys::PDESystem`: The validated PDE system.
- `disc::AbstractDiscretization`: The discretization being applied.
- `boundarymap`: The parsed boundary map.

# Returns
- `Bool`: `true` when the transformation hook should run, otherwise `false`.
"""
should_transform(pdesys::PDESystem, disc::AbstractDiscretization, boundarymap) = false

"""
    transform_pde_system!(v::VariableMap, boundarymap, pdesys::PDESystem, disc::AbstractDiscretization)

Transform a PDE system into the representation required by `disc`.

This hook is called only when `should_transform` returns `true`. Update or
replace the system consistently with `v` and `boundarymap`; the returned system
is used for all subsequent interface calls. A transformation must preserve the
meaning of the equations and boundary conditions.

# Arguments
- `v::VariableMap`: The variable map built before transformation.
- `boundarymap`: The parsed boundary map associated with `pdesys`.
- `pdesys::PDESystem`: The system to transform.
- `disc::AbstractDiscretization`: The discretization that will consume it.

# Returns
- `PDESystem`: The transformed system passed to the construction hooks.

# Examples
```julia
function PDEBase.transform_pde_system!(v, boundarymap, pdesys, disc::MyDiscretization)
    return rebuild_for_my_backend(pdesys)
end
```
"""
function transform_pde_system!(
        v::VariableMap, boundarymap,
        pdesys::PDESystem, disc::AbstractDiscretization
    )
    return nothing
end

"""
    split_ics_bcs(boundarymap, v::VariableMap, disc::AbstractDiscretization)

Separate the parsed boundary map into initial conditions and spatial boundary
conditions.

The default collects `boundarymap[u][t]` over every dependent variable `u`, where
`t = get_time(disc)`, as the initial conditions (none when `t === nothing`), and
builds `bcmap[u][x] = boundarymap[u][x]` for every spatial variable `x` in
`indvars(v)`. For an `AbstractMeshDiscretization` it builds `bcmap[u][Γ] =
boundarymap[u][Γ]` for every [`BoundaryPiece`](@ref) key `Γ` instead. The initial conditions are passed to `generate_ic_defaults` and
`bcmap` to `discretize_equation!`. Extend it when `boundarymap` is keyed by
something other than the spatial independent variables.

# Arguments
- `boundarymap`: The parsed boundary map, after any transformation.
- `v::VariableMap`: The variable map for the PDE system.
- `disc::AbstractDiscretization`: The discretization being applied.

# Returns
- `(ics, bcmap)`: A vector of initial-condition boundaries and a nested `Dict`
  of spatial boundaries keyed by dependent-variable operation and then by key.

# Examples
```julia
function PDEBase.split_ics_bcs(boundarymap, v::VariableMap, disc::MyDiscretization)
    t = PDEBase.get_time(disc)
    ics = mapreduce(u -> boundarymap[u][t], vcat, keys(boundarymap))
    bcmap = Dict(u => filter(kv -> !isequal(first(kv), t), boundarymap[u])
        for u in keys(boundarymap))
    return ics, bcmap
end
```
"""
function split_ics_bcs(boundarymap, v::VariableMap, disc::AbstractDiscretization)
    t = get_time(disc)
    ics = t === nothing ? [] :
        mapreduce(u -> boundarymap[u][t], vcat, operation.(depvars(v)))

    bcmap = Dict(
        map(collect(keys(boundarymap))) do u
            u => Dict(
                map(indvars(v)) do x
                    x => boundarymap[u][x]
                end
            )
        end
    )
    return ics, bcmap
end

# A mesh boundary map is keyed by boundary pieces instead of spatial variables.
function split_ics_bcs(boundarymap, v::VariableMap, disc::AbstractMeshDiscretization)
    t = get_time(disc)
    ics = t === nothing ? [] :
        mapreduce(u -> boundarymap[u][t], vcat, operation.(depvars(v)))
    bcmap = Dict(
        map(collect(keys(boundarymap))) do u
            u => Dict(filter(kv -> first(kv) isa BoundaryPiece, collect(boundarymap[u])))
        end
    )
    return ics, bcmap
end

"""
    construct_disc_state(disc::AbstractDiscretization)

Construct the mutable state that accumulates discretized equations.

The returned object is passed to `discretize_equation!` for every PDE and then
to `generate_system`. It must be mutable itself or contain mutable fields that
can hold the generated interior and boundary equations.

# Arguments
- `disc::AbstractDiscretization`: The discretization being applied.

# Returns
- `AbstractDiscretizationState`: State consumed by the discretization hooks.

# Examples
```julia
struct MyState <: PDEBase.AbstractDiscretizationState
    equations::Vector
end

PDEBase.construct_disc_state(::MyDiscretization) = MyState(Any[])
```
"""
construct_disc_state(::AbstractDiscretization) = []

construct_disc_state(::AbstractEquationSystemDiscretization) = EquationState()

"""
    construct_discrete_space(v::VariableMap, disc::AbstractDiscretization)

Construct the discrete spatial representation used by a discretizer.

# Arguments
- `v::VariableMap`: The variable map containing the dependent and independent variables.
- `disc::AbstractDiscretization`: The discretization scheme.

# Returns
- `AbstractDiscreteSpace`: A grid or mesh representation and any associated
  discrete variables.

# Examples
```julia
struct MySpace <: PDEBase.AbstractCartesianDiscreteSpace
    grid::Vector{Float64}
end

PDEBase.construct_discrete_space(v, ::MyDiscretization) = MySpace(collect(0:0.1:1))
```
"""
construct_discrete_space(v::VariableMap, disc::AbstractDiscretization) = nothing

"""
    construct_discrete_space(v::VariableMap, pdesys::PDESystem, disc::AbstractOptimizationSystemDiscretization)

Create the space on which an optimization-system discretization evaluates its residuals.
Unlike the grid of an equation-system discretization it may depend on the boundary
conditions of `pdesys` (for example to keep interior collocation points off the
boundaries), so the system is passed along with the variable map. The default forwards to
`construct_discrete_space(v, disc)`.

# Arguments
- `v::VariableMap`: The variable map constructed from `pdesys`.
- `pdesys::PDESystem`: The normalized PDE system being discretized.
- `disc::AbstractOptimizationSystemDiscretization`: The discretization being applied.

# Returns
- The discrete space passed as `s` to the remaining hooks.
"""
function construct_discrete_space(
        v::VariableMap, pdesys::PDESystem, disc::AbstractOptimizationSystemDiscretization
    )
    return construct_discrete_space(v, disc)
end

"""
    construct_var_equation_mapping(pdeeqs, bmap, s::AbstractDiscreteSpace, disc::AbstractDiscretization)

Construct the mapping from each PDE to the discrete variable it determines.

The mapping is also the place to precompute information needed when boundary
conditions are applied. `get_eqvar` must be able to retrieve the variable for
each equation in `pdeeqs` from the returned object.

# Arguments
- `pdeeqs`: The PDE equations after any optional transformation.
- `bmap`: The parsed boundary map.
- `s::AbstractDiscreteSpace`: The discrete space returned by
  `construct_discrete_space`.
- `disc::AbstractDiscretization`: The discretization scheme.

# Returns
- `AbstractVarEqMapping`: A mapping consumed by `get_eqvar` and
  `discretize_equation!`.
"""
function construct_var_equation_mapping(
        pdeeqs, bmap,
        s::AbstractDiscreteSpace,
        disc::AbstractDiscretization
    )
    return nothing
end

"""
    construct_differential_discretizer(pdesys::PDESystem, s::AbstractDiscreteSpace,
                                       disc::AbstractDiscretization, orders)

Construct the derivative approximation data for `disc`.

`orders` maps each independent variable to the derivative orders appearing in
the PDE and boundary conditions. This hook is where a discretizer should do
precomputation such as finite-difference weights or element-local operators.

# Arguments
- `pdesys::PDESystem`: The transformed PDE system.
- `s::AbstractDiscreteSpace`: The discrete space on which derivatives act.
- `disc::AbstractDiscretization`: The discretization scheme.
- `orders`: A map from independent variables to required derivative orders.

# Returns
- `AbstractDifferentialDiscretizer`: Derivative data consumed by
  `discretize_equation!`.
"""
construct_differential_discretizer(pdesys, s, discretization, orders) = nothing

"""
    discretize_equation!(disc_state::AbstractDiscretizationState, pde::Equation, vareqmap::AbstractVarEqMapping,
                         eqvar, bcmap, depvars, s::AbstractDiscreteSpace,
                         derivweights::AbstractDifferentialDiscretizer, indexmap, discretization::AbstractDiscretization)
Add one PDE and its associated spatial boundary conditions to `disc_state`.

This hook is called once for every PDE after the equation-to-variable mapping
has been built. It should append or mutate the discrete interior and boundary
equations in `disc_state`; it should not construct the final symbolic system.

# Arguments
- `disc_state::AbstractDiscretizationState`: Accumulator returned by
  `construct_disc_state`.
- `pde::Equation`: The PDE equation currently being discretized.
- `vareqmap::AbstractVarEqMapping`: Mapping returned by
  `construct_var_equation_mapping`.
- `eqvar`: The discrete variable solved for by `pde`.
- `bcmap`: Boundary map with conditions on the time variable removed.
- `depvars`: Dependent variables appearing in `pde`.
- `s::AbstractDiscreteSpace`: The discrete space.
- `derivweights::AbstractDifferentialDiscretizer`: Derivative data returned by
  `construct_differential_discretizer`.
- `indexmap`: Map from each independent variable in `eqvar` to its index.
- `discretization::AbstractDiscretization`: The active discretization.

# Returns
- `nothing`: The state is updated in place.
"""
function discretize_equation!(
        disc_state::AbstractDiscretizationState, pde::Equation,
        vareqmap::AbstractVarEqMapping,
        eqvar, bcmap, depvars, s::AbstractDiscreteSpace,
        derivweights::AbstractDifferentialDiscretizer, indexmap,
        discretization::AbstractDiscretization
    )
    return nothing
end

"""
    discretize_equation!(disc_state, eq::Equation, kind::Symbol, s, derivweights, disc::AbstractOptimizationSystemDiscretization)

Lower one equation of the `PDESystem` for an optimization-system discretization and
record the result in `disc_state`. `kind` is `:pde` for the equations and `:bc` for the
boundary conditions; each is lowered on its own residual domain, so boundary conditions
are not required to lie on a boundary of the domain.

# Arguments
- `disc_state`: The mutable state created by `construct_disc_state`.
- `eq::Equation`: The equation, cardinalized to `lhs - rhs ~ 0`.
- `kind::Symbol`: `:pde` or `:bc`.
- `s`: The discrete space created by `construct_discrete_space`.
- `derivweights`: The value returned by `construct_differential_discretizer`.
- `disc::AbstractOptimizationSystemDiscretization`: The discretization being applied.

# Returns
- `nothing`: The state is mutated in place.
"""
function discretize_equation!(
        disc_state, eq::Equation, kind::Symbol, s, derivweights,
        disc::AbstractOptimizationSystemDiscretization
    )
    return nothing
end

"""
    generate_ic_defaults(ics, s::AbstractDiscreteSpace, disc::AbstractDiscretization)

Generate discrete initial-condition values and validate the input conditions.

# Arguments
- `ics`: Conditions on the time variable extracted from the boundary map.
- `s::AbstractDiscreteSpace`: The discrete space.
- `disc::AbstractDiscretization`: The discretization.

# Returns
- `u0`: Initial values in the representation expected by `generate_system`.

# Examples
```julia
PDEBase.generate_ic_defaults(ics, space::MySpace, ::MyDiscretization) =
    [evaluate_initial_condition(ic, space) for ic in ics]
```
"""
generate_ic_defaults(ics, s::AbstractDiscreteSpace, disc::AbstractDiscretization) = nothing

"""
    generate_metadata(s, discretization, pdesys, boundarymap, complexmap, u0)

Create metadata needed to construct a problem or reconstruct a solution.

The metadata is passed unchanged to `generate_system`. Returning a subtype of
`SciMLBase.AbstractDiscretizationMetadata` allows the generic
`add_metadata!` hook to attach the original symbolic system.

# Arguments
- `s`: The discrete space.
- `discretization`: The active discretization.
- `pdesys`: The transformed PDE system.
- `boundarymap`: The parsed boundary map.
- `complexmap`: The map used while normalizing complex equations.
- `u0`: The discrete initial conditions returned by `generate_ic_defaults`.

# Returns
- `AbstractDiscretizationMetadata`: Metadata consumed by `generate_system`.
"""
generate_metadata(s, discretization, pdesys, boundarymap, complexmap, u0 = []) = nothing

"""
    generate_system(disc_state::AbstractDiscretizationState, s::AbstractDiscreteSpace, u0, tspan, metadata::AbstractDiscretizationMetadata, discretization::AbstractDiscretization)

Build the final symbolic system from the accumulated discretization results.

The returned value is passed back through `SciMLBase.symbolic_discretize`; in
the usual time-dependent case it is a tuple `(system, tspan)`. Set `tspan` on the
system as well, so that a caller who keeps only the system can still build a
time-dependent problem without supplying a time span. Use `checks` when
constructing the symbolic system and forward it to the constructor rather than
silently changing validation behavior.

# Arguments
- `disc_state::AbstractDiscretizationState`: State populated by
  `discretize_equation!`.
- `s::AbstractDiscreteSpace`: The discrete space.
- `u0`: Initial values returned by `generate_ic_defaults`.
- `tspan`: The time interval, or `nothing` for a steady-state system.
- `metadata::AbstractDiscretizationMetadata`: Metadata returned by
  `generate_metadata`.
- `discretization::AbstractDiscretization`: The active discretization.

# Keywords
- `checks::Bool = true`: Enable symbolic-system consistency checks.

# Returns
- `Tuple`: The generated symbolic system and its time span.

# Examples
```julia
function PDEBase.generate_system(state, space, u0, tspan, metadata,
        disc::MyDiscretization; checks = true)
    system = ODESystem(state.equations, get_time(disc); checks)
    return system, tspan
end
```
"""
function generate_system(
        disc_state::AbstractDiscretizationState, s::AbstractDiscreteSpace,
        u0, tspan, metadata::AbstractDiscretizationMetadata,
        discretization::AbstractDiscretization;
        checks = true
    )
    return nothing
end

############################################################################################
# Default interface functions for `AbstractDiscretization`
############################################################################################
"""
    get_time(discretization::AbstractDiscretization)

Return the independent time variable used by `discretization`.

Return `nothing` for steady-state discretizations. The result is used to
separate initial conditions from spatial boundary conditions and to determine
whether the generated system is time dependent.

# Arguments
- `discretization::AbstractDiscretization`: The discretization being applied.

# Returns
- A symbolic time variable, or `nothing` for a steady-state discretization.
"""
get_time(discretization::AbstractDiscretization) = nothing
get_time(_) = nothing

############################################################################################
# Default interface functions for `AbstractDiscreteSpace`
############################################################################################

"""
    get_discvars(s::AbstractDiscreteSpace)

Return the discrete variables associated with `s`.

The result is consumed by the final system constructor to identify the unknowns
created by the discretizer. Discretizations without a separate variable map
may return an empty collection.

# Arguments
- `s::AbstractDiscreteSpace`: The discrete space.

# Returns
- A collection or map of symbolic discrete variables.
"""
get_discvars(s::AbstractDiscreteSpace) = []

"""
    get_system_unknowns(s::AbstractDiscreteSpace)

Return the unknowns of the system `generate_system` builds from `s`, as a vector of
symbolic variables.

The default flattens [`get_discvars`](@ref) into its scalar entries, one unknown per
discrete variable. A discretizer that represents each dependent variable as one
array-valued symbolic, such as `u(t)[1:n]`, can override this to return those arrays
so that the system has one unknown per dependent variable. The scalar entries of
`get_discvars` then index into these arrays and remain what the discretized
equations, boundary conditions and initial conditions are written in terms of.

# Arguments
- `s::AbstractDiscreteSpace`: The discrete space.

# Returns
- A vector of symbolic unknowns.
"""
function get_system_unknowns(s::AbstractDiscreteSpace)
    discvars = get_discvars(s)
    return vec(reduce(vcat, vec(unique(reduce(vcat, vec.(values(discvars)))))))
end

############################################################################################
# Default interface functions for `AbstractVarEqMapping`
############################################################################################
"""
    get_eqvar(vareqmap::AbstractVarEqMapping, eq)

Return the discrete variable solved for by equation `eq`.

Every equation passed to `discretize_equation!` must have a corresponding
variable. The returned variable must use the same independent-variable ordering
that `ivs` and `x2i` expect for the associated `VariableMap`.

# Arguments
- `vareqmap::AbstractVarEqMapping`: Mapping returned by
  `construct_var_equation_mapping`.
- `eq`: The PDE equation being visited.

# Returns
- The discrete dependent variable associated with `eq`.
"""
get_eqvar(vareqmap::AbstractVarEqMapping, eq) = nothing

############################################################################################
# Default interface functions for `AbstractDiscretizationMetadata`
############################################################################################
"""
    add_metadata!(metadata::AbstractDiscretizationMetadata, value)

Attach `value` to the mutable metadata slot in `metadata`.

This hook is used by PDEBase to retain the original symbolic system after a
compiled system is built. Metadata implementations should provide a mutable
`metadata` field, typically a `Ref`, when they use the default method.

# Arguments
- `metadata::AbstractDiscretizationMetadata`: Metadata object to update.
- `value`: Value to store in its metadata slot.

# Returns
- The assigned `value`.
"""
function add_metadata!(metadata::AbstractDiscretizationMetadata, value)
    return metadata.metadata[] = value
end

############################################################################################
# Supported conditions: documentation and introspection only
############################################################################################
const CONDITION_KINDS = (
    :dirichlet => "Dirichlet", :neumann => "Neumann", :robin => "Robin",
    :periodic => "Periodic", :interface => "Interface", :time_dependent => "Time-dependent",
)

"""
    SupportedConditions(; geometries = (), dirichlet = nothing, neumann = nothing,
        robin = nothing, periodic = nothing, interface = nothing,
        time_dependent = nothing, notes = "")

Machine-readable description of the geometries and boundary-condition kinds a
discretization supports, returned by [`supported_conditions`](@ref PDEBase.supported_conditions).

Each condition field is `true` (supported), `false` (not supported) or `nothing`
(not declared). An empty `geometries` tuple likewise means not declared.

# Fields
- `geometries::Tuple`: Supported geometries, as symbols such as `:interval`,
  `:box` or `:mesh`, or as domain types.
- `dirichlet`, `neumann`, `robin`: Value, normal-derivative and mixed conditions.
- `periodic`: Periodic conditions linking opposite boundaries.
- `interface`: Conditions coupling variables across an interface between domains.
- `time_dependent`: Boundary data that varies with time.
- `notes::String`: Free-text restrictions, such as required derivative orders or
  unsupported combinations.

`show(io, MIME"text/markdown"(), sc)` renders the description as a Markdown table
for discretizer documentation.

# Examples
```julia
PDEBase.supported_conditions(::MyDiscretization) = PDEBase.SupportedConditions(
    geometries = (:box,), dirichlet = true, neumann = true, robin = true,
    periodic = true, interface = true
)
```
"""
struct SupportedConditions
    geometries::Tuple
    dirichlet::Union{Bool, Nothing}
    neumann::Union{Bool, Nothing}
    robin::Union{Bool, Nothing}
    periodic::Union{Bool, Nothing}
    interface::Union{Bool, Nothing}
    time_dependent::Union{Bool, Nothing}
    notes::String
end

function SupportedConditions(;
        geometries = (), dirichlet = nothing, neumann = nothing, robin = nothing,
        periodic = nothing, interface = nothing, time_dependent = nothing, notes = ""
    )
    return SupportedConditions(
        Tuple(geometries), dirichlet, neumann, robin, periodic, interface,
        time_dependent, notes
    )
end

function Base.isempty(sc::SupportedConditions)
    return isempty(sc.geometries) && isempty(sc.notes) &&
        all(kind -> getfield(sc, first(kind)) === nothing, CONDITION_KINDS)
end

support_label(::Nothing) = "unknown"
support_label(x::Bool) = x ? "yes" : "no"

function Base.show(io::IO, ::MIME"text/markdown", sc::SupportedConditions)
    geometries = isempty(sc.geometries) ? "unknown" : join(sc.geometries, ", ")
    println(io, "| Condition | Supported |")
    println(io, "|:--- |:--- |")
    println(io, "| Geometries | ", geometries, " |")
    for (field, label) in CONDITION_KINDS
        println(io, "| ", label, " | ", support_label(getfield(sc, field)), " |")
    end
    isempty(sc.notes) || print(io, "\n", sc.notes, "\n")
    return nothing
end

"""
    supported_conditions(disc::AbstractDiscretization) -> SupportedConditions

Declare the geometries and boundary-condition kinds supported by `disc`.

Discretizer packages extend this function so that their supported conditions are
documented in one machine-readable place. It is documentation and introspection
only: PDEBase never calls it while discretizing, so it does not replace
validation in [`interface_errors`](@ref PDEBase.interface_errors) or
[`check_boundarymap`](@ref PDEBase.check_boundarymap).

# Arguments
- `disc::AbstractDiscretization`: The discretization being described.

# Returns
- `SupportedConditions`: The default is `SupportedConditions()`, in which nothing
  is declared. This means unknown, not that nothing is supported.

# Examples
```julia
PDEBase.supported_conditions(::MyDiscretization) = PDEBase.SupportedConditions(
    geometries = (:box,), dirichlet = true, neumann = true, robin = true,
    periodic = true, interface = true, notes = "Uniform and nonuniform grids."
)
```
"""
supported_conditions(disc::AbstractDiscretization) = SupportedConditions()
