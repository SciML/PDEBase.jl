"""
    VariableMap

A container that maps dependent and independent variables from a `PDESystem` along with their properties.

`VariableMap` is created during `SciMLBase.symbolic_discretize` and is passed to
discretizer hooks. `VariableMap(pdesys)` is useful when a package needs to
inspect the symbolic variables without selecting a discretization.

# Fields
- `ū`: Collection of all dependent variables (filtered to exclude boundary values)
- `x̄`: Collection of spatial independent variables (excludes time)
- `ps`: Parameters of the system
- `time`: The time variable (if the system is time-dependent)
- `intervals`: Dictionary mapping each independent variable to its domain bounds (min, max)
- `args`: Dictionary mapping each dependent variable operation to its arguments
- `depvar_ops`: Operations (function symbols) of the dependent variables
- `x2i`: Dictionary mapping independent variables to their dimension indices
- `i2x`: Dictionary mapping dimension indices to independent variables
- `replaced_vars`: Dictionary of variable substitutions/replacements
- `domains`: Dictionary mapping each independent variable to the raw domain pairing that
  declares it. A tuple pairing `(x, y) ∈ Ω` is entered once for each of its variables.
- `spatial_axes`: The spatial independent variables in the order the domain list declares
  them, which is the axis order of the spatial box (see [`spatial_domain`](@ref)).

# Constructors
- `VariableMap(pdesys)`: Build a map without a discretization.
- `VariableMap(pdesys, disc; replaced_vars = Dict())`: Build a map for `disc`.

# Returns
- `VariableMap`: The normalized variable and domain information used by PDEBase.

# Examples
```julia
v = VariableMap(pdesys)
depvars(v)
indvars(v)
```
"""
struct VariableMap
    ū::Any
    x̄::Any
    ps::Any
    time::Any
    intervals::Any
    args::Any
    depvar_ops::Any
    x2i::Any
    i2x::Any
    replaced_vars::Any
    domains::Any
    spatial_axes::Any
end

function VariableMap(ū, x̄, ps, time, intervals, args, depvar_ops, x2i, i2x, replaced_vars)
    return VariableMap(
        ū, x̄, ps, time, intervals, args, depvar_ops, x2i, i2x, replaced_vars,
        Dict{Any, Any}(), Any[]
    )
end

# The variables of `x̄` in the order the pairings of `domain` name them.
function _declared_order(x̄, domain)
    axes = Any[]
    for pairing in domain
        coords = _pairing_coordinates(pairing)
        for c in (coords === nothing ? [safe_unwrap(pairing.variables)] : coords)
            i = findfirst(x -> isequal(safe_unwrap(x), c), x̄)
            i === nothing || any(a -> isequal(a, x̄[i]), axes) || push!(axes, x̄[i])
        end
    end
    return axes
end

# The coordinates named by a tuple pairing `(x, y) ∈ Ω`, or `nothing` for a scalar pairing.
function _pairing_coordinates(pairing)
    vars = safe_unwrap(pairing.variables)
    vars isa Union{Tuple, AbstractVector} && return collect(map(safe_unwrap, vars))
    if iscall(vars) && operation(vars) === tuple
        return collect(arguments(vars))
    end
    return nothing
end

function _check_mesh_pairing(pairing, ncoords)
    d = pairing.domain
    d isa MeshDomain || return nothing
    n = meshdim(d.mesh)
    n == ncoords || throw(
        ArgumentError(
            "domain MeshDomain of a $(typeof(d.mesh)) mesh has $n coordinates but the pairing " *
                "names $ncoords; write one coordinate per mesh dimension, e.g. (x, y) ∈ Ω for a 2D mesh"
        )
    )
    return nothing
end

_is_interval_product(d) =
    d isa DomainSets.ProductDomain && all(c -> c isa DomainSets.AbstractInterval, DomainSets.components(d))

# The pairing declaring `x` and, for a tuple pairing, the position of `x` in the tuple.
function _find_pairing(x, domain)
    x = safe_unwrap(x)
    for pairing in domain
        coords = _pairing_coordinates(pairing)
        if coords === nothing
            isequal(x, safe_unwrap(pairing.variables)) && return pairing, nothing
        else
            i = findfirst(c -> isequal(x, c), coords)
            i === nothing || return pairing, i
        end
    end
    return nothing, nothing
end

function VariableMap(pdesys, disc; replaced_vars = Dict())
    time = safe_unwrap(get_time(disc))
    eqs = get_eqs(pdesys)
    bcs = get_bcs(pdesys)
    depvars = get_dvs(pdesys)
    domain = get_domain(pdesys)
    ps = get_ps(pdesys)

    if ps isa SciMLBase.NullParameters
        ps = []
    end
    depvar_ops = get_ops(depvars)
    # Flatten equations - in MTK v11, complex equations may be nested vectors (real/imag splits)
    flat_eqs = _flatten_eqs(eqs)
    flat_bcs = _flatten_bcs(bcs)
    # Get all dependent variables in the correct type
    # A trace `u(t, Γ)` has a domain in place of its spatial arguments; it declares nothing.
    alldepvars = filter(
        u -> !any(x -> unwrap_const(x) isa DomainSets.Domain, arguments(u)),
        get_all_depvars([flat_eqs; flat_bcs], depvar_ops)
    )
    # Filter out boundaries
    ū = filter(u -> !any(x -> unwrap_const(x) isa Number, arguments(u)), alldepvars)
    # Get all independent variables in the correct type
    # `alldepvars` comes out of a `Set`, whose order differs between Julia versions; the order of
    # `x̄` follows the order of `depvars`, a call with all its arguments symbolic before a boundary
    # call of the same variable, and the printed arguments last.
    ivorder = sort(
        alldepvars; by = u -> (
            something(findfirst(op -> isequal(op, operation(u)), depvar_ops), length(depvar_ops) + 1),
            count(x -> unwrap_const(x) isa Number, arguments(u)), string(arguments(u)),
        )
    )
    allivs = collect(filter(x -> !(unwrap_const(x) isa Number), reduce(union, map(arguments, ivorder))))
    x̄ = remove(allivs, time)
    for pairing in domain
        coords = _pairing_coordinates(pairing)
        _check_mesh_pairing(pairing, coords === nothing ? 1 : length(coords))
    end
    intervals = Dict()
    domains = Dict{Any, Any}()
    for x in allivs
        pairing, i = _find_pairing(x, domain)
        if pairing === nothing
            error("Could not find domain for variable $x. Available domains: $([d.variables for d in domain])")
        end
        domains[x] = pairing
        xdomain = pairing.domain
        if i === nothing
            xdomain isa MeshDomain && continue
            intervals[x] = (DomainSets.infimum(xdomain), DomainSets.supremum(xdomain))
        elseif _is_interval_product(xdomain)
            axis = DomainSets.components(xdomain)[i]
            intervals[x] = (DomainSets.infimum(axis), DomainSets.supremum(axis))
        end
    end
    nspace = length(x̄)
    args = [operation(u) => arguments(u) for u in ū]
    x̄2dim = [x̄[i] => i for i in 1:nspace]
    dim2x̄ = [i => x̄[i] for i in 1:nspace]
    return VariableMap(
        ū, x̄, ps, time, Dict(intervals), Dict(args),
        depvar_ops, Dict(x̄2dim), Dict(dim2x̄), replaced_vars, domains,
        _declared_order(x̄, domain)
    )
end

VariableMap(pdesys) = VariableMap(pdesys, nothing)

"""
    spatial_domain(v::VariableMap)

Return the domain of the spatial coordinates `indvars(v)`, or `nothing` when there are none.

When one pairing declares exactly the spatial coordinates, as `(x, y) ∈ Ω` or
`x ∈ Interval(0, 1)`, its domain is returned as given. Otherwise, when every spatial
coordinate has an interval (its own interval pairing, or its axis of a box pairing that also
names time, as in `(t, x) ∈ Rectangle(...)`), the result is the interval of a single
coordinate or the `DomainSets.ProductDomain` of the intervals.

The axes of that box follow the user's declaration, never the order of `indvars(v)`, which
comes from the argument order of the dependent variables: for one tuple pairing they are its
coordinates in tuple order with time dropped, and for split pairings
`[t ∈ It, x ∈ Ix, y ∈ Iy]` they are the spatial pairings in the order of the domain list, so
the box is `Rectangle(Ix, Iy)`.
Throws an `ArgumentError` for a `MeshDomain` pairing that also names a non-spatial
coordinate, and when the spatial coordinates are declared by more than one pairing and at
least one of them is not an interval or interval product, for example a `MeshDomain` for
`(x, y)` and an interval for `z`.
"""
function spatial_domain(v::VariableMap)
    isempty(v.x̄) && return nothing
    axes = _box_axes(v)
    isempty(axes) && throw(
        ArgumentError("this VariableMap was built without its domain pairings, so it has no spatial domain")
    )
    pairings = unique(x -> objectid(v.domains[x]), axes)
    if length(pairings) == 1
        pairing = v.domains[only(pairings)]
        coords = _pairing_coordinates(pairing)
        (coords === nothing || length(coords) == length(axes)) && return pairing.domain
        pairing.domain isa MeshDomain && throw(
            ArgumentError(
                "the MeshDomain pairing $(pairing) also names coordinates that are not spatial " *
                    "variables of the system (such as time); pair only the spatial coordinates with a MeshDomain"
            )
        )
    end
    if all(x -> haskey(v.intervals, x), axes)
        intervals = map(x -> _axis_domain(v.domains[x], x), axes)
        return length(intervals) == 1 ? only(intervals) : DomainSets.ProductDomain(intervals...)
    end
    throw(
        ArgumentError(
            "the spatial variables $(v.x̄) are declared by more than one domain pairing " *
                "($([v.domains[x] for x in pairings])) and not all of them are intervals; " *
                "declare the spatial coordinates with one pairing, e.g. (x, y) ∈ Ω"
        )
    )
end

# The spatial variables in box axis order: axis `i` of `spatial_domain(v)` is `_box_axes(v)[i]`.
_box_axes(v::VariableMap) = v.spatial_axes

function _axis_domain(pairing, x)
    coords = _pairing_coordinates(pairing)
    coords === nothing && return pairing.domain
    return DomainSets.components(pairing.domain)[findfirst(c -> isequal(safe_unwrap(x), c), coords)]
end

function update_varmap!(v, newdv)
    push!(v.ū, newdv)
    merge!(v.args, Dict(operation(newdv) => arguments(newdv)))
    return push!(v.depvar_ops, operation(safe_unwrap(newdv)))
end

"""
    ivs(u, v::VariableMap)

Return the independent variables of dependent variable `u`, excluding the time
variable stored in `v`.
"""
ivs(u, v::VariableMap) = remove(v.args[operation(u)], v.time)

Base.ndims(u, v::VariableMap) = length(ivs(u, v))

"""
    all_ivs(v::VariableMap)
    all_ivs(u, v::VariableMap)

Return all independent variables tracked by `v`, including time when present.
With a dependent variable `u`, return the full argument list for `u`.
"""
all_ivs(v::VariableMap) = v.time === nothing ? v.x̄ : v.x̄ ∪ [v.time]

all_ivs(u, v::VariableMap) = v.args[operation(u)]

"""
    depvar(u, v::VariableMap)

Reconstruct dependent variable `u` with the argument ordering stored in `v`.
"""
depvar(u, v::VariableMap) = operation(u)(v.args[operation(u)]...)

"""
    depvars(v::VariableMap)

Return the dependent variables tracked by `v`.
"""
depvars(v::VariableMap) = v.ū

"""
    replaced_vars(v::VariableMap)

Return the variable replacements made while normalizing array-valued dependent
variables. The dictionary maps each flattened scalar variable to its original
indexed variable, allowing discretizers to translate solution queries between
the normalized and user-facing forms. The map is empty when no array-valued
dependent variables are present.
"""
replaced_vars(v::VariableMap) = v.replaced_vars

"""
    indvars(v::VariableMap)

Return the spatial independent variables tracked by `v`, excluding time.
"""
indvars(v::VariableMap) = v.x̄

"""
    x2i(v::VariableMap, u, x)

Return the spatial dimension index of independent variable `x` in the arguments
of dependent variable `u`.
"""
x2i(v::VariableMap, u, x) = findfirst(isequal(x), remove(v.args[operation(u)], v.time))

@inline function axiesvals(v::VariableMap, u_, x_, I)
    u = depvar(u_, v)
    return map(ivs(u, v)) do x
        x => (I[x2i(v, u, x)] == 1 ? v.intervals[x][1] : v.intervals[x][2])
    end
end

# Helper functions to flatten nested equation structures
# In MTK v11, complex equations may result in nested Vector{Equation}
function _flatten_eqs(eqs)
    result = Equation[]
    for eq in eqs
        if eq isa Equation
            push!(result, eq)
        elseif eq isa AbstractVector
            append!(result, _flatten_eqs(eq))
        end
    end
    return result
end

function _flatten_bcs(bcs)
    result = Any[]
    for bc in bcs
        if bc isa Equation || bc isa Pair
            push!(result, bc)
        elseif bc isa AbstractVector
            append!(result, _flatten_bcs(bc))
        end
    end
    return result
end
