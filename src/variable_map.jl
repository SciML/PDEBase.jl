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
    alldepvars = get_all_depvars([flat_eqs; flat_bcs], depvar_ops)
    # Filter out boundaries
    ū = filter(u -> !any(x -> unwrap_const(x) isa Number, arguments(u)), alldepvars)
    # Get all independent variables in the correct type
    allivs = collect(filter(x -> !(unwrap_const(x) isa Number), reduce(union, map(arguments, alldepvars))))
    x̄ = remove(allivs, time)
    intervals = Dict(
        map(allivs) do x
            idx = findfirst(d -> isequal(safe_unwrap(x), safe_unwrap(d.variables)), domain)
            if idx === nothing
                error("Could not find domain for variable $x. Available domains: $([d.variables for d in domain])")
            end
            xdomain = domain[idx]
            x => (DomainSets.infimum(xdomain.domain), DomainSets.supremum(xdomain.domain))
        end
    )
    nspace = length(x̄)
    args = [operation(u) => arguments(u) for u in ū]
    x̄2dim = [x̄[i] => i for i in 1:nspace]
    dim2x̄ = [i => x̄[i] for i in 1:nspace]
    return VariableMap(
        ū, x̄, ps, time, Dict(intervals), Dict(args),
        depvar_ops, Dict(x̄2dim), Dict(dim2x̄), replaced_vars
    )
end

VariableMap(pdesys) = VariableMap(pdesys, nothing)

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
